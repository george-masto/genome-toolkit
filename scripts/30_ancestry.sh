#!/usr/bin/env bash
# ============================================================================
# 30_ancestry.sh - genetic-ancestry composition by PCA projection.
#
# Estimates your super-population make-up (European, East Asian, African,
# South Asian, Middle Eastern, Admixed American, Oceanian) and wires it into the
# MCP as get_ancestry (analysis/calls/ancestry.json -> build_insights.py -> tool).
#
# WHY PROJECTION (you can't PCA one genome):
#   A single sample has nothing to be "principal" relative to. So we compute PCA
#   on a diverse REFERENCE panel (HGDP+1kGP - the same panel scripts/21_prs.sh
#   uses), then PROJECT your genome onto those fixed axes using the reference's
#   per-SNP allele weights, and see which reference samples you land nearest
#   (k-NN over the top PCs -> super-population vote shares).
#
# SCALE TRICK: we score the reference panel through the SAME --score step as the
#   subject, so both sets of PC coordinates share one scale and k-NN distances
#   are meaningful (raw .eigenvec and projected scores are scaled differently).
#
# Reuses 21_prs's reference/pgsc panel + reference FASTA, and the gVCF->dense
# reconstruction (fills hom-ref sites so common-SNP dosages are complete).
#
# NOTE: like the other heavy scripts (10_align, 21_prs, second-pass callers),
# this needs the reference panel + your genome and runs for a while; it is not
# exercised by the quick-start run_all.sh.
# ============================================================================
set -uo pipefail
source "$(dirname "$0")/lib/common.sh"

require_tool plink2
require_tool bcftools
require_vcf_snp

REF_PGSC="${REF_PGSC:-$REF_DIR/pgsc}"
REFP="${REFP:-$REF_PGSC/GRCh38_HGDP+1kGP_ALL}"     # .pgen / .pvar.zst / .psam
FA="${REF_FA:-$REF_DIR/GRCh38.primary.fa}"
WORK="$OUT_DIR/ancestry.work"
OUTJSON="$CALLS_DIR/ancestry.json"
NPC="${NPC:-10}"
KNN="${KNN:-25}"
mkdir -p "$WORK" "$CALLS_DIR"
LOG="$WORK/log"

# 0. reference panel present? (same one PRS downloads)
[[ -f "$REFP.pgen" && ( -f "$REFP.pvar.zst" || -f "$REFP.pvar" ) && -f "$REFP.psam" ]] || die \
"HGDP+1kGP reference panel not found at $REFP.*
   It's the same panel scripts/21_prs.sh uses. Fetch it once:
     mkdir -p '$REF_PGSC' && cd '$REF_PGSC'
     curl -fLO 'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/resources/pgsc_HGDP+1kGP_v1.tar.zst'
     tar --use-compress-program=unzstd -xf pgsc_HGDP+1kGP_v1.tar.zst"
[[ -f "$FA" ]] || die "Reference FASTA $FA not found (scripts/21_prs.sh downloads it; set REF_FA=)."

# 1. reference PCA on LD-pruned, common, biallelic autosomal SNPs (cached)
if [[ ! -f "$WORK/ref_pca.eigenvec.allele" || ! -f "$WORK/ref_pca.acount" ]]; then
  log "Reference PCA on HGDP+1kGP (LD-prune -> $NPC PCs with allele weights)…"
  plink2 --pfile "$REFP" vzs --autosome --snps-only --max-alleles 2 --maf 0.05 \
    --set-all-var-ids '@:#' --rm-dup exclude-all \
    --indep-pairwise 1000 50 0.2 --out "$WORK/prune" >>"$LOG" 2>&1 \
    || die "plink2 LD-prune failed (see $LOG)"
  plink2 --pfile "$REFP" vzs --autosome --set-all-var-ids '@:#' --rm-dup exclude-all \
    --extract "$WORK/prune.prune.in" \
    --freq counts --pca "$NPC" allele-wts \
    --out "$WORK/ref_pca" >>"$LOG" 2>&1 \
    || die "plink2 reference PCA failed (see $LOG)"
  ok "reference PCA built: $(wc -l < "$WORK/prune.prune.in" | tr -d ' ') SNPs, $NPC PCs"
fi

# 1b. project the REFERENCE samples through the same scoring (scale-matched coords)
if [[ ! -f "$WORK/ref_proj.sscore" ]]; then
  log "Projecting reference samples onto their own PCs (scale reference)…"
  plink2 --pfile "$REFP" vzs --set-all-var-ids '@:#' --rm-dup exclude-all \
    --read-freq "$WORK/ref_pca.acount" \
    --score "$WORK/ref_pca.eigenvec.allele" 2 5 header-read no-mean-imputation variance-standardize \
    --score-col-nums 6-$((5 + NPC)) \
    --out "$WORK/ref_proj" >>"$LOG" 2>&1 \
    || die "plink2 reference projection failed (see $LOG)"
fi

# 2. subject dense genotypes at the PCA SNP positions (hom-ref filled from FASTA)
awk -F: 'NF>=2{print $1"\t"$2}' "$WORK/prune.prune.in" | sort -k1,1V -k2,2n -u > "$WORK/pca_sites.txt"
log "Reconstructing dense genotypes at $(wc -l < "$WORK/pca_sites.txt" | tr -d ' ') PCA sites…"
bcftools view -R "$WORK/pca_sites.txt" "$VCF_SNP" -Ou 2>>"$LOG" \
  | bcftools convert --gvcf2vcf -f "$FA" -Oz -o "$WORK/subject.dense.vcf.gz" 2>>"$LOG" \
  || die "dense reconstruction failed (see $LOG)"
tabix -f -p vcf "$WORK/subject.dense.vcf.gz" 2>>"$LOG"

# 3. project the SUBJECT onto the same PCs (identical --score recipe)
plink2 --vcf "$WORK/subject.dense.vcf.gz" \
  --set-all-var-ids '@:#' --rm-dup exclude-all \
  --read-freq "$WORK/ref_pca.acount" \
  --score "$WORK/ref_pca.eigenvec.allele" 2 5 header-read no-mean-imputation variance-standardize \
  --score-col-nums 6-$((5 + NPC)) \
  --out "$WORK/subject_proj" >>"$LOG" 2>&1 \
  || die "plink2 subject projection failed (see $LOG)"

# 4. k-NN classify subject against reference in PC space -> composition JSON
log "k-NN ($KNN) super-population assignment…"
SAMPLE="$(sample_name)" REFPSAM="$REFP.psam" NPC="$NPC" KNN="$KNN" GEN="$(now_utc)" \
python3 - "$WORK/ref_proj.sscore" "$WORK/subject_proj.sscore" "$OUTJSON" <<'PY' 2>>"$LOG" || die "classification failed (see $LOG)"
import csv, json, math, os, sys
ref_f, sub_f, out_f = sys.argv[1], sys.argv[2], sys.argv[3]
npc, knn = int(os.environ["NPC"]), int(os.environ["KNN"])
SUPERPOP = {"AFR":"African","AMR":"Admixed American","CSA":"Central/South Asian",
            "EAS":"East Asian","EUR":"European","MID":"Middle Eastern","OCE":"Oceanian",
            "SAS":"South Asian"}

def score_cols(header):
    return [c for c in header if c.startswith("SCORE") and c.endswith("_AVG")][:npc]

def load(path):
    r = list(csv.DictReader(open(path), delimiter='\t'))
    if not r: sys.exit(f"empty score file: {path}")
    iid = next(c for c in r[0] if c.lstrip('#') == 'IID' or c == '#IID' or c.endswith('IID'))
    cols = score_cols(list(r[0].keys()))
    return r, iid, cols

ref, riid, rcols = load(ref_f)
sub, siid, scols = load(sub_f)
# super-population label per reference IID (psam col index 2, matching 21_prs)
pop = {}
for line in open(os.environ["REFPSAM"]):
    if line.startswith('#'): continue
    p = line.rstrip('\n').split('\t')
    if len(p) >= 3: pop[p[0]] = p[2]

s = sub[0]
svec = [float(s[c]) for c in scols]
dists = []
for row in ref:
    sp = pop.get(row[riid])
    if not sp: continue
    d = math.dist(svec, [float(row[c]) for c in rcols])
    dists.append((d, sp))
dists.sort(key=lambda t: t[0])
neigh = dists[:knn]
if not neigh: sys.exit("no reference neighbours with population labels")

tally = {}
for _, sp in neigh: tally[sp] = tally.get(sp, 0) + 1
comp = sorted(({"superpopulation": sp, "label": SUPERPOP.get(sp, sp),
                "proportion": round(n / len(neigh), 4)} for sp, n in tally.items()),
              key=lambda x: -x["proportion"])
top = comp[0]

out = {
    "kind": "ancestry", "status": "ok", "generated_utc": os.environ["GEN"],
    "sample": os.environ.get("SAMPLE", "unknown"),
    "method": f"PCA projection onto HGDP+1kGP + k-NN (k={knn}) over {npc} PCs (plink2)",
    "reference_panel": "HGDP+1kGP (PGS Catalog resources), autosomal LD-pruned common SNPs",
    "n_pcs": npc,
    "assigned_superpopulation": top["superpopulation"],
    "assigned_label": top["label"],
    "confidence": top["proportion"],
    "composition": comp,
    "pcs": {f"PC{i+1}": round(v, 4) for i, v in enumerate(svec)},
    "caveats": ("Ancestry is a continuum, not discrete boxes. These are k-NN vote "
                "shares against one reference panel — reference-relative estimates, "
                "NOT admixture fractions or identity. For formal local-ancestry / "
                "admixture use RFMix or ADMIXTURE. Recent-admixed genomes may split "
                "across several groups."),
}
json.dump(out, open(out_f, "w"), indent=2)
print(f"  assigned={top['superpopulation']} ({top['label']}) conf={top['proportion']:.2f}  "
      f"composition={[c['superpopulation']+':'+str(c['proportion']) for c in comp]}")
PY

ok "Ancestry written -> $OUTJSON"
echo "  Rebuild + serve:  scripts/50_serve_mcp.sh --rebuild"
