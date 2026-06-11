# Second pass: what the reads know that the gVCF doesn't

The quick-start pipeline (`run_all.sh`) reads a **flat, unphased SNP/indel
gVCF**. That's the right tool for ~99% of your genome — but a handful of the
most clinically important loci are *structurally* hard: they sit next to a
near-identical pseudogene, vary in copy number, are an HLA gene with thousands
of alleles, or are a tandem repeat that no SNP caller even looks at. For those,
a generic caller returns a dead-end ("copy-number not callable") or, worse, a
confidently wrong answer.

The fix is the same one clinical labs use: **re-align the raw 30× reads to a
CRAM, then run a purpose-built, paralog-aware caller per hard locus.** The CRAM
(`scripts/10_align_from_fastq.sh`) is the one-time cost; it unlocks the whole
panel below, not just CYP2D6.

```
  FASTQ ──10_align_from_fastq.sh──▶ CRAM ──┬─ 11_cyp2d6.sh   (Cyrius)        ─┐
   (hours, ~80 GB, one time)               ├─ 13_hla.sh      (HLA caller)    │   analysis/calls/*.json
                                           ├─ 14_repeats.sh  (ExpansionHunter)├──▶  build_insights.py
                                           ├─ 15_smn.sh      (SMNCNCaller)   │      (ingests + translates)
                                           ├─ 16_gba.sh      (Gauchian)      │            │
                                           └─ 17_mtdna.sh    (Haplogrep+RNR1)┘            ▼
                                                                                   MCP self-upgrades
```

Each caller writes a **normalized `analysis/calls/<name>.json`**. `build_insights.py`
ingests whatever is present and degrades gracefully when a caller hasn't run —
so the MCP shows real calls the moment you rebuild, and "pending" until then
(see the `get_second_pass_status` tool).

---

## The panel — why each locus needs a specialist

| Script | Tool | Locus | Why the gVCF fails | What you get | Clinically because… |
|---|---|---|---|---|---|
| `11_cyp2d6.sh` | **Cyrius** | CYP2D6 | CYP2D7 paralog + copy-number/hybrid alleles | full diplotype + CN → activity-score phenotype | codeine/tramadol, **tamoxifen**, atomoxetine, many antidepressants/antipsychotics |
| `13_hla.sh` | HLA caller | HLA-A/B/C… | thousands of alleles; SNPs don't resolve 2-field types | `HLA-B*57:01`, `B*15:02`, `A*31:01`, `B*58:01` status | **abacavir** & **carbamazepine** & **allopurinol** hypersensitivity (CPIC contraindications) |
| `14_repeats.sh` | **ExpansionHunter** | STR loci | SNP callers never look at tandem repeats | repeat counts at HTT, C9orf72, FMR1, DM1, SCAs, FXN | Huntington's, ALS/FTD, Fragile-X, myotonic dystrophy, ataxias |
| `15_smn.sh` | **SMNCopyNumberCaller** | SMN1 / SMN2 | SMN1≈SMN2 paralog; it's a copy-number question | SMN1 CN (0=affected, 1=**carrier**, 2=typical) + SMN2 CN | spinal muscular atrophy — one of the most common lethal recessive carrier states |
| `16_gba.sh` | **Gauchian** | GBA | GBAP1 pseudogene + recombinant alleles | GBA genotype incl. RecNciI, p.L444P/p.N370S | Gaucher disease (recessive) **and** the strongest common genetic Parkinson's risk |
| `17_mtdna.sh` | Haplogrep + bcftools | chrM | ploidy/heteroplasmy ≠ diploid SNP model | haplogroup, **MT-RNR1 m.1555A>G/m.1494C>T**, heteroplasmy | aminoglycoside-induced hearing loss (CPIC); maternal lineage; contamination QC |

Optional cross-check layer (not gap-filling, but corroborating):

- **Aldy 4** — re-calls CYP2D6 *and* a dozen other PGx genes (CYP2A6, 2B6, 2C8,
  3A4/5, DPYD, TPMT, SLCO1B1, UGT1A1, NUDT15…) copy-number-aware. Run it to
  *agree-or-flag* against Cyrius and the SNP panel.
- **PharmCAT** — can ingest the Cyrius CYP2D6 diplotype as an "outside call" and
  emit the full official CPIC report. We keep the hand-rolled engine for the
  SNP-callable genes (it's transparent and already yours) and use the specialty
  callers for the genes it *can't* do; PharmCAT is the optional capstone if you
  want the canonical report format.

---

## Why CYP2D6 is the headliner (and the honest caveat)

CYP2D6 is the single highest-value upgrade because it governs the largest drug
list and is the one gene the SNP panel openly can't type. The activity-score
translation (CPIC 2019 consensus) lives in
`assets/pharmacogenes/cyp2d6_allele_functions.tsv`; `build_insights.py` turns a
Cyrius diplotype like `*4/*41` into "Intermediate metabolizer (AS 0.25)".

Caveat worth keeping attached: for **paroxetine specifically**, CYP2D6 genotype
carries less prospective weight than its reputation suggests, because paroxetine
inhibits CYP2D6 and phenoconverts most people toward poor-metabolizer at steady
state regardless of starting genotype. The genotype matters far more for
**tamoxifen, codeine/tramadol, atomoxetine**, and several antidepressants/
antipsychotics. So we do this for the *whole* CYP2D6 drug list, not paroxetine.

---

## Status & honesty

Everything here is **research-grade**, predicted from consumer 30× WGS. Reported
concordance with clinical/long-read assays is high (Cyrius/SMN ~the high-90s%),
but a prescriber writing an actual script still confirms in a CLIA lab. The MCP
labels every call with its `method` (`SNP panel` vs `Cyrius`/etc.) so an agent
never mistakes a predicted call for a clinical one.
