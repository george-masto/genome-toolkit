#!/usr/bin/env python3
"""
build_embeddings.py — a LOCAL semantic-search index over your genome findings.

Embeddings are computed ON-DEVICE (fastembed / ONNX, CPU). Only the small
PUBLIC embedding model is downloaded once; your genome text is never sent
anywhere. With ~1-2k findings we skip a vector DB entirely — a normalized
NumPy matrix + a dot product is millisecond-fast.

Outputs (gitignored, in analysis/):
  genome_embeddings.npz     float32 vectors (one row per finding)
  genome_embeddings.jsonl   parallel metadata (what each vector points to)

Usage:  python3 mcp/build_embeddings.py
Needs:  pip install -r mcp/requirements-embeddings.txt ; build_insights.py first.
"""
import json, os, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("OUT_DIR", REPO / "analysis"))
INSIGHTS = OUT / "genome_insights.json"
CLINVAR = OUT / "20_clinvar_all_categorized.tsv"
MODEL = os.environ.get("EMBED_MODEL", "BAAI/bge-small-en-v1.5")
INCLUDE_BENIGN = os.environ.get("EMBED_INCLUDE_BENIGN", "0") == "1"


def build_corpus():
    """One 'document' per finding: a short text to embed + metadata to return."""
    docs = []
    data = json.loads(INSIGHTS.read_text()) if INSIGHTS.exists() else {}

    for t in data.get("traits", []):
        txt = f"{t['gene']} — {t['trait']}. {t['effect']}"
        docs.append((txt, {"source": "trait", "gene": t["gene"], "rsid": t["rsid"],
                            "genotype": t["genotype"], "trait": t["trait"],
                            "detail": t["effect"]}))

    for g in data.get("pharmacogenomics", []):
        txt = (f"{g['gene']} {g['phenotype']}. Affects drugs: "
               f"{', '.join(g['drugs'])}. {g['implication']}")
        docs.append((txt, {"source": "pharmacogenomics", "gene": g["gene"],
                            "genotype": g["genotype"], "phenotype": g["phenotype"],
                            "drugs": g["drugs"], "detail": g["implication"]}))

    if CLINVAR.exists():
        for line in CLINVAR.read_text().splitlines():
            p = line.split("\t")
            if len(p) < 10:
                continue
            bucket, stars, chrom, pos, gene, rsid, geno, zyg, sig, cond = p[:10]
            if bucket == "9_Benign" and not INCLUDE_BENIGN:
                continue
            condition = cond.replace("_", " ")
            txt = f"{gene} — {condition} ({sig.replace('_', ' ')})"
            docs.append((txt, {"source": "clinvar", "tier": bucket.split("_", 1)[-1],
                               "stars": int(stars), "gene": gene, "rsid": rsid,
                               "genotype": geno, "zygosity": zyg,
                               "significance": sig, "condition": condition}))
    return docs


def main():
    if not INSIGHTS.exists():
        sys.exit("Missing analysis/genome_insights.json — run mcp/build_insights.py first.")
    try:
        import numpy as np
        from fastembed import TextEmbedding
    except ImportError:
        sys.exit("Install embedding deps:  pip install -r mcp/requirements-embeddings.txt")

    docs = build_corpus()
    if not docs:
        sys.exit("No findings to embed. Run scripts/run_all.sh (+ 20_clinvar_annotate.sh) first.")
    texts = [d[0] for d in docs]
    print(f"Embedding {len(texts)} findings locally with {MODEL} "
          f"(model downloads once, then offline)…", file=sys.stderr)

    model = TextEmbedding(model_name=MODEL)
    vecs = np.asarray(list(model.embed(texts)), dtype=np.float32)
    vecs /= (np.linalg.norm(vecs, axis=1, keepdims=True) + 1e-9)  # normalize → cosine via dot

    np.savez_compressed(OUT / "genome_embeddings.npz", vectors=vecs, model=MODEL)
    with open(OUT / "genome_embeddings.jsonl", "w") as f:
        for _, meta in docs:
            f.write(json.dumps(meta) + "\n")
    print(f"Wrote {OUT/'genome_embeddings.npz'}  ({vecs.shape[0]} findings × {vecs.shape[1]} dims)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
