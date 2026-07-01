# ============================================================================
#  genome-toolkit configuration  —  THE ONE FILE MOST USERS NEED TO EDIT.
#
#  1.  cp config.example.sh config.sh
#  2.  edit the paths below to point at your data
#  3.  run any script — they all read this automatically (via scripts/lib/common.sh)
#
#  `config.sh` is gitignored (it may hold personal paths), so your settings stay
#  local and never get committed. This template is safe to publish.
#
#  Each setting uses  : "${VAR:=value}"  which means: "use VAR from the
#  environment if it's already set, otherwise use this value." So precedence is
#     environment variable  >  this file  >  built-in default
#  — handy for one-off runs (`DATA_DIR=/other ./scripts/run_all.sh`) without
#  editing the file.
# ============================================================================

# --- WHERE YOUR DATA LIVES --------------------------------------------------
# The folder that holds your whole-genome download (the *.vcf.gz / *.fq.gz
# files). This is the main — often the only — thing you need to set.
# Keep it OUTSIDE the repo so your genome never sits next to the code.
: "${DATA_DIR:=$HOME/genome-data}"

# --- OPTIONAL: where results and big downloads go ---------------------------
# Generated reports/indexes (gitignored). Default: <repo>/analysis
# : "${OUT_DIR:=$DATA_DIR/analysis}"
#
# Reference-genome cache (large; downloaded on demand, gitignored).
# Default: <repo>/reference
# : "${REF_DIR:=$DATA_DIR/reference}"

# --- OPTIONAL: point directly at specific files -----------------------------
# Auto-detection handles sequencing.com / DRAGEN / GATK / DeepVariant layouts.
# Only uncomment these if your filenames are unusual and a script can't find them.
# : "${VCF_SNP:=$DATA_DIR/your.snp-indel.genome.vcf.gz}"   # small-variant VCF/gVCF
# : "${VCF_CNV:=$DATA_DIR/your.cnv.vcf.gz}"                 # copy-number calls
# : "${VCF_SV:=$DATA_DIR/your.sv.vcf.gz}"                   # structural variants
# : "${FASTQ_R1:=$DATA_DIR/your_R1.fastq.gz}"              # raw reads, mate 1
# : "${FASTQ_R2:=$DATA_DIR/your_R2.fastq.gz}"              # raw reads, mate 2
