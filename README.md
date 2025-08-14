# GWAS → GO Enrichment Pipeline

This repository provides a **one-command, end-to-end pipeline** to perform GO enrichment analysis given your SNP(s) and an annotated genome file. The script produces ready-to-use enrichment tables and plots.

## Prerequisites

- **Input Genome Annotation**:  
  You must have an annotated GFF3 file (e.g., `bcar_annotated.gff3`) that includes `GO=` attributes for genes.  
  *If your file is differently named, adjust the `-g` argument in the command.*

- **R environment**:  
  The script automatically installs required R packages as needed.

## Quick Start

### 1. Save the Pipeline Script

Copy the contents below into a file named `run_gwas_go.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

# ================================
# GWAS → GO Enrichment (One-Run)
# Inputs:
#   - bcar_annotated.gff3  (with GO terms)
#   - Either a single SNP (CHR + POS) or a TSV of SNPs (chr\tpos)
# Outputs:
#   - gwas_genes_unique.txt, all_genes_around_snps.tsv
#   - background_genes.txt, gene2go.map
#   - go_enrichment_* TSVs and PNG plots
# ================================

usage() {
  cat <<EOF
Usage:
  $(basename "$0") -g bcar_annotated.gff3 -w WINDOW_BP [-c CHR -p POS] [-s snps.tsv] [-o OUTDIR]

Required:
  -g    Annotated GFF3 file (e.g., bcar_annotated.gff3 with GO= attributes)
  -w    Window size in bp on each side (e.g., 100000 for ±100kb)

Choose ONE of:
  -c    Chromosome name (e.g., CM081020.1) AND -p position (integer)
  -s    SNP TSV with 2 columns: chr<TAB>pos (header allowed)

Optional:
  -o    Output directory (default: ./go_enrichment_run)

Examples:
  # Single SNP ±100kb:
  ./run_gwas_go.sh -g bcar_annotated.gff3 -w 100000 -c CM081020.1 -p 31630794

  # Many SNPs from a file ±50kb:
  ./run_gwas_go.sh -g bcar_annotated.gff3 -w 50000 -s snps.tsv
EOF
}

# ---------- Parse args ----------
GFF=""; WIN=""; CHR=""; POS=""; SNPFILE=""; OUTDIR="./go_enrichment_run"
while getopts ":g:w:c:p:s:o:h" opt; do
  case $opt in
    g) GFF="$OPTARG" ;;
    w) WIN="$OPTARG" ;;
    c) CHR="$OPTARG" ;;
    p) POS="$OPTARG" ;;
    s) SNPFILE="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    h) usage; exit 0 ;;
    \?) echo "Unknown option: -$OPTARG" >&2; usage; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; usage; exit 1 ;;
  esac
done

[[ -z "$GFF" || -z "$WIN" ]] && { echo "ERROR: -g and -w are required."; usage; exit 1; }
[[ -s "$GFF" ]] || { echo "ERROR: GFF not found: $GFF" >&2; exit 1; }

if [[ -n "${CHR}" && -n "${POS}" && -z "${SNPFILE}" ]]; then
  MODE="single"
elif [[ -z "${CHR}" && -z "${POS}" && -n "${SNPFILE}" ]]; then
  MODE="file"
  [[ -s "$SNPFILE" ]] || { echo "ERROR: SNP file not found: $SNPFILE" >&2; exit 1; }
else
  echo "ERROR: Provide either (-c & -p) OR (-s), not both." >&2
  usage; exit 1
fi

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Timing helpers
ts() { date +"%Y-%m-%d %H:%M:%S"; }
STEP_START=$SECONDS
echo "[`ts`] Starting run. Output dir: $PWD"
echo "[`ts`] Inputs: GFF=$GFF | WIN=$WIN | MODE=$MODE"

# ---------- Step 1: Build SNP list ----------
if [[ "$MODE" == "single" ]]; then
  echo -e "${CHR}\t${POS}" > snps.tsv
  echo "[`ts`] Wrote snps.tsv with 1 SNP."
else
  # normalize to two-column CHR TAB POS without header
  awk 'BEGIN{OFS="\t"} NR==1{if($2 ~ /^[0-9]+$/) print $1,$2; next} {if($2 ~ /^[0-9]+$/) print $1,$2}' "$SNPFILE" > snps.tsv
  echo "[`ts`] Wrote snps.tsv with " $(wc -l < snps.tsv) " rows."
fi

# ---------- Step 2: Extract genes around SNPs ----------
echo "[`ts`] Extracting genes within ±$WIN bp of SNPs..."
GENE_LIST="gwas_genes_unique.txt"
FULL_TSV="all_genes_around_snps.tsv"
: > "$GENE_LIST"
echo -e "snp_chr\tsnp_pos\tchr\tstart\tend\tgene_id\tlocus_tag\tortholog\tdescription\tGO\tGO_names" > "$FULL_TSV"

while IFS=$'\t' read -r SNP_CHR SNP_POS; do
  [[ -z "$SNP_CHR" || -z "$SNP_POS" ]] && continue
  [[ "$SNP_POS" =~ ^[0-9]+$ ]] || continue
  START=$(( SNP_POS - WIN )); (( START < 0 )) && START=0
  END=$(( SNP_POS + WIN ))

  awk -v c="$SNP_CHR" -v s="$START" -v e="$END" -v sc="$SNP_CHR" -v sp="$SNP_POS" -F'\t' '
    $1==c && $3=="gene" && $4<=e && $5>=s {
      split($9,a,";");
      id=""; lt=""; ort=""; desc=""; go=""; gon="";
      for(i in a){
        split(a[i],kv,"=");
        k=kv[1]; v=kv[2];
        gsub(/%([0-9A-Fa-f]{2})/,"\\x\\1",v);  # URL-decode
        if(k=="ID") id=v;
        else if(k=="locus_tag") lt=v;
        else if(k=="Ortholog")  ort=v;
        else if(k=="Description") desc=v;
        else if(k=="GO")        go=v;
        else if(k=="GO_names")  gon=v;
      }
      if(id!=""){
        print id >> "__tmp_genes.txt";
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
               sc, sp, $1, $4, $5, id, lt, ort, desc, go, gon >> "__tmp_full.tsv";
      }
    }
  ' "$GFF"
done < snps.tsv

if [[ -f "__tmp_genes.txt" ]]; then
  sort -u "__tmp_genes.txt" > "$GENE_LIST"
  if [[ -f "__tmp_full.tsv" ]]; then cat "__tmp_full.tsv" >> "$FULL_TSV"; fi
  rm -f __tmp_genes.txt __tmp_full.tsv
else
  echo "[`ts`] WARNING: No genes found within windows."
fi
echo "[`ts`] Genes for enrichment: $(wc -l < "$GENE_LIST")"
echo "[`ts`] Full region table: $FULL_TSV"

# ---------- Step 3: Build background and gene→GO ----------
echo "[`ts`] Building background and gene→GO mapping..."
grep -P "\tgene\t" "$GFF" \
| awk -F'\t' '
  {
    split($9,a,";");
    id=""; go="";
    for(i in a){
      split(a[i],kv,"=");
      if(kv[1]=="ID") id=kv[2];
      else if(kv[1]=="GO") go=kv[2];
    }
    if(id!="" && go!="") print id;
  }
' | sort -u > background_genes.txt

grep -P "\tgene\t" "$GFF" \
| awk -F'\t' '
  {
    split($9,a,";");
    id=""; go="";
    for(i in a){
      split(a[i],kv,"=");
      if(kv[1]=="ID") id=kv[2];
      else if(kv[1]=="GO") go=kv[2];
    }
    if(id!="" && go!="") print id "\t" go;
  }
' > gene2go.map

echo "[`ts`] Universe size (with GO): $(wc -l < background_genes.txt)"

# ---------- Step 4: Enrichment (R script is generated on the fly) ----------
cat > 3_go_enrichment.R <<'RSCRIPT'
... # [Omitted for brevity: see script in this repository]
RSCRIPT

echo "[`ts`] Running enrichment (R)..."
Rscript 3_go_enrichment.R

# ---------- Done & timings ----------
TOTAL_SEC=$(( SECONDS - STEP_START ))
echo "[`ts`] All done. Total runtime: ${TOTAL_SEC}s"
echo "Outputs in: $PWD"
ls -1 go_enrichment_*.* || true
```

Make it executable:
```bash
chmod +x run_gwas_go.sh
```

---

## 2. Example Runs

**Single SNP, ±100 kb:**
```bash
./run_gwas_go.sh -g /path/to/bcar_annotated.gff3 -w 100000 -c CM081020.1 -p 31630794
```

**Multiple SNPs file, ±50 kb:**
```bash
# snps.tsv should have two columns: chr <TAB> pos
./run_gwas_go.sh -g /path/to/bcar_annotated.gff3 -w 50000 -s snps.tsv
```

---

## 3. Output Files

Inside the output folder (default: `./go_enrichment_run`):

- `gwas_genes_unique.txt` — unique genes near your SNP(s)
- `all_genes_around_snps.tsv` — full decoded table of all genes in windows
- `background_genes.txt` — universe of genes with GO
- `gene2go.map` — gene → GO mapping
- `go_enrichment_topGO_BP/MF/CC.tsv` + `_barplot.png`
- `go_enrichment_clusterProfiler_ALL/BP/MF/CC.tsv` + `_dotplot.png`

All plots and tables are robust, readable, and ready for downstream analysis.

---

## Notes

- **Runtime:**  
  Most time will be spent on first-time R package installs. Subsequent runs are fast (seconds to a couple minutes).

- **Warnings:**  
  If no genes are found in your window, the script will warn but still produce available outputs.  
  *Increase `-w` window size if you expect more overlap.*

- **Extending:**  
  If you want to automate the annotation-building step (e.g., Arabidopsis RBH → `bcar_annotated.gff3`) as part of this pipeline, open an Issue or PR!

---

## License

MIT

---

## Citation

If you use this pipeline, please cite this repository or its author.

```
Absolutely—here’s a one-command, end-to-end pipeline you can run that takes your SNP(s) + window and produces the GO enrichment tables and plots.
```
