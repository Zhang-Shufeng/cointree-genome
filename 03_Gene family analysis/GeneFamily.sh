#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status

# ==== Configuration ====
DATA_DIR="01_data"
THREADS=30
MAFFT_THREADS=12
IQTREE_THREADS=14
CAFE_INPUT="gene_families_filter.txt"
CAFE_TREE="modified_tree.txt"
CAFE_OUTPUT="k2p"

echo "==> Running OrthoFinder..."
orthofinder -f ${DATA_DIR} -t ${THREADS} -a ${THREADS} -M msa -S blast -T iqtree

echo "==> Running MAFFT alignment (example with single file)..."
# If you have multiple files, use a loop like this:
# for file in input_folder/*.fasta; do
#     output_file="aligned/$(basename "$file")"
#     mafft --auto --thread ${MAFFT_THREADS} "$file" > "$output_file"
# done

# Example for a single file alignment
file="example_input.fasta"
output_file="example_aligned.fasta"
mafft --auto --thread ${MAFFT_THREADS} "$file" > "$output_file"

echo "==> Running trimAl to trim alignments..."
trimal -in "$file" -out "$output_file" -automated1 -keepheader

echo "==> Running IQ-TREE phylogenetic analysis..."
iqtree2 -s supergene.fasta -T ${IQTREE_THREADS} --alrt 1000 -B 5000

echo "==> Running MCMCTree..."
mcmctree mcmctree.ctl

echo "==> Processing Orthogroups gene count file..."
awk -v OFS='\t' '{$NF=null; print $1,$0}' Orthogroups.GeneCount.tsv | \
    sed -E -e 's/Orthogroup/desc/' -e 's/_[^\t]+//g' > gene_families.txt

echo "==> Filtering gene families..."
awk 'NR==1 || ($3<100 && $4<100 && $5<100 && $6<100 && $7<100 && $8<100 && $9<100 && $10<100 && $11<100) {print $0}' gene_families.txt > ${CAFE_INPUT}

echo "==> Running CAFE5 for gene family expansion/contraction analysis..."
./cafe5 -i ${CAFE_INPUT} -t ${CAFE_TREE} -p -k 2 -o ${CAFE_OUTPUT}