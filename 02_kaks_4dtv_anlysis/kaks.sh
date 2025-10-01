#!/bin/bash

# ======================= CONFIG ===========================
INPUT_LIST="sample_list.txt"  
DATA_DIR="../01_data"
BLOCK_DIR="../03_blocks"
OUTPUT_DIR="../05_kaks_resoult"

THREADS=10
ALIGN_METHOD="mafft"
KAKS_METHOD="YN"

# ======================= START PIPELINE ====================
set -e 

command -v ParaAT.pl >/dev/null 2>&1 || { echo "ParaAT.pl not found. Please check your environment."; exit 1; }
command -v KaKs_Calculator >/dev/null 2>&1 || { echo "KaKs_Calculator not found."; exit 1; }
command -v seqkit >/dev/null 2>&1 || { echo "seqkit not found."; exit 1; }

mkdir -p "${OUTPUT_DIR}"

while read -r sample; do
    echo "==> Processing: ${sample}"

    cp "${DATA_DIR}/${sample}_cds.fa" ./

    sed '/^>/s/\(>[^locus]*\)locus=.*\s*/\1/' "${sample}_cds.fa" > "${sample}_cds.fas"
    mv "${sample}_cds.fas" "${sample}_cds.fa"

    seqkit translate --trim "${sample}_cds.fa" > "${sample}_pro.fa"

    awk '{print $1}' "${sample}_cds.fa" > "${sample}_new.cds.fa"
    awk '{print $1}' "${sample}_pro.fa" > "${sample}_new.pro.fa"

    sed -i "s/\./_/g" "${sample}_new.cds.fa"
    sed -i "s/\./_/g" "${sample}_new.pro.fa"

    cp "${BLOCK_DIR}/${sample}.i1_resoult.blocks" ./
    sed -i "s/\./_/g" "${sample}.i1_resoult.blocks"

    echo "${THREADS}" > proc.txt

    # Run ParaAT
    ParaAT.pl -h "${sample}.i1_resoult.blocks" \
              -a "${sample}_new.pro.fa" \
              -n "${sample}_new.cds.fa" \
              -p proc.txt -o "${sample}_align_out" \
              -m "${ALIGN_METHOD}" -f axt

    cat "${sample}_align_out"/*.axt > "${sample}_merge_align.axt"

    KaKs_Calculator -m "${KAKS_METHOD}" \
                    -i "${sample}_merge_align.axt" \
                    -o "${sample}_resoult.txt"

    mv "${sample}_resoult.txt" "${OUTPUT_DIR}/"

    rm -rf "${sample}"*

    echo "==> ${sample} done."

done < "${INPUT_LIST}"

echo "==> All samples processed!"
