#!/bin/bash

# Reference genome
FA_REF="./00_DATA/hap.fa"

# Hi-C reads (paired-end)
FQ1="./00_DATA/HiC.1.clean.fq.gz"
FQ2="./00_DATA/HiC.2.clean.fq.gz"

# Thread settings
THREADS=50
SAMTOOLS_THREADS=40

# Juicer tools path
JUICER_TOOLS="/work/user/juicertools/juicer_tools_1.22.01.jar"

# Output prefix
OUT_PREFIX="hap"

# ========== BEGIN PIPELINE ==========

set -e  # exit on error

ln -sf "${FA_REF}" hap.fa
samtools faidx hap.fa
chromap -i -r hap.fa -o contigs.index

chromap --preset hic -r hap.fa -x contigs.index --remove-pcr-duplicate \
    -1 "${FQ1}" -2 "${FQ2}" --SAM -o aligned.sam -t ${THREADS}

samtools view -@ ${SAMTOOLS_THREADS} -bh aligned.sam > aligned.bam 2> samtools.log

yahs --no-contig-ec --no-scaffold-ec hap.fa aligned.bam

juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final.agp hap.fa.fai

java -Xmx100G -jar "${JUICER_TOOLS}" pre out_JBAT.txt out_JBAT.hic <(echo "assembly") &

echo "==> All steps completed successfully!"
