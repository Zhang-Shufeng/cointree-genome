#!/bin/bash
set -e  

THREADS=40
GENOME_FASTA="hapHifi.final.fas"
GFF3_FILE="hap.final.gff3"
GTF_FILE="hap.final.gtf"
INDEX_DIR="./index"
RNA_READ1="RNA.clean.R1.fq.gz"
RNA_READ2="RNA.clean.R2.fq.gz"
INTRLENGTH_OUTPUT="1.txt"
STAR_PREFIX="hap"
RSEM_PREFIX="hap.Hifi.final"

mkdir -p ${INDEX_DIR}

echo "Converting GFF3 to GTF..."
python gff3togtf.py ${GFF3_FILE} > ${GTF_FILE}

echo "Generating STAR genome index..."
STAR --runThreadN ${THREADS} \
     --runMode genomeGenerate \
     --genomeDir ${INDEX_DIR} \
     --genomeFastaFiles ${GENOME_FASTA} \
     --sjdbGTFfile ${GTF_FILE} \
     --sjdbOverhang 149

echo "Filtering GFF3 annotations..."
grep -v -e "exon" -e "mRNA" ${GFF3_FILE} > 1.gff3

cp /work/01_RNA_seq/03_STAR/${PYTHON_SCRIPT_INTRLENGTH} ./

echo "Calculating intron lengths..."
python ${PYTHON_SCRIPT_INTRLENGTH} 1.gff3 ${INTRLENGTH_OUTPUT}

echo "Intron length max:"
grep -v "0" ${INTRLENGTH_OUTPUT} | sort -k1,1nr | tail -1

echo "Intron length min:"
grep -v "0" ${INTRLENGTH_OUTPUT} | sort -k1,1nr | head -1

echo "Running STAR two-pass alignment..."
STAR --twopassMode Basic \
     --quantMode TranscriptomeSAM GeneCounts \
     --runThreadN 20 \
     --genomeDir ${INDEX_DIR} \
     --alignIntronMin 43 \
     --alignIntronMax 89634 \
     --outSAMtype BAM SortedByCoordinate \
     --sjdbOverhang 149 \
     --outSAMattrRGline ID:sample SM:sample PL:ILLUMINA \
     --outFilterMismatchNmax 1 \
     --outSAMmultNmax -1 \
     --outFileNamePrefix ${STAR_PREFIX} \
     --readFilesCommand gunzip -c \
     --readFilesIn ${RNA_READ1} ${RNA_READ2}

echo "Preparing RSEM reference..."
rsem-prepare-reference --gtf ${GTF_FILE} ${GENOME_FASTA} ${RSEM_PREFIX} -p 8

echo "All steps finished."