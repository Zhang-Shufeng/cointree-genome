#!/bin/bash
set -e  

JAVA_MEM="8000M"
NLR_ANNOTATOR_JAR="/work/software/NLR-Annotator/NLR-Annotator-v2.1b.jar"
MOTIF_TXT="/work/software/NLR-Annotator/src/mot.txt"
STORE_TXT="/work/software/NLR-Annotator/src/store.txt"
INPUT_FASTA="cds.fas"
THREADS=40

OUTPUT_PREFIX="output"

echo "==> Running NLR-Annotator..."
java -Xmx${JAVA_MEM} -jar ${NLR_ANNOTATOR_JAR} \
    -i ${INPUT_FASTA} \
    -x ${MOTIF_TXT} \
    -y ${STORE_TXT} \
    -o ${OUTPUT_PREFIX}.txt \
    -g ${OUTPUT_PREFIX}.gff \
    -b ${OUTPUT_PREFIX}.bed \
    -m ${OUTPUT_PREFIX}.motifs.bed \
    -a genome_NLR.nbarcMotifAlignment.fasta \
    -t ${THREADS}

echo "==> Extracting unique Hap IDs from GFF..."
grep "Hap" ${OUTPUT_PREFIX}.gff | cut -f1 | sort | uniq > R-gene_pre.list

echo "==> Filtering gene lists from annotation file..."
grep -f R-gene_pre.list ../hap.ipr | grep "NB-ARC" | cut -f1 | sort | uniq > NB-ARC.gene.list
grep -f R-gene_pre.list ../hap.ipr | grep "TIR"    | cut -f1 | sort | uniq > TIR.gene.list
grep -f R-gene_pre.list ../hap.ipr | grep "RPW8"  | cut -f1 | sort | uniq > RPW8.gene.list

echo "==> Combining R-gene result lists..."
cat NB-ARC.gene.list TIR.gene.list RPW8.gene.list | sort | uniq > R-gene_resoult.list

echo "==> Generating R-gene result bed file..."
grep -f R-gene_resoult.list ${OUTPUT_PREFIX}.gff | grep -v "#" | grep "nlr1" | \
    cut -f9 | sed "s/name=//g" | awk -F'[;=]' '{print $1 "\t" $3}' | sed "s/_nlr1//g" > R-gene_resoult.bed

echo "==> Splitting R-gene results into categories..."

grep $'\tTIR-NBARC$' R-gene_resoult.bed > TN.list
grep $'\tCC-NBARC$' R-gene_resoult.bed > CN.list
grep $'\tNBARC-LRR$' R-gene_resoult.bed > NL.list
grep $'\tCC-NBARC-LRR$' R-gene_resoult.bed > CNL.list
grep $'\tNBARC$' R-gene_resoult.bed > NB.list
grep $'\tTIR-NBARC-LRR$' R-gene_resoult.bed > TNL.list

echo "==> All done."