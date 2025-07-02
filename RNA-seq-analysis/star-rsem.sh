mkdir index
python gff3togtf.py hap.final.gff3 > hap.final.gtf
STAR --runThreadN 40 --runMode genomeGenerate --genomeDir ./index  --genomeFastaFiles hap
ifi.final.fas --sjdbGTFfile hap.final.gtf --sjdbOverhang 149
cat hap.final.gff3 |grep -v "exon" |grep -v "mRNA" > 1.gff3
cp /data/02_work/zhsf/wenjian/12_D_R_hap_new/04_RNA_seq/03_STAR/01_DZ/tq_intron_length.py ./
python tq_intron_length.py 1.gff3 1.txt
more 1.txt |grep -v "0" |sort -k1,1nr |tail -1
more 1.txt |grep -v "0" |sort -k1,1nr |head -1
STAR --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --runThreadN 20 --genomeDir index --alignIntronMin 43 --alignIntronMax 89634 --outSAMtype BAM SortedByCoordinate --sjdbOverhang 149 --outSAMattrRGline ID:sample SM:sample PL:ILL
UMINA --outFilterMismatchNmax 1 --outSAMmultNmax -1 --outFileNamePrefix hap --readFilesCommand gunzip -c --readFilesIn RNA.clean.R1.fq.gz RNA.clean.R2.fq.gz
rsem-prepare-reference --gtf hap.final.gtf hap.Hifi.final.fas hap.Hifi.final -p 8
