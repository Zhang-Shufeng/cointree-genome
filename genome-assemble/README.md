# T2T genome assemble and annotation 
## Step1. genome assemble 
```
hifiasm -o TQS_hap -t 40 --ul-cut 15000 -D10 --hom-cov 71 --h1 Hic_clean.R1.fastq.gz --h2 Hic_clean.R1.fastq.gz --ul ONT_pass_reads.fasta.gz Hifi.ccs.fastq.gz >hap.log
# the above command generated two preliminary haplotypes: hap.hic.hap1.p_ctg.gfa and hap.hic.hap2.p_ctg.gfa
```

## Step2. The purge for contig genome
```
bash purge.sh input.list
```

## Step3. Anchoring to chromosomes
```
ln -s ../00_DATA/hap.fa
samtools faidx hap.fa
chromap -i -r hap.fa -o contigs.index
chromap --preset hic -r hap.fa -x contigs.index --remove-pcr-duplicate -1 HiC.1.clean.fq.gz -2 HiC.2.clean.fq.gz --SAM -o aligned.sam -t 50
samtools view -@ 40 -bh aligned.sam > aligned.bam 2>2.log 
yahs --no-contig-ec --no-scaffold-ec hap.fa aligned.bam
juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final.agp hap1_400k.fa.fai
nohup java -Xmx100G -jar /data/02_work/zhsf/chengxu/juicertools/juicer_tools_1.22.01.jar pre out_JBAT.txt out_JBAT.hic <(echo "assembly 383783947") &

#Manual correction in Juicerbox

```

## Step4. Fill the gap
```
python quartet.py  GapFiller  -d  hap.fa  -g HiFi.ccs.fasta -t 30 --minimapoption \'-x map-hifi\'

```

## Step5. Telomere and centromere identification
```
perl search_motif.pl -i hap.fa -m motif.txt -o duanli.bed
python3 quartet.py TeloExplorer -i hap.fa -c plant -o TQS_tm

```

## Step6. Genome annotation
```
bash genome_annotation.sh

```