# T2T genome assemble and annotation 
## Step1. genome assemble 
```
hifiasm -o TQS_hap -t 40 --ul-cut 15000 -D10 --hom-cov 71 --h1 Hic_clean.R1.fastq.gz --h2 Hic_clean.R1.fastq.gz --ul ONT_pass_reads.fasta.gz Hifi.ccs.fastq.gz >hap.log
# the above command generated two preliminary haplotypes: hap.hic.hap1.p_ctg.gfa and hap.hic.hap2.p_ctg.gfa
```

## Step2. Eliminating redundant or overlapping contigs from the assembled genome
```
bash purge.sh 
```

## Step3. Anchoring to chromosomes
```
bash yahs.sh

#Manual correction in Juicerbox

```

## Step4. Fill the gap
```
python quartet.py  GapFiller  -d  hap.fa  -g HiFi.ccs.fasta -t 30 --minimapoption \'-x map-hifi\'

```

## Step5. Telomere and centromere identification
```
perl search_motif.pl -i input.fa -m motif.txt -o telomere.bed
python quartet.py TeloExplorer -i input.fa -c plant -o centromere

```

## Step6. Genome annotation
```
bash genome_annotation.sh

```