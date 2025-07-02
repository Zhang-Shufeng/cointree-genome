# NLR genes annotation 
## Step1. Searching for the NLR genes
```
java -Xmx8000M -jar /data/02_work/zhsf/chengxu/R-gene/NLR-Annotator/NLR-Annotator-v2.1b.jar -i TQS_hap1.cds.fas  -x /data/02_work/zhsf/chengx
u/R-gene/NLR-Annotator/src/mot.txt -y /data/02_work/zhsf/chengxu/R-gene/NLR-Annotator/src/store.txt -o output.txt -g output.gff -b output.bed -m ou
tput.motifs.bed -a genome_NLR.nbarcMotifAlignment.fasta -t 40 
```

## Step2. Classification of NLR genes
```

grep -f R-gene_pre.list ../hap.ipr |grep "NB-ARC"|cut -f 1|sort|uniq > NB-ARC.gene.list
grep -f R-gene_pre.list ../hap.ipr |grep "TIR"|cut -f 1|sort|uniq > TIR.gene.list
grep -f R-gene_pre.list ../hap.ipr |grep "RPW8"|cut -f 1|sort|uniq > RPW8.gene.list
cat NB-ARC.gene.list TIR.gene.list RPW8.gene.list |sort|uniq > R-gene_resoult.list
grep -f R-gene_resoult.list output.gff |grep -v "#"|grep "nlr1" |cut -f 9|sed "s/name=//g"|awk -F'[;=]' '{print $1 "\t" $3}'|sed "s/_nlr1//g" > R-g
ene_resoult.bed

grep $'\tTIR-NBARC$' R-gene_resoult.bed > TN.list
grep $'\tCC-NBARC$' R-gene_resoult.bed > CN.list
grep $'\tNBARC-LRR$' R-gene_resoult.bed > NL.list
grep $'\tCC-NBARC-LRR$' R-gene_resoult.bed > CNL.list
grep $'\tNBARC$' R-gene_resoult.bed > NB.list
grep $'\tTIR-NBARC-LRR$' R-gene_resoult.bed > TNL.list

```

