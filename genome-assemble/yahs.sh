ln -s ../00_DATA/hap.fa
samtools faidx hap.fa
chromap -i -r hap.fa -o contigs.index
chromap --preset hic -r hap.fa -x contigs.index --remove-pcr-duplicate -1 HiC.1.clean.fq.gz -2 HiC.2.clean.fq.gz --SAM -o aligned.sam -t 50
samtools view -@ 40 -bh aligned.sam > aligned.bam 2>2.log 
yahs --no-contig-ec --no-scaffold-ec hap.fa aligned.bam
juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final.agp hap1_400k.fa.fai
nohup java -Xmx100G -jar /data/02_work/zhsf/chengxu/juicertools/juicer_tools_1.22.01.jar pre out_JBAT.txt out_JBAT.hic <(echo "assembly 383783947") &

