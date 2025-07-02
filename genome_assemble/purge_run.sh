more run/run2.sh 
#!/bin/sh
input_list=$1
cat $input_list |while read -r col1;
do
        echo ${col1}
	mkdir ${col1}
	cd ${col1}
	ln -s /data/02_work/yangm/03_project/01_pangenome26/00_Hifi_reassembly/01_HiFi_seq/02-Hifi_Revio_fastq/26-RG.fastq.gz
	ln -s ../../${col1}_finally_rename_all.fas 
	minimap2 -x map-hifi -c -o ${col1}.paf --MD -t 80 --secondary=no  ${col1}_finally_rename_all.fas 26-RG.fastq.gz 2> pbcns.log
	gzip ${col1}.paf
	pbcstat *.paf.gz
	calcuts PB.stat > cutoffs 2>calcults.log
	split_fa ${col1}_finally_rename_all.fas > asm.split
	minimap2 -t 80 -xasm5 -DP asm.split asm.split | pigz -c > asm.split.self.paf.gz
	purge_dups -2 -T cutoffs -c PB.base.cov asm.split.self.paf.gz > ${col1}_dups.bed 2> purge_dups.log
	cd ..
done

