#!/bin/sh
do
	cp ../01_data/${col1}_cds.fa ./
	sed '/^>/s/\(>[^locus]*\)locus=.*\s*/\1/' ${col1}_cds.fa > ${col1}_cds.fas
	mv ${col1}_cds.fas ${col1}_cds.fa
	seqkit translate --trim ${col1}_cds.fa > ${col1}_pro.fa
	awk '{print $1}' ${col1}_cds.fa > ${col1}_new.cds.fa
	awk '{print $1}' ${col1}_pro.fa > ${col1}_new.pro.fa
	sed -i "s/\./_/g" ${col1}_new.cds.fa
	sed -i "s/\./_/g" ${col1}_new.pro.fa
	cp ../03_blocks/${col1}.i1_resoult.blocks ./
	sed -i "s/\./_/g" ${col1}.i1_resoult.blocks
	echo "10" > proc.txt
	ParaAT.pl -h ${col1}.i1_resoult.blocks -a ${col1}_new.pro.fa -n ${col1}_new.cds.fa -p proc.txt -o ${col1}_align_out -m mafft -f axt
	cat ${col1}_align_out/*.axt > ${col1}_merge_align.axt
	KaKs_Calculator -m YN -i ${col1}_merge_align.axt -o  ${col1}_resoult.txt
	mv ${col1}_resoult.txt ../05_kaks_resoult/
	rm -rf ${col1}*
done
