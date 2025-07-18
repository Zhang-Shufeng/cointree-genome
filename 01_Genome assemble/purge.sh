#!/bin/sh

for fas in *_finally_rename_all.fas; do
    col1=${fas%%_finally_rename_all.fas}  
    echo "${col1}"
    mkdir "${col1}"
    cd "${col1}"

    ln -s /work/00_Hifi_reassembly/01_HiFi_seq/02-Hifi_Revio_fastq/26-RG.fastq.gz
    ln -s ../${fas}

    minimap2 -x map-hifi -c -o ${col1}.paf --MD -t 80 --secondary=no  ${fas} 26-RG.fastq.gz 2> pbcns.log
    gzip ${col1}.paf

    pbcstat *.paf.gz
    calcuts PB.stat > cutoffs 2> calcults.log

    split_fa ${fas} > asm.split
    minimap2 -t 80 -xasm5 -DP asm.split asm.split | pigz -c > asm.split.self.paf.gz

    purge_dups -2 -T cutoffs -c PB.base.cov asm.split.self.paf.gz > ${col1}_dups.bed 2> purge_dups.log

    cd ..
done
