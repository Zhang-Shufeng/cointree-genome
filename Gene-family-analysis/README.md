# Gene family analysis for cointree
## Step1. Identification of homologous genes 
```

orthofinder -f 01_data -t 30 -a 30 -M msa -S blast -T iqtree
```

## Step2. Divergence time tree construction using homologous gene sequences
```

mafft --auto --thread 12 "$file" > "$output_file"
trimal -in "$file" -out "$output_file" -automated1 -keepheader
iqtree2 -s supergene.fasta -T 14 --alrt 1000 -B 5000
mcmctree mcmctree.ctl
```

## Step3. Detection of gene family expansions and contractions
```
awk -v OFS="\t" '{$NF=null;print $1,$0}' Orthogroups.GeneCount.tsv |sed -E -e 's/Orthogroup/desc/' -e 's/_[^\t]+//g' >gene_families.txt
awk 'NR==1 || $3<100 && $4<100 && $5<100 && $6<100 && $7<100 && $8<100 && $9<100 && $10<100 && $11<100  {print $0}' gene_families.txt >gene_families_filter.txt
cafe5 -i gene_families_filter.txt -t modified_tree.txt -p -k 2 -o k2p
```
