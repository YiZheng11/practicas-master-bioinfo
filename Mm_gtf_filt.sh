#!/bin/bash

gtf=$1
filt_genes_bed=$2
assembly_report=$3

# Create TSV file of chromosome to RefSeq ID
echo -e "Creating file chr_to_refseq.tsv..."
awk '{
	if ($1 !~ /^#/) {
		print $10"\t"$7
	}
}' ./data/$assembly_report> ./data/chr_to_refseq.tsv
sed -i 's/\r//g' ./data/chr_to_refseq.tsv
sed -i 's/v2//g' ./data/chr_to_refseq.tsv

echo -e "Finished!"

# Modifying MitoCarta bed file, chr1 -> RefSeq ID
echo -e "Modifying $filt_genes_bed using chr_to_refseq.tsv..."
awk '
	NR==FNR{
		map[$1]=$2;
		next
	}
	{ if ($1 in map) $1=map[$1]; print}
' OFS='\t' ./data/chr_to_refseq.tsv ./data/$filt_genes_bed > ./data/Mouse.MitoCarta3.0_mod.bed

echo -e "Finished modifying, created Mouse.MitoCarta3.0_mod.bed file"

# Filtering GTF file using modified MitoCarta bed file
echo -e "Filtering $gtf file..."
bedtools intersect -a ./data/$gtf -b ./data/Mouse.MitoCarta3.0_mod.bed -wa > GTF_mod.bed

