#!/usr/bin/env bash

BASE_VCF="/home/maycon/projects/turtles/work/bcftools/All.filtered.sorted.renamed.normalised.vcf.gz"
POP_FILE="/home/maycon/projects/turtles/work/pi/pures_popmap.tsv"

echo "pop,avg_pi" > summary.csv

for spp in Cc_BR Cm Ei Lo Cc_ERG; do
	grep "$spp" "$POP_FILE" | cut -f 1 > "$spp"_samples.txt
	bcftools view -S "$spp"_samples.txt -Oz -o "$spp".All.vcf.gz "$BASE_VCF"
	tabix "$spp".All.vcf.gz
	vcftools --gzvcf "$spp".All.vcf.gz --site-pi --out "$spp".All
	avg_pi=$(awk 'NR>1 { if ($3 != "-nan") { sum+=$3; sites+=1 } } END { print sum/sites }' "$spp".All.sites.pi)
	echo "$spp,$avg_pi" >> summary.csv
done
