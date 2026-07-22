export REF="/home/maycon/projects/turtles/data/ref/rDerCor1.pri.v4.fna"
export REF_INDEX="/home/maycon/projects/turtles/data/ref/rDerCor1.pri.v4.fna.fai"
export BAM_LIST="/home/maycon/projects/turtles/work/bcftools/samples_depth5_nodups_dc.txt"
export OUTDIR="/home/maycon/projects/turtles/work/bcftools"
export REPEATS_BED="/home/maycon/projects/turtles/data/ref/rDerCor1.pri.v4.repeats.bed"
export WINDOWS_BED="/home/maycon/projects/turtles/data/ref/rDerCor1.pri.v4.20mb_windows.5mb_chroms.bed"

# VCFtools filters
export MIN_DP=10
export MAX_DP=100
export MAX_MISS=0.8
export MIN_Q=20

do_call() {
	jobname="${1}_${2}.dc"

	echo -e "$1\t$2\t$3" > "$OUTDIR/$jobname.bed"

	bcftools mpileup -Ou -b "$BAM_LIST" -f "$REF" -R "$OUTDIR/$jobname.bed" -a DP,AD |
		bcftools call -m -Ou -f GQ,GP |
		bcftools view -T ^"$REPEATS_BED" \
		-Oz -o "$OUTDIR/$jobname.no_repeats.vcf.gz"

	vcftools \
		--gzvcf "$OUTDIR/$jobname.no_repeats.vcf.gz" \
		--min-meanDP "$MIN_DP" \
		--max-meanDP "$MAX_DP" \
		--max-missing "$MAX_MISS" \
		--minQ "$MIN_Q" \
		--remove-indels \
		--max-alleles 2 \
		--recode \
		--stdout | bgzip -c > "$OUTDIR/$jobname.no_repeats.filtered.vcf.gz"
}
export -f do_call

# BEWARE: this is a tab chacter on colsep
parallel --csv --colsep '	' -j 30 do_call {1} {2} {3} :::: "$WINDOWS_BED"

bcftools concat -a -Oz -o "$OUTDIR/All.dc.filtered.vcf.gz" --threads 32 "$OUTDIR"/*.filtered.vcf.gz
bcftools index "$OUTDIR/All.dc.filtered.vcf.gz"
