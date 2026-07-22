#!/usr/bin/env bash

IN="/home/maycon/projects/turtles/work/plink2/LD_pruned.MN10.MX100.MS80.MQ20.admixture_input.bed"
K_MIN=1
K_MAX=10
N_RUNS=20
THREADS=8

OUT_ROOT=$(basename $IN .bed)

echo "Saving results to 'summary.csv'"
echo -e "k,run,ll,cv" > summary.csv

for each_k in $(seq $K_MIN $K_MAX); do
	mkdir K${each_k}
	for n_run in $(seq 1 $N_RUNS); do
		admixture -s time --cv $IN $each_k > out.log 
		
		ll=$(grep ^Loglikelihood out.log | cut -d: -f 2 | tr -d ' ')
		cv=$(grep ^CV out.log | cut -d: -f 2 | tr -d ' ')
		
		echo -e "${each_k},${n_run},${ll},${cv}" >> summary.csv

		mv $OUT_ROOT.${each_k}.P K${each_k}/$OUT_ROOT.${each_k}.${n_run}.P
		mv $OUT_ROOT.${each_k}.Q K${each_k}/$OUT_ROOT.${each_k}.${n_run}.Q
	done
done

rm out.log
