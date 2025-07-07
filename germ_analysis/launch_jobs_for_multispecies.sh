#!/usr/bin/env bash
RUN_PATH=$1
cd $RUN_PATH
for file in $(ls $RUN_PATH/*.fa)
do
    SAMPLE=`basename $file`
	echo $SAMPLE
    sbatch -N 1 --mem=32GB -c8 -J $SAMPLE -t 2:00:00 -o logs/entropy.%A.log --wrap="Rscript --vanilla /nemo/lab/ulej/home/users/farawar/GASR/germs/code/multispecies/calculate_entropy.r ${SAMPLE}"
done