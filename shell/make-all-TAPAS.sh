#!/bin/bash

datadir=/projects/academic/naokimas/neil/brains-ns50/

for file in ${datadir}*; do
    infile=\'${file}\' # because the glob already gives the full path
    jobname=${file/${datadir}/""} # just want the specific file name

    # echo $infile
    # echo $jobname
    # break
    sbatch --job-name=${jobname} --output=./output/TAPAS/${jobname}.out request-TAPAS-fmri.sh ${infile}
done
