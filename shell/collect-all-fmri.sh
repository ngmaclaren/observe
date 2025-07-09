#!/bin/bash

dynamics=(doublewell mutualistic SIS genereg wilsoncowan)

for dynamic in ${dynamics[@]}; do
    jobname=fmri-$dynamic

    sbatch --job-name=$jobname --output=./output/${jobname}.out collect-fmri.sh $dynamic
done
