#!/bin/bash

dynamics=(doublewell mutualistic SIS genereg wilsoncowan) # no SIS when cparam is u

for dynamic in ${dynamics[@]}; do
    jobname=fmri-TAPAS-$dynamic

    sbatch --job-name=$jobname --output=./output/${jobname}.out collect-TAPAS-fmri.sh $dynamic
done
