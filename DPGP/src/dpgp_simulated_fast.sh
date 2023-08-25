#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate dpgp

EXPR="DPGP/data/simulated/*.txt"
out_path="DPGP/results/simulated/fast_shape"

for shape in $(seq 12 -2 2)
do
    mkdir -p ${out_path}${shape}
    DP_GP_cluster.py -i $EXPR -o ${out_path}${shape}/case --fast --sigma_n2_shape ${shape} -p png --plot
done

conda deactivate
