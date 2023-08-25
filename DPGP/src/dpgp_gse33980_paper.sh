#!/bin/bash

ura3="DPGP/data/GSE33980/GSE33980_dpgp_geo_[u]*[AB].txt"
d258="DPGP/data/GSE33980/GSE33980_dpgp_geo_[d]*[AB].txt"

awk -f DPGP/src/gse33980_replicates.awk ${d258} > DPGP/data/GSE33980/GSE33980_dpgp_geo_d258_C.txt
d258="DPGP/data/GSE33980/GSE33980_dpgp_geo_[d]*[AC].txt"

eval "$(conda shell.bash hook)"
conda activate dpgp

out_path="DPGP/results/GSE33980/paper"

mkdir -p ${out_path}

for EXPR_NAME in ura3 d258
do
    eval "EXPR=\$$EXPR_NAME"
    DP_GP_cluster.py -i ${EXPR} -o ${out_path}/geo_${EXPR_NAME} --fast --sigma_n2_shape 6 --sigma_n2_rate 2 -p png --plot
done

conda deactivate
