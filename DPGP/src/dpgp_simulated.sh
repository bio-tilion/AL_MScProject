#!/bin/bash

python DPGP/src/simulated_data_preprocess.py
mkdir -p DPGP/results/simulated

eval "$(conda shell.bash hook)"
conda activate dpgp

# each of the two simulated datasets
for folder in linear const_impulse_linear
do
    OUTDIR="DPGP/results/simulated/"${folder}
    mkdir -p ${OUTDIR}

    # on full dataset and DE only
    for pattern in "[^D]" "DE_"
    do
        EXPR=DPGP/data/simulated/${folder}/${pattern}*.txt
        if [ $pattern == DE_ ]
        then
            out_path=${OUTDIR}/DE_shape_
        else
            out_path=${OUTDIR}/shape_
        fi

        # with different shape values
        for shape in 12 9 6 3
        do
            mkdir -p ${out_path}${shape}

            # run DPGP
            DP_GP_cluster.py -i $EXPR -o ${out_path}${shape}/case --fast --true_times --sigma_n2_shape ${shape} -p png --plot
        done
    done
done

conda deactivate
