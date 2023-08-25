import pandas as pd
from dpgp_preprocess import ExpressionMatrix


for dataset in ["linear", "const_impulse_linear"]:
    # read data
    data_path = f"DPGP/data/simulated/{dataset}/matObservedCounts.csv"
    df_data = pd.read_csv(data_path, index_col=0)
    df_data.insert(0, "gene", df_data.index)

    # filter for DE genes as per ImpulseDE2 results
    results_path = f"DPGP/data/simulated/{dataset}/dfImpulseDE2Results.csv"
    df_results = pd.read_csv(results_path)
    DE_genes = df_results[df_results["padj"] < 0.05]["Gene"]
    df_dataDE = df_data[df_data["gene"].isin(DE_genes)]

    # read metadata
    metadata_path = f"DPGP/data/simulated/{dataset}/dfAnnotation.csv"
    df_metadata = pd.read_csv(metadata_path)
    df_metadata["replicate"] = df_metadata["Sample"].str.extract(r"([rR]ep\d)")

    data = ExpressionMatrix(df_data, df_metadata)
    dataDE = ExpressionMatrix(df_dataDE, df_metadata)

    for data_obj, prefix in zip([data, dataDE], ["", "DE_"]):
        # split into replicates
        for rep in data_obj.metadata["replicate"].unique():
            data_rep = data_obj.filter_by_metadata_value("replicate", rep)
            data_rep.convert_labels_to_time(label_column="Sample", time_column="Time")

            out_path = f"DPGP/data/simulated/{dataset}/{prefix}expression_counts_{rep}"
            data_rep.save_expression(out_path)
