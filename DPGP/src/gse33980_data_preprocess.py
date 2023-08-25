import pandas as pd
from dpgp_preprocess import ExpressionMatrix


#################################
#                               #
#   Read datasets               #
#                               #
#################################

geo_path = "DPGP/data/GSE33980/GSE33980_average_data.txt.gz"
suppl_path = "DPGP/data/GSE33980/GSE33980_data_h2o2_suppl.csv"

# read list of DEG
with open("DPGP/data/GSE33980/GSE33980_DGE_h2o2_suppl.txt", "r") as f:
    deg = f.read().splitlines()

# read datasets
df_geo = pd.read_csv(geo_path, sep="\t")
df_suppl = pd.read_csv(suppl_path)


#################################
#                               #
#   Get metadata                #
#                               #
#################################

## get metadata table for dataset dowloaded from GEO Datasets
# split column names but the gene column
metadata_geo = pd.Series(df_geo.columns[1:]).str.split("_", expand=True)

# rename columns
column_labels = {
    0: "treatment",
    1: "strain",
    2: "time",
    3: "_",
    4: "reference",
    5: "replicate",
}
metadata_geo = metadata_geo.rename(columns=column_labels)

# add full label column names
metadata_geo.insert(0, "label", pd.Series(df_geo.columns[1:]))

# clean up
metadata_geo = metadata_geo.drop(columns=["_"])

# format time column
metadata_geo["time"] = metadata_geo["time"].str.rstrip("min").replace("n", "-", regex=True).astype("int")

# save to file
metadata_geo_path = "DPGP/data/GSE33980/GSE33980_metadata.csv"
metadata_geo.to_csv(metadata_geo_path, index=False)


## get metadata table for dataset obtained from supplementary material
# split column names but the gene column
metadata_suppl = pd.Series(df_suppl.columns[1:]).str.split("_", expand=True)

# rename columns
column_labels = {
    0: "time",
    1: "strain",
    2: "replicate",
}
metadata_suppl = metadata_suppl.rename(columns=column_labels)

# format time column
metadata_suppl["time"] = metadata_suppl["time"].astype("int")

# add full label column names
metadata_suppl.insert(0, "label", pd.Series(df_suppl.columns[1:]))

# save to file
metadata_suppl_path = "DPGP/data/GSE33980/GSE33980_h2o2_suppl_metadata.csv"
metadata_suppl.to_csv(metadata_suppl_path, index=False)


#################################
#                               #
#   Filter datasets             #
#                               #
#################################

geo_data = ExpressionMatrix(df_geo, metadata_geo)
suppl_data = ExpressionMatrix(df_suppl, metadata_suppl)

# filter for only DEG in datasets
geo_data_deg = geo_data.filter_by_gene(deg)
suppl_data_deg = suppl_data.filter_by_gene(deg)

# keep only H2O2 tretment for GEO dataset
geo_data_deg_h2o2 = geo_data_deg.filter_by_metadata_value("treatment", "H2O2")

# save filtered data for clustering with DPGP (and relative metadata)
# tab delimited text files
geo_data_deg_h2o2_path = "DPGP/data/GSE33980/GSE33980_dpgp_geo"
suppl_data_deg_path = "DPGP/data/GSE33980/GSE33980_dpgp_suppl"

geo_data_deg_h2o2.save_expression(geo_data_deg_h2o2_path, metadata=True)
suppl_data_deg.save_expression(suppl_data_deg_path, metadata=True)


#################################
#                               #
#   Split datasets              #
#                               #
#################################

datasets = [
    [
        geo_data_deg_h2o2,
        geo_data_deg_h2o2_path
    ],
    [
        suppl_data_deg,
        suppl_data_deg_path
    ],
]

for data, path in datasets:
    # remove negative timepoints
    metadata_selection = data.metadata[data.metadata["time"] >= 0]

    # remove measurements that only have one replicate
    # metadata_selection = metadata_selection.groupby(["strain", "time"]).filter(lambda x: x['replicate'].nunique() > 1)

    data = data.filter_by_metadata(metadata_selection)

    # Divide each dataset into control and mutant strains for independent analysis
    for strain in data.metadata["strain"].unique():
        data_strain = data.filter_by_metadata_value("strain", strain)

        # Ulteriorly split into replicate
        for rep in data_strain.metadata["replicate"].unique():
            data_strain_rep = data_strain.filter_by_metadata_value("replicate", rep)

            # rename labels as timepoint values
            data_strain_rep.convert_labels_to_time()
            
            # save for clustering
            out_path = f"{path}_{strain}_{rep}"
            data_strain_rep.save_expression(out_path)
