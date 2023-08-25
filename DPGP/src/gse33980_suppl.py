import pandas as pd


# Excel file clusters from supplementary material
# (sheet numer, cluster name)
clusters = [
    (2, "cluster1a"),
    (3, "cluster1b"),
    (4, "cluster2a"),
    (5, "cluster2b_early"),
    (6, "cluster2b_late"),
    (7, "cluster3"),
    (8, "cluster_dynamic"),
    (9, "cluster4a"),
    (10, "cluster4b"),
]

df_clusters = pd.DataFrame()

for sheet_num, cluster_name in clusters:
    # read sheet in Excel file
    df = pd.read_excel(
        "DPGP/papers/GSE33980 supplementary material/1471-2164-13-351-S5.xls", 
        sheet_num, 
        dtype={"genbank.gi": "object"}
    )

    # check if it has notes at the end of the table
    if any(df["ORF"].str.contains("*", regex=False)):
        # cleanup
        df = df[:-2]
        df["ORF"] = df["ORF"].str.replace("*", "", regex=False)
    
    # add column for cluster
    df.insert(1, "cluster", [cluster_name] * len(df))

    # merge dataframes
    df_clusters = pd.concat([df_clusters, df], ignore_index=True)

# save merged dataframe
csv_path = "DPGP/data/GSE33980/GSE33980_DGE_h2o2_suppl_clusters.csv"
df_clusters.to_csv(csv_path, index=False)

# save list of DGE as text
dge_path = "DPGP/data/GSE33980/GSE33980_DGE_h2o2_suppl.txt"
with open(dge_path, "w") as f:
    f.write(df_clusters["ORF"].to_string(index=False))

# save processed expression data from supplementary material
exp_csv_path = "DPGP/data/GSE33980/GSE33980_data_h2o2_suppl.csv"
df = pd.read_excel("DPGP/papers/GSE33980 supplementary material/1471-2164-13-351-S5.xls", 1)
df.to_csv(exp_csv_path, index=False)
