import pandas as pd
import numpy as np
from dpgp_analysis import cocluster
from dpgp_analysis import get_contingency_tables
from dpgp_analysis import get_clusters_fishers_exact_test
from scipy.stats import pearsonr


mean_expr_ura3 = pd.read_csv("DPGP/results/GSE33980/paper/geo_ura3_mean_expression.csv")
mean_expr_d258 = pd.read_csv("DPGP/results/GSE33980/paper/geo_d258_mean_expression.csv")

# calculate pearson correlation
corr_matrix = list()
for ctr_cluster in mean_expr_ura3.columns:
    temp = list()
    x = mean_expr_ura3[ctr_cluster]
    for mut_cluster in mean_expr_d258.columns:
        y = mean_expr_d258[mut_cluster]
        r, p = pearsonr(x, y)
        temp.append((r, p))
    corr_matrix.append(temp)

# convert to numpy matrix and save
corr_matrix = np.array(corr_matrix)
np.save("DPGP/results/GSE33980/paper/geo_correlation_matrix.npy", corr_matrix)

# split matrix into two: r and p value
corr_df_r = pd.DataFrame(np.split(corr_matrix, 2, axis=2)[0].reshape((6,6)))
corr_df_p = pd.DataFrame(np.split(corr_matrix, 2, axis=2)[1].reshape((6,6)))

# relabel
idx_labels = {i: f"ctrl cluster {i+1}" for i in range(6)}
col_labels = {i: f"mut cluster {i+1}" for i in range(6)}
corr_df_r.rename(columns=col_labels, index=idx_labels, inplace=True)
corr_df_p.rename(columns=col_labels, index=idx_labels, inplace=True)

# export to latex
corr_df_r.style.set_table_styles(
    [
        # horizontal lines
        {'selector': 'midrule', 'props': ':hline;'},
    ], overwrite=False
).highlight_max(
    axis=1, props='cellcolor:{lightgray}; bfseries: ;'
).format(
    escape="latex", precision=4
).to_latex(
    "report/src/pearson_r.tex",
    column_format="l|" + ("r" * 6)
)

# equivalent clusters
eq_clusters = [(ctrl, mut) for ctrl, mut in corr_df_r.idxmax(axis=1).items()]


clusters_control = pd.read_csv("DPGP/results/GSE33980/paper/geo_ura3_optimal_clustering.txt", sep="\t")
clusters_mutant = pd.read_csv("DPGP/results/GSE33980/paper/geo_d258_optimal_clustering.txt", sep="\t")

grouped_clusters_control = clusters_control.groupby("cluster")
grouped_clusters_mutant = clusters_mutant.groupby("cluster")

# get table of co-clustered genes
control_vs_mutant = cocluster(grouped_clusters_mutant, grouped_clusters_control)

# relabel
idx_labels = {i: f"ctrl cluster {i}" for i in control_vs_mutant.index}
col_labels = {i: f"mut cluster {i}" for i in control_vs_mutant.columns}
control_vs_mutant.rename(columns=col_labels, index=idx_labels, inplace=True)

# restyle equivalent clusters cells
style_df = pd.DataFrame(
    "", 
    index=[f"ctrl cluster {i}" for i in range(1,7)],
    columns=[f"mut cluster {i}" for i in range(1,7)],
)
for row, col in eq_clusters:
    style_df.loc[row, col] = 'cellcolor:{lightgray}; bfseries: ;'

# export to latex
control_vs_mutant.style.apply(
    lambda _: style_df, axis=None
).set_table_styles(
    [
        # horizontal lines
        {'selector': 'midrule', 'props': ':hline;'},
    ], overwrite=False
).format(
    escape="latex"
).to_latex(
    "report/src/cocluster.tex", 
    column_format="l|" + ("r" * 6)
)


# calculate contingency tables and Fisher's exact test
contingency_tables = get_contingency_tables(control_vs_mutant, eq_clusters)
fet = get_clusters_fishers_exact_test(contingency_tables)

col_labels = [
    "control_strain", "mutant_strain",
    "pearson_r", "pearson_p",
    "fet_odds_ratio", "fet_p",
    "diff_dyn",
]

# create summary dataframe and save
summary_df = pd.DataFrame(
    [
        (
            ctrl, mut, 
            corr_df_r.loc[ctrl, mut], corr_df_p.loc[ctrl, mut], 
            fet[(ctrl, mut)][0], fet[(ctrl, mut)][1],
            f"{contingency_tables[(ctrl, mut)][0][1]:3} ({contingency_tables[(ctrl, mut)][0][1]/sum(contingency_tables[(ctrl, mut)][0]):.1%})",
        ) for ctrl, mut in eq_clusters
    ],
    columns=col_labels
)
summary_df.to_csv("DPGP/results/GSE33980/paper/equivalent_clusters.csv", index=False)

# workaround for latex column formatting
summary_df["temp"] = pd.Series(None, dtype="object")

# export to latex
relabels = [
    "Control strain", "Mutant strain",
    "Pearson corr r", "Pearson p value",
    "Odds ratio", "FET p value",
    "Different dynamic genes", ""
]

summary_df.rename(
    columns={old: new for old, new in zip(summary_df.columns, relabels)}
).style.hide(
    axis=0
).format(
    # formatter for numerical values
    {
        "Pearson corr r": "{:.4f}",
        "Pearson p value": "{:.2e}",
        "Odds ratio": "{:.2f}",
        "FET p value": "{:.2e}"
    },
    # workaround: empty last column
    na_rep = "",
    escape="latex",
).format_index(
    "\\textbf{{{}}}", escape="latex", axis=1
).to_latex(
    "report/src/summary.tex", 
    column_format = ">{\centering}m{1.7cm}" +
        ">{\centering}m{1.8cm}" +
        ">{\centering}m{1.2cm}S[table-format = 1.2e2]" +
        ">{\centering}m{1cm}S[table-format = 1.2e2]" +
        ">{\centering}m{2.2cm}p{0pt}"
)


# print a summary on screen
tot_diff_gene = int()

for idx, col in sorted(eq_clusters, key=lambda x: x[0]):
    pearson_r = corr_df_r.loc[idx, col]
    pearson_p = corr_df_p.loc[idx, col]
    fishers_p = fet[(idx, col)][1]
    diff_gene = contingency_tables[(idx, col)][0][1]
    tot_ctrl_gene = sum(contingency_tables[(idx, col)][0])
    tot_diff_gene += diff_gene
    print(
        idx, " -> ", col, "\t",
        f"(r={pearson_r:.4}, p={pearson_p:.2e})\n",
        "\t" * 5,
        f"Fisher's exact test: p={fishers_p:.2e}",
        "-> significant relationship\n" if fishers_p < 0.05 else "\n",
        "\t" * 5,
        f"genes with different dynamic: {diff_gene:3} ({diff_gene/tot_ctrl_gene:.2%})"
    )

print(f"Total number of genes with different dynamic: {tot_diff_gene}")

# cluster 5
print("\nCluster 5: upregulation after 40'")
tot_cluster5_genes = control_vs_mutant.loc['ctrl cluster 5', :].sum()
inv_dynamic_genes = control_vs_mutant.loc["ctrl cluster 5",:][corr_df_r.loc["ctrl cluster 5",:] < -0.5]
tot_inv_dynamic_genes = inv_dynamic_genes.sum()

print(
    "Total number of genes in control cluster:".rjust(63),
    f"{tot_cluster5_genes:4}\n",
    "Genes with inverted dynamic in corresponding mutant clusters:".rjust(62),
    f"{tot_inv_dynamic_genes:4}", f"({tot_inv_dynamic_genes/tot_cluster5_genes:.2%})"
)
for i, (cluster, n) in enumerate(inv_dynamic_genes.items()):
    txt = "of which" if i==0 else ""
    print(
        txt.rjust(62),
        f"{n:5} in {cluster}",
    )

print(
    "Paper results", "\n"
    "Total number of genes with different dynamic: 372", "\n"
    "Cluster upregulated after 40'", "\n",
    "Total number of genes in control cluster:  141".rjust(67), "\n",
    "Genes with inverted dynamic in corresponding mutant clusters:   89".rjust(67), f"({89/141:.2%})", "\n",
    "of which   72".rjust(67), "in mutant cluster 3", "\n",
    "17".rjust(67), "in mutant cluster 5"
)
