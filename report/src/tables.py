import pandas as pd


# Table adapted from Oh2021 paper
table = pd.read_csv("report/src/oh2021.csv", sep=";")

table.style \
.hide(axis=0) \
.format_index("\\textbf{{{}}}", escape="latex", axis=1) \
.format(escape="latex") \
.to_latex(buf="report/src/oh2021.tex")

# Table obtained with entrez.py script
types = {
    'PMID': 'object',
    'geo_dataset_id': 'object',
    'bioproject_id': 'object'
}
table = pd.read_csv("report/src/oh2021_metric.csv", sep=",", dtype=types)

# Filter table for readability: subset of columns and relative new labels
subset_cols = [
    "Name", "Time course",
    "cited_by", "avg_year_cit",
    "geo_dataset_id", "bioproject_id"
]

relabels = [
    "Name", "Time course",
    "Cited by", "Average yearly citations",
    "GEO Dataset ID", "Bioproject ID"
]

col_format_latex = [
    "l", "l",
    "c", ">{\centering}m{2.3cm}",
    "c", "c"
]

new_labels = dict()
for i, j in zip(subset_cols, relabels):
    new_labels[i] = j

col_latex = "".join(col_format_latex)

# Sort table before export it in LaTeX
table[subset_cols].sort_values("avg_year_cit", ascending=False).rename(new_labels, axis=1).style \
.hide(axis=0) \
.relabel_index(relabels, axis=1) \
.format_index("\\textbf{{{}}}", escape="latex", axis=1) \
.format(escape="latex", na_rep="", precision=2) \
.to_latex(buf="report/src/oh2021_metric.tex", column_format=col_latex)
