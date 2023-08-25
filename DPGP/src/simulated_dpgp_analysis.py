import pandas as pd

def is_function(n: str) -> str:
    n = int(n)
    if n <= 200:
        simulation_function = "constant"
    elif n <= 400:
        simulation_function = "impulse"
    elif n <= 600:
        simulation_function = "linear"
    else:
        simulation_function = "unknown"
    return simulation_function


results_path_prefix = "DPGP/results/simulated/const_impulse_linear/"

for DE_analysis in ["", "DE_"]:
    for shape in [12, 9, 6, 3]:
        results_path = f"{results_path_prefix}{DE_analysis}shape_{shape}/"
        clusters = pd.read_csv(f"{results_path}case_optimal_clustering.txt", sep="\t")
        
        # count gene from each function
        clusters["function"] = clusters["gene"].str.extract("(\d+)", expand=False).apply(is_function)
        table = clusters.pivot_table(
            index='function', columns='cluster',
            aggfunc='count', values='gene',
        ).fillna(0).astype(int)

        # save table
        table.to_csv(f"{results_path}genes_function_clustered.csv")

        # export to LaTeX
        table.style.background_gradient(
            axis=1,
            cmap="PuBu"
        ).set_table_styles(
            [
                # horizontal lines
                {'selector': 'midrule', 'props': ':hline;'},
            ], overwrite=False
        ).format_index(
            "\\textbf{{{}}}", escape="latex", axis=1
        ).format_index(
            "\\textbf{{{}}}", escape="latex", axis=0
        ).format(escape="latex").to_latex(
            f"report/src/cil_{DE_analysis}{shape}.tex",
            convert_css=True,
            column_format=f"r|{len(table.columns)*'r'}"
        )
