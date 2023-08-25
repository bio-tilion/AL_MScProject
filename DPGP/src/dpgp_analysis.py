import pandas as pd
from scipy.stats import fisher_exact

def cocluster(df_one_grouped, df_two_grouped) -> pd.DataFrame:
    co_cluster_dict = dict()
    for one_name, one_group in df_one_grouped:
        co_cluster_dict[one_name] = dict()
        for two_name, two_group in df_two_grouped:
            co_cluster_dict[one_name][two_name] = one_group["gene"].isin(two_group["gene"]).sum()

    return pd.DataFrame(co_cluster_dict)

def get_contingency_tables(df_cocluster: pd.DataFrame, corresponding_cluster: list[tuple]) -> dict[tuple,list]:
    n = df_cocluster.values.sum()
    contingency_dict = dict()

    for ctrl_cluster, mut_cluster in corresponding_cluster:
        ctrl = df_cocluster.loc[ctrl_cluster, :]
        mut = df_cocluster.loc[:, mut_cluster]
        in_ctrl_in_mut = ctrl.pop(mut_cluster)
        in_ctrl_no_mut = ctrl.sum()
        no_ctrl_in_mut = mut.drop(ctrl_cluster).sum()
        no_ctrl_no_mut = n - (in_ctrl_in_mut + in_ctrl_no_mut + no_ctrl_in_mut)

        contingency_dict[(ctrl_cluster, mut_cluster)] = [
            [in_ctrl_in_mut, in_ctrl_no_mut],
            [no_ctrl_in_mut, no_ctrl_no_mut]
        ]
    
    return contingency_dict

def get_clusters_fishers_exact_test(contingency_dict: dict[tuple,list]) -> dict[tuple,tuple]:
    clusters_fet = dict()

    for corresponding_clusters, cont_table in contingency_dict.items():
        fet = fisher_exact(cont_table)
        clusters_fet[corresponding_clusters] = fet

    return clusters_fet
