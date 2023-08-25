from DP_GP import core
from DP_GP import cluster_tools
from DP_GP.plot import adjust_spines
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from collections import defaultdict
import GPy

np.random.seed(1234)


def read_gene_expression(expression):
    gene_expression_matrix, gene_names, t, t_labels = core.read_gene_expression_matrices(
        expression
    )
    return gene_expression_matrix, gene_names, t, t_labels


def get_optimal_clusters(gene_expression_matrix, t, inputs_path_head):
    # clusterings
    clusterings_path = inputs_path_head + "_clusterings.txt"
    sampled_clusterings = pd.read_csv(clusterings_path, delim_whitespace=True)
    gene_names = list(sampled_clusterings.columns)

    # log likelihoods
    log_likelihoods_path = inputs_path_head + "_log_likelihoods.txt"
    with open(log_likelihoods_path, 'r') as f:
        log_likelihoods = [float(line.strip()) for line in f]

    # select best clustering by maximum a posteriori estimate
    optimal_clusters = cluster_tools.best_clustering_by_log_likelihood(
        np.array(sampled_clusterings), log_likelihoods
    )

    # combine gene_names and optimal_cluster info
    optimal_cluster_labels = defaultdict(list)
    optimal_cluster_labels_original_gene_names = defaultdict(list)
    for gene, (gene_name, cluster) in enumerate(zip(gene_names, optimal_clusters)):
        optimal_cluster_labels[cluster].append(gene)
        optimal_cluster_labels_original_gene_names[cluster].append(gene_name)

    # optimize GP model for best clustering
    optimal_clusters_GP = {}
    for cluster, genes in optimal_cluster_labels.iteritems():
        optimal_clusters_GP[cluster] = core.dp_cluster(
            members=genes, 
            X=np.vstack(t), 
            Y=np.array(np.mat(gene_expression_matrix[genes,:])).T
        )
        optimal_clusters_GP[cluster] = optimal_clusters_GP[cluster].update_cluster_attributes(gene_expression_matrix)

    return optimal_clusters_GP


def plot_cluster_expression(optimal_clusters_GP, gene_expression_matrix, gene_names, t, t_labels, output_path_head, strain):
    gene_expression_matrix = pd.DataFrame(gene_expression_matrix, index=gene_names, columns=t)
    # predict mean expression
    Xgrid = np.vstack(np.linspace(min(t), max(t), num=500))
    for id, cluster in optimal_clusters_GP.iteritems():
        mu, v = cluster.model.predict(Xgrid, full_cov=False, kern=cluster.model.kern)
        mu = np.hstack(mu.mean(axis=1))
        v = v[:,0]

        fig = plt.figure(num=None, figsize=(4,4), dpi=300, facecolor='w', edgecolor='k')
        ax = fig.add_subplot(111)
        GPy.plotting.matplot_dep.base_plots.gpplot(Xgrid, mu, mu - 2*v**(0.5),  mu + 2*v**(0.5), ax=ax)
        ax.set_xlim((min(t),max(t)))
        ax.set_ylim((-3,3))

        # plot an x-axis at zero
        plt.axhline(0, color='black', ls='--', alpha=0.5)
        # plot the expression of each gene in the cluster
        for gene in list(cluster.members):
            ax.plot(t, np.array(gene_expression_matrix.ix[gene]), color='red', alpha=0.1)

        # plot mean expression of cluster
        ax.plot(Xgrid, mu, color='blue')
        # create legend
        light_blue_patch = mpatches.Rectangle([0, 0], 1, 1, facecolor='#33CCFF', edgecolor='blue', lw=1, alpha=0.3)
        red_line = mlines.Line2D([], [], color='red', label='individual gene trajectory')
        ax.legend(
            [ax.lines[0], light_blue_patch, red_line],
            ['cluster mean', u'cluster mean \u00B1 2 x std. dev.', 'individual gene trajectory'],
            loc=4, frameon=False, prop={'size':8}
        )

        # prettify axes
        adjust_spines(ax, ['left', 'bottom'])
        # label axis
        ax.set_xlabel("Time")
        ax.set_xticks(t)
        ax.set_xticklabels(t_labels)
        ax.set_ylabel('Gene expression')
        ax.set_title('%s cluster %s (n=%d)'%(strain, id, len(cluster.members)))

        file_out = output_path_head + "_fig_" + str(id) + ".png"

        plt.tight_layout()
        plt.savefig(file_out)

        print "Cluster figure saved to " + file_out


def save_mean_expression(optimal_clusters_GP, output_path_head):
    # predict mean expression
    Xgrid = np.vstack(np.linspace(min(t), max(t), num=500))
    df_out = pd.DataFrame()
    for id, cluster in optimal_clusters_GP.iteritems():
        # calculate mean and variance at grid of x values
        mu, v = cluster.model.predict(Xgrid, full_cov=False, kern=cluster.model.kern)
        mu = pd.Series(np.hstack(mu))
        df_out[id] = mu

    try:
        df_out_path = output_path_head + "_mean_expression.csv"
        df_out.to_csv(df_out_path, index=False)
        txt_out = "Mean expression per cluster matrix saved to " + df_out_path
    except:
        txt_out = "Something went wrong"
    
    print txt_out


def plot_similarity_matrix(similarity_matrix, gene_index, output_path_head, strain):
    similarity_matrix = pd.read_csv(similarity_matrix, sep="\t", index_col=0)
    with open(gene_index, "r") as f:
        gene_index = f.read().splitlines()

    similarity_matrix = similarity_matrix.loc[gene_index,:].loc[:,gene_index]

    file_out = output_path_head + "_heatmap.png"

    fig = plt.figure(figsize=(5,5))
    ax = fig.add_axes([0.05, 0.1, 0.75, 0.75])

    heatmap = ax.matshow(similarity_matrix, aspect='auto', origin='lower', cmap="cubehelix", vmax=1, vmin=0)
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])

    plt.subplots_adjust(bottom=0.1, left=0.1, right=0.9, top=0.9)
    ax_bar = fig.add_axes([0.82, 0.1, 0.02, 0.75])
    plt.colorbar(heatmap, cax=ax_bar)
    ax_bar.set_ylabel(
        'Proportion of Gibbs samples in which \nrow i and column j were co-clustered',
        rotation=270, labelpad=25
    )

    ax.set_title("%s strain"%(strain), fontweight = "bold", fontsize="xx-large")
    plt.savefig(file_out)

    print "Heatmap figure saved to " + file_out


ura3_expression = [
    "DPGP/data/GSE33980/GSE33980_dpgp_geo_ura3_A.txt",
    "DPGP/data/GSE33980/GSE33980_dpgp_geo_ura3_B.txt"
]
d258_expression = [
    "DPGP/data/GSE33980/GSE33980_dpgp_geo_d258_A.txt",
    "DPGP/data/GSE33980/GSE33980_dpgp_geo_d258_C.txt"
]

results_path = "DPGP/results/GSE33980/paper/"

for expression, prefix, strain in zip([ura3_expression, d258_expression], ["geo_ura3", "geo_d258"], ["Control", "Mutant"]):
    gene_expression_matrix, gene_names, t, t_labels = read_gene_expression(expression)
    optimal_clusters = get_optimal_clusters(gene_expression_matrix, t, results_path + prefix)
    save_mean_expression(optimal_clusters, results_path + prefix)
    plot_cluster_expression(
        optimal_clusters, gene_expression_matrix,
        gene_names, t, t_labels,
        results_path + prefix, strain
    )

ura3_similarity_matrix = "DPGP/results/GSE33980/paper/geo_ura3_posterior_similarity_matrix.txt"
d258_similarity_matrix = "DPGP/results/GSE33980/paper/geo_d258_posterior_similarity_matrix.txt"
ura3_gene_index = "DPGP/results/GSE33980/paper/geo_ura3_posterior_similarity_matrix_heatmap_key.txt"
d258_gene_index = "DPGP/results/GSE33980/paper/geo_d258_posterior_similarity_matrix_heatmap_key.txt"

for similarity_matrix, gene_index, prefix, strain in [
    (ura3_similarity_matrix, ura3_gene_index, "geo_ura3", "Control"),
    (d258_similarity_matrix, d258_gene_index, "geo_d258", "Mutant"),
    (d258_similarity_matrix, ura3_gene_index, "geo_d258_alt", "Mutant"),
]:
    plot_similarity_matrix(similarity_matrix, gene_index, results_path + prefix, strain)
