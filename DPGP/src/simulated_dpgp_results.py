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


def plot_cluster_gene_expression(clusters, gene_expression_matrix, t, t_labels, output_path_prefix, column_number=3):
    ''' 
    Plot gene expression over a time course with a panel for each cluster. Each panel contains
    transparent red lines for the expression of each individual gene within the cluster, the
    cluster mean, and a ribbon twice the standard deviation about the cluster mean.  This is
    essentially a wrapper function for GPy.plotting.matplot_dep.base_plots.gpplot.
    
    :param clusters: dictionary of dp_cluster objects
    :type clusters: dict
    :param gene_expression_matrix: expression over timecourse of dimension |genes|x|timepoints|
    :type gene_expression_matrix: pandas dataframe
    :param t: sampled timepoints
    :type t: numpy array of floats
    :param output_path_prefix: absolute path to output
    :type output_path_prefix: str
    
    :rtype: None (output is saved to file(s))
    
    ''' 
    # cluster IDs:
    IDs = sorted(clusters)
    # one panel per cluster:
    total_subplots = len(IDs)
    # max of 6 panels per figure or page
    subplots_per_fig = 12
    IDs_split = [IDs[i:i+subplots_per_fig] for i in xrange(0, len(IDs), subplots_per_fig)]
    index = 1
    for c, IDs in enumerate(IDs_split):
        total_cols = column_number # generate this many columns of subplots in the figure.
        total_rows = np.ceil(float(len(IDs))/float(total_cols)) # each figure generate will have this many rows.
        fig = plt.figure(num=None, figsize=(4*total_cols,4*total_rows), dpi=300, facecolor='w', edgecolor='k') #figsize=(12,8),
        for i, ID in enumerate(IDs):
            ax = fig.add_subplot(total_rows, total_cols, i+1)
            # create a range of values at which to evaluate the covariance function
            Xgrid = np.vstack(np.linspace(min(t), max(t), num=500))
            # calculate mean and variance at grid of x values
            mu, v = clusters[ID].model.predict(Xgrid, full_cov=False, kern=clusters[ID].model.kern)
            mu = np.hstack(mu.mean(axis=1))
            v = v[:,0]
            GPy.plotting.matplot_dep.base_plots.gpplot(Xgrid, mu, mu - 2*v**(0.5),  mu + 2*v**(0.5), ax=ax)
            ax.set_xlim((min(t),max(t)))
            ax.set_ylim((-3,3))
            
            # plot an x-axis at zero
            plt.axhline(0, color='black', ls='--', alpha=0.5)
            # plot the expression of each gene in the cluster
            for gene in list(clusters[ID].members):
                ax.plot(t, np.array(gene_expression_matrix.ix[gene]), color='red', alpha=0.1)
            
            # plot mean expression of cluster
            ax.plot(Xgrid, mu, color='blue')
            # create legend
            light_blue_patch = mpatches.Rectangle([0, 0], 1, 1, facecolor='#33CCFF', edgecolor='blue', lw=1, alpha=0.3)
            red_line = mlines.Line2D([], [], color='red', label='individual gene trajectory')
            ax.legend([ax.lines[0], light_blue_patch, red_line], \
                      ['cluster mean', u'cluster mean \u00B1 2 x std. dev.', 'individual gene trajectory'], 
                      loc=4, frameon=False, prop={'size':8})
            # prettify axes
            adjust_spines(ax, ['left', 'bottom'])
            # label x-axis
            ax.set_xlabel("Time")
            ax.set_xticks(t)
            ax.set_xticklabels(t_labels)
            ax.set_ylabel('Gene expression')
            ax.set_title('Cluster %s'%(index))
            index+=1
        
        file_out = output_path_prefix + '_new_expression_fig_' + str(c+1) + '.png'

        plt.tight_layout()
        plt.savefig(file_out)

        print "Expression figure saved to " + file_out


results_path_prefix = "DPGP/results/simulated/"

for simulated_dataset in ["linear", "const_impulse_linear"]:
    
    for DE_analysis in ["", "DE_"]:
        expression_matrices = [
            "DPGP/data/simulated/" + simulated_dataset + "/" + DE_analysis + "expression_counts_Rep1.txt",
            "DPGP/data/simulated/" + simulated_dataset + "/" + DE_analysis + "expression_counts_Rep2.txt",
            "DPGP/data/simulated/" + simulated_dataset + "/" + DE_analysis + "expression_counts_Rep3.txt",
        ]

        gene_expression_matrix, gene_names, t, t_labels = read_gene_expression(expression_matrices)

        for shape in [12, 9, 6, 3]:
            results_path = results_path_prefix + simulated_dataset + "/" + DE_analysis + "shape_" + str(shape) + "/case"

            clusterings = pd.read_csv(results_path + "_optimal_clustering.txt", sep="\t")
            if len(clusterings["cluster"].unique()) % 4 == 0 or len(clusterings["cluster"].unique()) > 9:
                n_cols = 4
            else:
                n_cols = 3
            
            optimal_clusters = get_optimal_clusters(gene_expression_matrix, t, results_path)
            plot_cluster_gene_expression(
                optimal_clusters, 
                pd.DataFrame(gene_expression_matrix, index=gene_names, columns=t),
                t, t_labels,
                results_path,
                column_number=n_cols
            )

print "Done."
