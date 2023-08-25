# DPGP - Dirichlet process Gaussian process mixture model
### Papers
- The full [DPGP paper](https://pubmed.ncbi.nlm.nih.gov/29337990/) is available as a [PDF file](papers/McDowell2018.pdf) (`McDowell2018.pdf`) in the `papers` folder. All supplementary materials are in the [corresponding subfolder](papers/supplementary%20material/).
- A second [paper](https://pubmed.ncbi.nlm.nih.gov/22846541/) is available as [PDF file](papers/GSE33980.pdf) (`GSE33980.pdf`). This is the original paper that produced the oxidative stress response dataset, which McDowell *et al.* used in their paper (figure 2). All supplementary materials are in the corresponding subfolder.

### Source repository
- The [code](code) file contains the github repository link.

### Datasets
- The text file [GEO_acc_list.txt](data/GEO_acc_list.txt) contains the list of all GEO accession records used in the original paper:
    - [GSE33980](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33980), oxidative stress response in archaebacterium *Halobacterium salinarum* exposed to H<sub>2</sub>O<sub>2</sub>
    - [GSE104714](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104714), glucocorticoid (dexamethasone) transcriptional response in human lung epithelial adenocarcinoma cell line A549

Running [geo_DPGP.py](src/geo_DPGP.py) will download those GEO Datasets in a corresponding subfolder in `data/`.

## DPGP installation
1. Create a conda environment with the following: ```$ conda create --name dpgp --file DPGP_conda_env.txt```
2. Activate the environment: ```$ conda activate dpgp```
3. Install DPGP: ```$ pip install git+https://github.com/PrincetonUniversity/DP_GP_cluster.git```
4. DPGP can be executed by calling: ```$ DP_GP_cluster.py -h``` for showing usage and additional arguments

### To Do:
- convert `conda` environment to `pipenv`
- upgrade DPGP code to python3 and latest libraries

## Usage
For a full description, see the [next section](#reproducing-results).

After installation, execute the following script:
```bash
$ ./data_gse33980.sh
```

Then run DPGP with the following commands:
```bash
$ conda activate dpgp
```
```bash
$ DP_GP_cluster.py -i DPGP/data/GSE33980/GSE33980_dpgp_geo_ura3_A.txt DPGP/data/GSE33980/GSE33980_dpgp_geo_ura3_B.txt -o DPGP/results/GSE33980/geo_ura3 --fast --sigma_n2_shape 6 --sigma_n2_rate 2 -p png --plot
```
```bash
$ DP_GP_cluster.py -i DPGP/data/GSE33980/GSE33980_dpgp_geo_d258_A.txt DPGP/data/GSE33980/GSE33980_dpgp_geo_d258_B.txt -o DPGP/results/GSE33980/geo_d258 --fast --sigma_n2_shape 6 --sigma_n2_rate 2 -p png --plot
```
```bash
$ conda deactivate
```

Or execute the following script:
```bash
$ ./dpgp_gse33980_paper.sh
```

## Reproducing results
### GSE33980 - oxidative stress response
GEO Datasets entry [GSE33980](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33980) contains two datasets listed as supplementary files:
1. Raw microarray data, `GSE33980_RAW.tar`
2. Processed microarray data, `GSE33980_average_data.txt.gz`

The latter was obtained as follows:
> Within the R Bioconductor m-array and limma packages, resultant raw data were background-subtracted using normexp, Loess normalized within each array, and quantile normalized between all arrays. Any of the 12 gene-specific probes for each gene lying outside the 99th % confidence interval were removed using Dixon’s test. Finally, remaining probe intensities for each gene were averaged and log2 ratios were calculated, yielding one expression ratio per gene.

In the original work, authors identified 626 differentially expressed genes (DEGs) in response to H<sub>2</sub>O<sub>2</sub>. Among these, 294 were also DEGs in response to the transcription factor rosR mutation, and were clustered with k-means clustering.

McDowell *et al.* ran DPGP on the 626 DEGs for both strains independently (parent and rosR mutated strains), obtaining similar results which 
> recapitulate previous knowledge of RosR-mediated gene regulation
in response to H2O2 with reduced user input.

#### Data collection
In order to gather and process the data for the clustering, I collected from the supplementary material the list of the 626 DEGs and the expression table that should be equivalent to the one uploaded on GEO Datasets.

Running [gse33980_suppl.py](src/gse33980_suppl.py) will extrapolate that information from the relative Excel table in the supplementary materials of the paper. The script will also save in the `data/GSE33980/` folder the files:
- `GSE33980_data_h2o2_suppl.csv`, complete microarray expression data table of strains treated with H<sub>2</sub>O<sub>2</sub>,
- `GSE33980_DGE_h2o2_suppl_clusters.csv`, single table containing information on all the clusters of DEGs from the original paper,
- `GSE33980_DGE_h2o2_suppl.txt`, text file containing only the list of the 626 DEGs.

#### Data preprocess
McDowell *et al.* independently ran DPGP on the two groups - parent/control and mutant strains - therefore, the datasets must be split first.

Running [gse33980_data_preprocess.py](src/gse33980_data_preprocess.py) will:
1. Read the datasets obtained in the previous step,
2. Generate a metadata table for each of those datasets to help in the filtering process,
3. Filter the datasets for the 626 DEGs,
4. Split each dataset into control and mutant expression tables for clustering.

The script will produce the following files in the `data/GSE33980/` folder:
- `GSE33980_metadata.csv`, metadata file for the GEO Dataset
- `GSE33980_h2o2_suppl_metadata.csv`, metadata file relative to the table obtained from the supplementary material
- `GSE33980_dpgp_[geo/suppl].txt`, full DEGs expression data in text format for DPGP clustering
- `GSE33980_dpgp_[geo/suppl]_metadata.csv`, metadata table relative to the DEG expression data above
- `GSE33980_dpgp_[geo/suppl]_[strain]_[replicate].txt`, DEGs expression data in text format split into control and mutant strains, and into each replicate for DPGP clustering. Negative timepoints have been excluded to mirror results in the paper. A 20min timepoint for the mutant strain has also been removed because it only has a single replicate (not two like all the others) and DPGP doesn't accept expression matrices of replicates with a different number of entries.

Note: once filtered and split, diffs of corresponding files reveal that the datasets originated from the GEO Dataset or the supplementary material are identical, therefore it is possible to focus on only one, *i.e.* the GEO Datasets one.

Running the script [data_gse33980.sh](src/data_gse33980.sh) 
```bash
$ ./data_gse33980.sh
```
will perform all these preliminary steps:
1. Download the data from GEO Dataset
2. Extrapolate information from the supplementary material
3. Data preprocesing

#### DPGP
To run DPGP, first activate the environment:
```bash
$ conda activate dpgp
```
Then run DPGP with the following parameters to mirror the analysis from the paper:
```bash
$ DP_GP_cluster.py -i EXPRESSION_MATRICES -o OUTPUT_FOLDER --fast --true_times --sigma_n2_shape 6 --sigma_n2_rate 2 -p png --plot
```
The argument `--fast` will run DPGP in *fast mode*, like it was run for the datasets in the paper.

The argument `--plot` will output images of clusters like the ones in the paper, in the format specified by `-p`.

The arguments `--sigma_n2_shape` and `--sigma_n2_rate` are the hyperparameters $\alpha$ and $\beta$ respectively, that can be varied to allow for greater variability by reducing the shape parameter. Because these are microarray data, McDowell *et al.* used $\alpha=6$ (default 12) and $\beta=2$ (default 2).

The argument `--true_times` assumes the headers of the expression matrix are numerical values corresponding to each timepoint, which do not need to be equally spaced. If not specified, the software assumes the timepoints are equally spaced.

Run the tool for each strain using all replicate expression matrices. Afterwards, deactivate the environment:
```bash
$ conda deactivate
```

Alternatively, DPGP can be run by simply executing the script [dpgp_gse33980_paper.sh](src/dpgp_gse33980_paper.sh)
```bash
$ ./dpgp_gse33980_paper.sh
```
which
1. deals with a missing time point in the second replicate od the mutant strain
2. activates the conda environment 
3. runs DPGP for both strains with the same parameters used in the paper
4. deactivates the conda environment
