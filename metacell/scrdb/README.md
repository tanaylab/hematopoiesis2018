### Brief overview of the package for programmers ###


####Notations ####

_i_ - cell

_g_ - gene

_c_ - cluster

_b_ - batch 

_u<sub>ig</sub>_  - the number of umi's in cell i gene g


#### scmat.r ####

This class is wrapping a sparse matrix of UMI with meta data (batch information usually). It mainly provides support for reading MARS-seq batches and pooling them together.

The class also provide support for very basic contamination filtering: given a user defined epsilon (\(\epsilon\)) parameter, it find for a gene _g_ the maximal threshold _T_ such that _\(\epsilon\) * N<sub>g</sub> > sum(ifelse(u<sub>ig</sub> <= T, u<sub>ig</sub>, 0))_ and clean umis of gene _g_ in cells with less than or equal to _T_ umis.

#### gene_select.r####

This is supporting the computation of multiple types of statistics on the umi matrix, giving each gene a set of scores that can be used to select markers for subsequent clustering.

#### cluster.r####

This class wraps a matrix object and support the generation of cell clustering using several methods (currently knn_graph clustering, kmeans). Clustering relies on getting a list of features from the user - it is not finding the features itself.

The cluster object also support the computation of post-process statistics, most importantly, _\@clust_fp_ which define for each cluster and gene the fold change enrichment of the gene compared to the median over all clusters, _\@gene_cov_ which tells us how many of the cells in each cluster have at least one umi per gene.

Using a clustering solution, it is also possible to apply two cleaning procedures on the data. First, there is a function for removing ambient noise (or cross-well contamination) by searching for clusters that express a gene at a level that is smaller than _\(\epsilon\) * N<sub>g</sub>_. Second, we can detect cells that are outliers by searching for genes that are experssed in a certain cell at levels that are much higher than expected by their cluster. Suspects outliers are being identified and visualized for further downstream analysis.

####knn_graph.r####

This is our current knn clustering algorithm.

####cluster_layout.r####

####pipeline.r####

Contains the main building blocks for scRNA analysis. These functions use an scdb object (a simple cache mechanism) to store and load their output. 

Here's a typical R script that performs a full cycle of analysis:

```
#!r

# If you have a copy of the scrdb repository then you'll want to load the scrdb library 
require(devtools)
load_all(Sys.getenv("SCRDB_HOME"))

# Otherwise, if you just installed scrdb, it is enough to require(scrdb) 

# parameter management package.
require(tgconfig)

# Override default parameters with your project specfic ones
override_params("/net/mraid14/export/data/users/atanay/proj/devscrdb/workdir/melanoma_cd45_cd3/test_new_pipe/test_mel_cd45_cd3_params.yaml", "scrdb")

# init the object cache manager
scdb = scdb_init()

# load MARS-seq batches and aggregate them to a single umi matrix
sc_raw_mat = sc_pipe_build_mat(scdb)

# clean raw umi matrix
sc_clean_mat = sc_pipe_mat_to_clean_mat(scdb, sc_raw_mat)

# generate cell modules from the cleaned umi matrix
sc_cl = sc_pipe_mat_to_cmods(scdb, sc_clean_mat)

# generate gene modules from cell modules 
sc_gm_fp = sc_pipe_cmods_to_gmods(scdb, sc_cl, feat_type = "clust_fp")
sc_gm_resid = sc_pipe_cmods_to_gmods(scdb, sc_cl, feat_type = "resid")

# generate plots
sc_2d = sc_pipe_plots(scdb, sc_cl)

```

