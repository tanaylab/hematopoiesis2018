# used by test_workflow.R but can also work on other datasets

t_scmat <- function(umis, cell_metadata=NULL) {

  if (is.null(cell_metadata)){
    cell_metadata = data.frame(a=round(1:ncol(umis), -2), b = (1:ncol(umis)) %% 10, row.names = colnames(umis))
  }

  context("tgScMat")
  # we first initialize a tgScMat object
  # (umis was added to data/umis.Rda )

  scmat = tgScMat(mat=as.matrix(umis), stat_type="umi", cell_metadata=cell_metadata)

  test_that("scmat created", {
    expect_false(is.null(scmat))
    expect_equal(dim(scmat@mat), dim(umis))
  })

  return(scmat)
}

t_feature_select <- function(scmat, outdir, n_marks) {

  #we next compute gene markers, while also dumping stats into a fig n gset.png
  #this stage is likley requiring intervention by tuning parameters
  context("feature selection")


  marks = tgscm_gene_select(scmat,n_marks = n_marks, fname="gene_set.png", outdir = outdir)

  test_that("markers selscted with tgscm_gene_select", {
    expect_gt(length(marks), 1)
    expect_lte(length(marks), n_marks)
  })
  return(marks)
}

# basic sanity tests for 1 clustering result
t_one_clustering <- function(desc, scl, scmat) {
  test_that(desc, {
    expect_equal(scl@scmat, scmat)
    expect_equal(length(scl@clusts), length(scmat@cells))
    expect_equal(ncol(scl@clust_fp), length(unique(scl@clusts)))
  })
}

# checks different clustering methods, returns a list of clustering results
t_clustering <- function(scmat, marks) {

  context("clustering")

  scl = tgScMatClust(scmat=scmat, feat_mat=scmat@mat[marks,], alg_type="knn_graph")
  t_one_clustering("check knn_graph clustering", scl, scmat)

  scl_kmeans <- tgScMatClust(scmat=scmat, feat_mat=scmat@mat[marks,], alg_type="kmeans")
  t_one_clustering("check kmeans clustering", scl_kmeans, scmat)

  cat("done testing clustering")

  return (list(kmeans = scl_kmeans, knn_graph = scl))
}

t_export_import <-function(scl_2d, outdir) {
  context("serialization (export/import)")

  # TODO: check partial export

  fname = paste0(outdir,"/data")

  scl = scl_2d@scl
  scmat = scl@scmat


  test_that("matrix serialization (export/import)", {
    tgscm_export(scmat, file=fname)
    imp_scmat = import_mat(fname)

    expect_equal(imp_scmat, scmat)
  })

  test_that("cluster serialization (export/import)", {
    # export all clustering
    fname2 = paste0(outdir,"/data2")
    tgscm_export(scl, file=fname2)
    imp_scl = import_clust(fname2)

    expect_equal(imp_scl, scl)
  })

  test_that("cluster serialization (export/import) when not exporting scmat", {
    # export without matrix clustering
    fname3 = paste0(outdir,"/data3")
    tgscm_export(scl, file=fname3, supress_mat=T)
    imp_scl3 = import_clust(fname3, scmat = scmat)

    expect_equal(imp_scl3, scl)
  })

  test_that("cluster layout serialization (export/import)", {
    # skip("graph_layout currently not saved to disk!")
    # export all clustering
    fname4 = paste0(outdir,"/data4")
    tgscm_export(scl_2d, file=fname4)
    imp_scl2d = import_layout(fname4)

    expect_equal(imp_scl2d, scl_2d)
  })
  print("done export/import")
}

# test visualization functions
t_visualize <- function(scl, outdir, marks, nplots = 5) {

  context("visualization")
  print("testing visualization")

  #we then compute the graph/cells 2d layout
  scl_2d = tgScMatClustLayout2D(scl=scl, alg_type="graph", T_edge=0.08)

  test_that("2D projection of the clusters was created", {
    ncells = length(scl@clusts)
    # expect_equal(length(scl_2d@x_cl), scl@nclust)  # currently not true....
    # expect_equal(length(scl_2d@y_cl), scl@nclust)
    # expect_equal(length(scl_2d@x), ncells)
    # expect_equal(length(scl_2d@y), ncells)

    expect_true(all(is.finite(scl_2d@x_cl)))
    expect_true(all(is.finite(scl_2d@y_cl)))
    expect_true(all(is.finite(scl_2d@x)))
    expect_true(all(is.finite(scl_2d@y)))
  })

  #plotting the clusters and cells around them - with some default coloring scheme
  tgscm_plot_clust_2d(scl_2d, fname = paste(outdir, "clusters.png",sep="/"))

  tgscm_plot_mat(scl, png_fn = paste(outdir, "heatmap.png",sep="/"))

  #we can now plot the density map for marker genes
  mk = head(marks, nplots)
for(nm in mk) {
    tgscm_plot_gene_2d(scl_2d,gene_nm=nm, reg_factor=5, w=2000,h=2000,base_dir=outdir)
  }

  return (scl_2d)
}

######################################################################################
# calls the testing functions, which use testthat package.
######################################################################################
call_workflow_tests <- function(umis, label = "toydata", basedir = "/tmp/", nplots = 5, n_marks = 600,
                                cell_metadata = NULL) {
  set.seed(1234)
  outdir=paste0(basedir,"/",label,"_",(format(Sys.time(), "%Y.%m.%d_%H.%M")))
  if (!dir.create(outdir)){
    stop()
  }
  message("writing output to ",outdir)

  scmat = t_scmat(umis, cell_metadata = cell_metadata )

  marks = t_feature_select(scmat, outdir,n_marks)

  # checks different clustering methods, returns a list of clustering results
  scls = t_clustering(scmat, marks)
  scl = scls[[1]]

  for (scl in scls) {
    scldir = paste0(outdir, "/", scl@alg_params$type)
    dir.create(scldir)
    scl_2d = t_visualize(scl, scldir, marks, nplots)
  }

  # export objects to files and read again
  t_export_import(scl_2d, outdir)

}
