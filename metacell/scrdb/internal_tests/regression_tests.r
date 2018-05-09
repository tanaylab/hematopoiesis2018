library(testthat)
library(scrdb)

#
# run this function to save the output before making changes to the code.
# This function saves objects to file, so that after changing the code you can run the regression tests.
# Objects are witten to a sub-directory under basedir.
#
# return: the directory where the data was saved
.save_before_changing = function(basedir = "/net/mraid14/export/data/users/atanay/proj/scmat/regression_tests/",
                                label = "base",
                                umis = NULL){

  outdir=paste0(basedir,"/",label,"_",(format(Sys.time(), "%Y.%m.%d_%H.%M")))
  if (!dir.create(outdir)){
    stop()
  }
  message("Creating objects")

  objs = .create_obj(umis)
  message("writing output to ",outdir)

  # save to disk
  for(name in names(objs)) {
    tgscm_export(objs[name][[1]], file=paste0(outdir, "/", name))
  }
  return(outdir)
}

# "/net/mraid14/export/data/users/atanay/proj/scmat/regression_tests//human_bm_all_2017.04.06_14.19"
.do_regression = function(dir, umis=NULL){
                          # check_mat = T, check_cls = T, check_2d = T){
  new_objs = .create_obj(umis)
  old_objs = .import_old_objs(dir, new_objs)

  test_that("Number of objects before and after", {
    expect_equal(length(new_objs), length(old_objs))
  })

  cat("Doing regression, comparing to directory ", dir)
  context("Regression tests")

  for(i in 1:length(new_objs)) {
    test_that(names(new_objs)[i], {
      expect_equal(new_objs[[i]], old_objs[[i]])
    })
  }
  cat("\nDone\n")
}



# used before and after changing the code
.create_obj = function(umis = NULL, rseed= 1) {

  set.seed(rseed)
  if(is.null(umis)) {
    message("Using toydata")
    data("toydata")
    umis = toydata
  }

  # create objects
  scmat = tgScMat(mat=as.matrix(umis), stat_type="umi")
  marks = tgscm_gene_select(scmat)
  scl_knn = tgScMatClust(scmat=scmat, feat_mat=scmat@mat[marks,], alg_type="knn_graph")
  scl_kmeans = tgScMatClust(scmat=scmat, feat_mat=scmat@mat[marks,], alg_type="kmeans")
  scl_2d = tgScMatClustLayout2D(scl=scl_knn, alg_type="graph", T_edge=0.08, remove_confused_clusts=T)

  # matrix & markers are included in clustering solution

  return(list(scl_knn=scl_knn, scl_kmeans=scl_kmeans, scl_2d=scl_2d))
}

.import_old_objs = function(dir, sample_objs){
  onames = names(sample_objs)
  funcs = c("tgScMat"=import_mat,"tgScMatClust"=import_clust, "tgScMatClustLayout2D" = import_layout)
  old_objs = list()
  for(i in 1:length(sample_objs)) {
    func=funcs[class(sample_objs[[i]])][[1]]
    obj = func(file = paste0(dir, "/", onames[i]))
    old_objs[[i]] = obj
  }
  names(old_objs) = onames
  return(old_objs)
}



