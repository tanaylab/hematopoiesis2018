

source("./tests/testthat//workflow.R")

test_all = function(label="base", # what are we testing now? will be part of the output pathe, should be a valid sub-directory name (no "/").
                    datasets = c("e9",  "hsc_tier3", "melanoma"), #all contain a matrix called umis
                    datadir = "/net/mraid14/export/data/users/atanay/proj/scmat/datasets/",
                    outdir = "/net/mraid14/export/data/users/atanay/proj/scmat/tests/"){
  for(ds in datasets) {
    one_test(ds, label, datadir, outdir)
  }
}

one_test = function(ds, label,
                    datadir = "/net/mraid14/export/data/users/atanay/proj/scmat/datasets/",
                    outdir = "/net/mraid14/export/data/users/atanay/proj/scmat/tests/", nplots=100) {
  load(paste(datadir, ds, ".Rda", sep=""))
  call_workflow_tests(umis,label = label, basedir = paste(outdir, ds, sep="/"), nplots = nplots)
}
