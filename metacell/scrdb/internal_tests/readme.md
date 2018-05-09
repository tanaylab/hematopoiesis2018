# testing
You can call `testthat::test_package("scrdb")` before committing to make sure at least the basic functions work.
Also, Use the functions for regression tests in internal_tests/regression_tests.r: 
#### Before changing the code, type:
`dir = save_before_changing(basedir = "/net/mraid14/export/data/users/atanay/proj/scmat/regression_tests/", label = "before_changing_x", umis = umis)`
to save the output of basic clustering and visualization in directory *basedir/label_timestamp*. 
#### to test:
`.do_regression(dir, umis)`
This will only work if the changes in code do not affect the created objects. (for example, the change was only to improve runtime). It might be useful also if you know the knn object will change, but want to make sure other objects are not affected.