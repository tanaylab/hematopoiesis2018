# run automatically using test_package("scrdb") or R CMD check.
# should be short.

library(testthat)
library(scrdb)

test_check("scrdb")
