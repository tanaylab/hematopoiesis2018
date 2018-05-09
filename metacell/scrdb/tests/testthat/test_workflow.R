# called using test_package("scrdb"), or R CMD check

source("./workflow.R")
data(toydata)
call_workflow_tests(toydata)
