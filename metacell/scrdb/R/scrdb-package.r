#' scrdb.
#'
#' @import dplyr
#' @import ggplot2
#' @import data.table
#' @import tgconfig
#' @import cowplot
#' @import tgstat
#' @import tgutil
#' @import igraph
#' @importFrom pdist pdist
#' @importFrom cluster silhouette
#' @name scrdb
#' @docType package

.onLoad <- function(libname, pkgname) {
	tgconfig::register_params(system.file('config/scrdb_params.yaml', package='scrdb'), package='scrdb')
	ggplot2::theme_set(ggplot2::theme_classic())
}
