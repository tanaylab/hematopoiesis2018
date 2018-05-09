#' Dimensionality reduction for gene expression
#'
#' Reduce a single cell matrix over a gene set, or through various dimensionality reduction strategies.
#'
#' @param Object .
#' @param type .
#' @param d_n .
#'
#' @export
#'
setGeneric("tgscm_reduce",
	function(.Object, type, d_n, ...) stnadrdGeneric("tgscm_reduce"))
setMethod("tgscm_reduce",
	signature = "tgScMat",
	definition =
	  function(.Object, type, d_n, ...) {
		if(type == "on_genes") {
			.tgscm_reduce_on_genes(.Object, d_n, ...)
		} else if(type == "pca") {
			.tgscm_reudce_pca(.Object, d_n, ...)
		}
	  }
)

.tgscm_reduce_on_genes = function(scmat, d_n, genes)
{
	if(length(setdiff(genes, rownames(scmat))) != 0) {
		ndiff = length(setdiff(genes, rownames(scmat)))
		message(ndiff, " genes out of ", length(genes), " missing in projection")
		return(NA)
	}
	tgScMat(scmat[genes,], stat_type=scmat@stat_type)
}

.tgscm_reduce_pca = function(scmat, d_n)
{
	#TBA
}
