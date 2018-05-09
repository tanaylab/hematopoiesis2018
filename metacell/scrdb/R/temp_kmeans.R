

############################################################################################
#' This is a temporay TGL kmeans wrapper, untill Aviezer builds a kmeans/utils package
#'
#' Effiecient implementation of kmeans++ algorithm
#' Allows keeping the output of the kmeans exec within the R return value
#'
#' @param data tgScMat
#' @param k integer
#' @param keep.log boolean
#' @param tmpdir string
#' @param tmpfn string
#' @param TGL.kmeans.bin string path to kmean binary file
#'
#' @return k-means clustering of the data
#' @export
#'
## ' @examples
TGL.kmeans2 <- function(data, k,
                        keep.log=FALSE,
                        tmpdir=NULL, tmpfn=NULL,
                        TGL.kmeans.bin = "/net/mraid14/export/data/tools/tgtools/TGLKMeans_static")
{
  tmpfile <- ''
  if (is.null(tmpdir)) {
    tmpdir <- tempdir()
  }
  if (is.null(tmpfn)) {
    tmpfile <- tempfile(pattern='tgl_kmeans.', tmpdir=tmpdir)
  } else {
    tmpfile <- file.path(tmpdir, tmpfn)
  }

  write.table(data, tmpfile, sep='\t', quote=FALSE, col.names=NA, row.names=TRUE)

  km_log <- system2(TGL.kmeans.bin, c(tmpfile, k, 'euclid', '-allow_nas=1'), stdout=keep.log, stderr=keep.log)
  kclust <- data.table::fread(paste0(tmpfile, '.kclust'), sep='\t', header=TRUE, data.table=FALSE)
  center <- data.table::fread(paste0(tmpfile, '.center'), sep='\t', header=FALSE, data.table=FALSE)

  kclust$clust <- kclust$clust + 1

  km <- list()

  km$cluster <- kclust$clust
  names(km$cluster) <- kclust$id

  km$centers <- as.matrix(center[,-1])
  colnames(km$centers) <- colnames(data)

    if (keep.log) {
    km$log <- km_log
  }

  # # reorder clusters
  if (min(apply(km$centers, 1, var)) == 0) {
    warning("standard deviation of kmeans center is 0")
    # a = km$centers; a[,1] = a[,1] + 5
    # centers_hc = hclust(dist(cor(t(a))), "ward.D2")
  } else {
    centers_hc = hclust(dist(cor(t(km$centers))), "ward.D2")
    tomap = rep(NA, k)
    tomap[centers_hc$order] = 1:k
    km$cluster = tomap[km$cluster]
    names(km$cluster) = kclust$id
    km$centers = km$centers[centers_hc$order,]
  }
  rownames(km$centers) = 1:k

  km$size <- tapply(km$clust, km$clust, length)
  
  return(km)
}
