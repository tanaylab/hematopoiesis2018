% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/temp_kmeans.R
\name{TGL.kmeans2}
\alias{TGL.kmeans2}
\title{This is a temporay TGL kmeans wrapper, untill Aviezer builds a kmeans/utils package}
\usage{
TGL.kmeans2(data, k, keep.log = FALSE, tmpdir = NULL, tmpfn = NULL,
  TGL.kmeans.bin = "/net/mraid14/export/data/tools/tgtools/TGLKMeans_static")
}
\arguments{
\item{data}{tgScMat}

\item{k}{integer}

\item{keep.log}{boolean}

\item{tmpdir}{string}

\item{tmpfn}{string}

\item{TGL.kmeans.bin}{string path to kmean binary file}
}
\value{
k-means clustering of the data
}
\description{
Effiecient implementation of kmeans++ algorithm
Allows keeping the output of the kmeans exec within the R return value
}
