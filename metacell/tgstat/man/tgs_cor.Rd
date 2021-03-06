\name{tgs_cor}
\alias{tgs_cor}
\title{
Calculates correlation between the matrix columns
}
\description{
Calculates correlation between the matrix columns.
}
\usage{
tgs_cor(x, pairwise.complete.obs = F, spearman = F, tidy = F,
        threshold = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric matrix}
  \item{pairwise.complete.obs}{see below}
  \item{spearman}{if 'TRUE' Spearman correlation is computed, otherwise Pearson}
  \item{tidy}{if 'TRUE' data is outputed in tidy format}
  \item{threshold}{absolute threshold above which values are outputed in tidy format}
}
\details{
'tgs_cor' is very similar to 'package:stats::cor'. Unlike the latter it uses all available CPU cores to compute the correlation in a much faster way. The basic implementation of 'pairwise.complete.obs' is also more efficient giving overall great run-time advantage.

Unlike 'package:stats::cor' 'tgs_cor' computes the correlation only between matrix columns and implements only two modes of treating data containing NA, which are equivalent to 'use="everything"' and 'use="pairwise.complete.obs". Please refer the documentation of this function for more details.

'tgs_cor(x, spearman = F)' is equivalent to 'cor(x, method = "pearson")'
'tgs_cor(x, spearman = T)' is equivalent to 'cor(x, method = "spearman")'
'tgs_cor(x, pairwise.complete.obs = T, spearman = T)' is equivalent to 'cor(x, use = "pairwise.complete.obs", method = "spearman")'
'tgs_cor(x, pairwise.complete.obs = T, spearman = F)' is equivalent to 'cor(x, use = "pairwise.complete.obs", method = "pearson")'

'tgs_cor' can output its result in "tidy" format: a data frame with three columns named 'col1', 'col2' and 'cor'. Only the correlation values which abs are equal or above the 'threshold' are reported. Correlation of column number X with column Y is reported only if X < Y.
}
\value{
If 'tidy' is 'FALSE' - the matrix of correlation values, where val[X,Y] represents the correlation between columns X and Y of the input matrix. If 'tidy' is 'TRUE' - a data frame, where each row represents correlation between two pairs of columns.
}

\examples{
set.seed(seed = 0)
rows <- 100
cols <- 1000
vals <- sample(1 : (rows * cols / 2), rows * cols, replace = TRUE)
m <- matrix(vals, nrow = rows, ncol = cols)
m[sample(1 : (rows * cols), rows * cols / 1000)] <- NA

r1 <- tgs_cor(m, spearman = F)
r2 <- tgs_cor(m, pairwise.complete.obs = T, spearman = T)
}
\keyword{ ~correlation }

