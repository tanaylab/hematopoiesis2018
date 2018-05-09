#' Gene selection
#'
#' Strategies for gene selection given a single cell matrix
#'
#' @param scmat tgScMat data matrix
#' @param type .
#' @param ... .
#'
#' @import dplyr
#'
#' @export
#'
setGeneric("select_markers",
	function(scmat, gstat = NULL, type=c("any", "niche", "varmean", "sizecor"), ...)
	  standardGeneric("select_markers"))

setMethod("select_markers",
	signature = "tgScMat",
	definition =
	  function(scmat, gstat=NULL, type="any", n_marks = 10^4,
	           blacklist = NULL,expand = F,
	           fname = NULL, outdir = ".",
	           sd.thresh = 2.5,
	           mark.top3_T = NULL,
	           mark.min_tot = NULL,
	           mark.niche_T = NULL,
	           mark.min_var_mean = NULL,
	           mark.sz_cor_norm_max = NULL,
	           rseed = 1
	           ) {

	    oldseed = .set_seed(rseed)

	    if(is.null(gstat)) {
	      gstat = umi_gene_stat(scmat)
	    }

	    ### Set default parameters
	    if (sd.thresh>10 | sd.thresh < 1) {
	      sd.thresh = 2.5
	      warning("sd.thresh must be between 1 and 10. Changed to default value ", sd.thresh)
	    }

	    # filter low genes
	    mask = .filter_genes(gstat, blacklist,
	                         sd.thresh,
	                         mark.top3_T,
	                         mark.min_tot)

	    # select markers
  		if(type == "niche") {
  			marks = .tgscm_select_gene_niche(gstat, mask, sd.thresh, mark.niche_T)
  		} else if(type == "varmean") {
  		  marks = .tgscm_select_gene_varmean(gstat, mask, sd.thresh, mark.min_var_mean)
  		} else if(type == "sizecor") {
  		  marks = .tgscm_select_gene_sz_cor(gstat, mask, sd.thresh, mark.sz_cor_norm_max)
  		} else if(type == "any") {
  		  marks = .tgscm_select_gene_any(gstat, mask, sd.thresh,
  		                                 mark.niche_T,
  		                                 mark.min_var_mean,
  		                                 mark.sz_cor_norm_max)
  		} else {
  		  message("invalid gene selection method: ", type)
  		  return(NULL)
  		}
	    marks = .postprocess_select_gene(scmat, marks= marks, gstat = gstat,
	                                     mark.n_marks = n_marks,
	                                     blacklist = blacklist, expand_markers = expand,
	                                     fname = fname, outdir = outdir)
	    .restore_seed(oldseed)

	    return(marks)
	  }
)

.filter_genes = function(gstat, blacklist,
                         sd.thresh,
                         mark.top3_T,
                         mark.min_tot
)
{
  gstat_over0 = gstat[gstat$var>0, ]

  if(is.null(mark.top3_T)) {
    mark.top3_T = max(4, quantile(gstat_over0$ds_top3, 0.8))
  }
  if(is.null(mark.min_tot)) {
    #TODO: tune this!
  	mark.min_tot = max(50, min(100, round(median(gstat_over0$tot))))
  }

  mask = gstat$ds_top3 >= mark.top3_T & gstat$tot >= mark.min_tot &
    !(gstat$name %in% blacklist)

  cat("Selecting marker genes: \nFiltered",sum(!mask), "genes with total expression under" , mark.min_tot,
      ", or 3rd highest count after down-sampling under", mark.top3_T ,"\n")
  return(mask)
}

# calls all methods of gene selection, returns the union of the markers,
# ordered by the sum of their rankings.
.tgscm_select_gene_any = function(gstat,mask,
                                  sd.thresh,
                                  mark.niche_T,
                                  mark.min_var_mean,
                                  mark.sz_cor_norm_max) {
  m = list()
  m[[1]] = .tgscm_select_gene_niche(gstat, mask, sd.thresh, mark.niche_T)
  m[[2]] = .tgscm_select_gene_varmean(gstat, mask, sd.thresh, mark.min_var_mean)
  m[[3]] = .tgscm_select_gene_sz_cor(gstat, mask, sd.thresh, mark.sz_cor_norm_max)
  marks = unique(unlist(m))
  ranks = matrix(0,length(m), length(marks), dimnames = list(NULL, marks))
  for (i in 1:length(m)) {
    genes = m[[i]]
    ranks[i,genes] = length(genes):1
  }
  return(marks[order(colSums(ranks), decreasing = T)])
}


.tgscm_select_gene_niche = function(gstat,mask,
			sd.thresh, # used for setting mark.niche_T
			mark.niche_T = NULL)
{
	# Set default parameters
	if(is.null(mark.niche_T)) {
		mark.niche_T = .calc_thresh(gstat$niche_norm, sd.thresh)
	}

  # select
	marks = gstat[mask,] %>%
    filter(niche_norm >= mark.niche_T) %>%
    arrange(desc(niche_norm)) %>%
    select(name)

	cat("found ", nrow(marks), " markers with niche score over", round(mark.niche_T,2), "\n")
	return(as.character(marks[,1]))
}


# returns marker with low correlation to cell size, ordered by score
.tgscm_select_gene_sz_cor = function(gstat, mask,
                                     sd.thresh,
                                     mark.sz_cor_norm_max = NULL) # -0.1 ,
{
  if(is.null(mark.sz_cor_norm_max)) {
  	mark.sz_cor_norm_max = .calc_thresh(gstat$sz_cor_norm, sd.thresh, low = T)
  }
  ### select and reorder markers
  marks = gstat[mask,] %>%
    filter(sz_cor_norm <= mark.sz_cor_norm_max) %>%
    arrange(sz_cor_norm) %>%
    select(name)
  cat("found ", nrow(marks), " markers with normalized size correlation under",
      round(mark.sz_cor_norm_max,2),"\n")
  return(as.character(marks[,1]))
}


.tgscm_select_gene_varmean = function(gstat, mask, sd.thresh, mark.min_var_mean)
{

  if(is.null(mark.min_var_mean)) {
    mark.min_var_mean = .calc_thresh(gstat$ds_vm_norm, sd.thresh)
  }
  ### select and reorder markers
  marks = gstat[mask,] %>%
    filter(ds_vm_norm >= mark.min_var_mean) %>%
    arrange(desc(ds_vm_norm)) %>%
    select(name)
  cat("found ", nrow(marks), " markers with normalized var/mean over", mark.min_var_mean,"\n")
  return(as.character(marks[,1]))
}

# Calculate a threshold for a score based on standard deviations from median,
# after removing zero scores.
.calc_thresh = function(score_norm, sd.thresh, lower.tail = F) {

	score_norm = score_norm[score_norm != 0]
	q95 = quantile(score_norm, c(0.95))
	q05 = quantile(score_norm, c(0.05))
	score_norm = score_norm[score_norm >= q05 & score_norm <= q95]
	sign = ifelse(lower.tail, -1, 1)
	return(median(score_norm) + sign * sd(score_norm) * sd.thresh)
}

#' Calculate basic statistics on a matrix
#'
#' @param mat The input matrix
#' @param niche_quantile A value between 0 and 1.
#'
#' @return Returns a dataset that contains statistic for all genes that passed the filtering stage.
#' Columns starting with ds contain UMI statistics after downsampling,
#' columns starting with n contain UMI statistics after normalizing UMIs so that the number of UMIs
#' per cell sums to 1.
#'
#' The columns are:
#'
#' \describe{
#'   \item{tot}{Total gene expression}
#'   \item{var}{Gene variance}
#'   \item{is_on_count}{Number of cells in which the gene is expressed}
#'   \item{sz_cor}{Correlation with cell size}
#'   \item{sz_cor_norm}{sz_cor after subtracting the trend}
#'   \item{niche_stat}{How many of the genen's umis are found in X\% of the most highly expressing cells. (regularized) }
#'   \item{niche_norm}{niche_stat after subtracting the niche_norm trend: median niche_norm value of genes with similar total expression}
#'   \item{n_mean}{Mean after normalization}
#'   \item{ds_top1}{Largest count, after downsampling}
#'   \item{ds_top2}{2nd largest count, after downsampling}
#'   \item{ds_top3}{3rd largest count, after downsampling}
#'   \item{ds_mean}{Mean on downsampled data}
#'   \item{ds_var}{Variance on downsampled data}
#'   \item{ds_log_varmean}{log2 of ds_var/ds_mean }
#'   \item{ds_vm_norm}{ds_log_varmean after subtracting the trend}
#'   \item{ds_is_on_count}{Number of cells in which the gene is expressed, after down sampling}
#'   \item{downsample_n}{Number of UMIs used for downsampling}
#'
#' }
# #'   \item{max_pk
# #' max_ratio
# #' var_meanpk
# #' ds_var_meanpk
#'
#'
#' @import zoo
#'
#' @export
#'
umi_gene_stat = function(scmat,
                         niche_quantile = 0.2, #TODO: change to 0.1?
                         rseed=321, downsample_n = NULL)
{
	oldseed = .set_seed(rseed)

  mat = as.matrix(scmat@mat)
  cat("Calculating gene statistics... ")
  K = 1000 # for normalizing by library size

  if (niche_quantile >1 | niche_quantile < 0 ) {
    stop("niche_quantile must be between 0 and 1, got ", niche_quantile)
  }

  # filtering cells with too many or too few umis
	f_oversize = colSums(mat) > quantile(colSums(mat),0.95)*2 |
							colSums(mat)<100
	mat = mat[,!f_oversize]

	# returns how many of the gene's umis are found in
	# X% of the most highly expressing cells. (regularized)
	quant_mean = function(x, k_reg) {
		n = length(x);
		up=sum(tail(sort(x), n=round(n*niche_quantile))); # sum umis in the top 20%
		return((k_reg+up)/(k_reg+sum(x)));
	}
	N = sum(!f_oversize)

	f_g = rowSums(mat) > 10
	mat_n = t(t(mat[f_g,])*(K/colSums(mat)))

  if(is.null(downsample_n)) {
    downsample_n = min(round(quantile(colSums(mat), 0.5)),
                  max(750, round(quantile(colSums(mat), 0.05))))
  }

	mat_ds = .downsamp(mat, downsample_n)
	mat_ds = mat_ds[f_g,]

  n_ds = ncol(mat_ds)
	n = ncol(mat)

	gene_stat = data.frame(stringsAsFactors = F,
	  name = rownames(mat[f_g,]),
		tot = rowSums(mat[f_g,]),
		var = round(matrixStats::rowVars(mat[f_g,]),3),
		# var = round(apply(mat[f_g,], 1, var),3),

		niche_stat =
		  apply(mat[f_g,], 1, quant_mean, k_reg = 10),

		n_mean = round(matrixStats::rowMeans2(mat_n),3),

		# ds_top1 = apply(mat_ds, 1, max),
		# ds_top2 = apply(mat_ds, 1,
		# 								function(x) sort(x, partial=n_ds-1)[n_ds-1]),
		# ds_top3 = apply(mat_ds, 1,
		# 											 function(x) sort(x, partial=n_ds-2)[n_ds-2]),

		ds_top1 = matrixStats::rowMaxs(mat_ds) ,
		ds_top2 = matrixStats::rowOrderStats(mat_ds, which = n_ds-1) ,
		ds_top3 = matrixStats::rowOrderStats(mat_ds, which = n_ds-2) ,

		# mean_pk = round(1000*apply(mat_n, 1, mean),3),
		# max_pk = round(1000*apply(mat_n, 1, max),3),
		# is_on_count = apply(mat[f_g,], 1, function(x) sum(as.numeric(x) >= 1)),
		is_on_count = matrixStats::rowCounts(mat[f_g,] > 0),

		ds_is_on_count = matrixStats::rowCounts(mat_ds > 0),
		ds_var = round(matrixStats::rowVars(mat_ds),3),
		ds_mean = round(matrixStats::rowMeans2(mat_ds),3)
	)

	ctot = colSums(mat)
	gene_stat$sz_cor = round(apply(mat[f_g,], 1, function(x) { cor(x, ctot) }),3)

	tot_ord = order(gene_stat$tot)
	cor_sz_ord = gene_stat$sz_cor[tot_ord]
	cmin = median(cor_sz_ord[1:101])
	cmax = median(cor_sz_ord[(length(cor_sz_ord)-101):length(cor_sz_ord)])
	sz_cor_trend = zoo::rollmedian(cor_sz_ord, 101, fill=c(cmin, NA, cmax))

	gene_stat$sz_cor_norm[tot_ord] =
		gene_stat$sz_cor[tot_ord] - sz_cor_trend

	# niche_stat =  how many of the genen's umis are found in
	# 20% of the most highly expressing cells. (regularized)
	niche_sz_ord = gene_stat$niche_stat[tot_ord]
	cmin = median(niche_sz_ord[1:101])
	cmax = median(niche_sz_ord[(length(niche_sz_ord)-101):length(niche_sz_ord)])
	niche_trend = zoo::rollmedian(niche_sz_ord, 101, fill=c(cmin, NA, cmax))

	gene_stat$niche_norm[tot_ord] =
		gene_stat$niche_stat[tot_ord] - niche_trend

	# var/mean
	ds_ord = order(gene_stat$ds_mean)

	gene_stat$ds_log_varmean = ifelse(gene_stat$ds_mean>0.01, log2((0.001+gene_stat$ds_var)/(gene_stat$ds_mean+0.001)), 0)
	vm_sz_ord = gene_stat$ds_log_varmean[ds_ord]
	cmin = median(vm_sz_ord[1:101])
	cmax = median(vm_sz_ord[(length(vm_sz_ord)-101):length(vm_sz_ord)])
	vm_trend = zoo::rollmedian(vm_sz_ord, 101, fill=c(cmin, NA, cmax))

	gene_stat$ds_vm_norm[ds_ord] =
		gene_stat$ds_log_varmean[ds_ord] - vm_trend

	gene_stat$downsample_n = downsample_n

	# currently not used
	# gene_stat$max_ratio = gene_stat$max_pk/gene_stat$mean_pk
	# gene_stat$var_meanpk = gene_stat$var/gene_stat$mean_pk
	# gene_stat$ds_var_meanpk = gene_stat$ds_var/gene_stat$mean_pk

	rownames(gene_stat) = gene_stat$name


	gene_stat = gene_stat[,c("name","tot","var","is_on_count", "sz_cor","sz_cor_norm",
	                      "niche_stat",  "niche_norm", "n_mean",
	                      "ds_top1", "ds_top2", "ds_top3",
	                      "ds_mean", "ds_var", "ds_log_varmean","ds_vm_norm", "ds_is_on_count", "downsample_n")]

	.restore_seed(oldseed)
	cat("..done\n")

	return(gene_stat)
}

# should incluse: black list genes, plot markers
# option: expand list to correlated genes
.postprocess_select_gene = function(scmat, marks, gstat = gstat, mark.n_marks,
                                    blacklist = NULL, expand_markers,
                                    fname = NULL, outdir = NULL,
                                    min_markers_per_cell = 2) 
{

  # expansion to correlated genes (here)
  if(expand_markers) {
    marks = setdiff(marks, blacklist)
    marks = expand_markers(scmat, marks, gstat)
  }

  # removing black list
  marks = setdiff(marks, blacklist)
  marks = head(marks,n=mark.n_marks)

  fmat=as.matrix(scmat@mat[marks,])

  # check
  count = colSums(fmat>0)
  if(min(count) < min_markers_per_cell) {
    warning("Chosen marker genes do not cover all cells. Try using other parameters; for example - lowering sd.thresh")
  }
  if(!is.null(fname) & !is.null(outdir)) {
    plot_markers(mat = as.matrix(scmat@mat), marks, gstat = gstat, fname=fname, outdir=outdir)
  }
  return(marks)
}

#' Expand markers list with correlated genes
#'
#' Expand markers list with correlated lowly expressed genes.
#'
#' @return the original list of marker, with the new markers added at the end
#'
#' @param scmat tgScMat
#' @param marks Original list of markers
#' @param cor_T Threshold for minimal correlation for addition
#'
######' @export
expand_markers = function(scmat, marks, cor_T = 0.1, sd.thresh = 2.5, mark.min_tot = 20) {

  warning("expand_markers not supported yet")
  # TODO: rewrite this whole function so it (1) calls gal's code,
  # (2) comutes only relevant corr: corr=cor(t(mat[marks,]),t(mat[cand,]), m="spearman")
  # Uses more lowly expressed genes

  new_marks = c()
  # mask = gstat$tot <= max_tot & gstat$is_on_count >= min_cells_expressing
  # sink("/dev/null")
  cand = select_markers(scmat, sd.thresh = sd.thresh,
  												 mark.min_tot = mark.min_tot,
                           expand = F) # otherwise we get eternal loop
  # sink()

  cand = setdiff(cand, marks)

  # use gal's code
  mat = as.matrix(scmat@mat)
  mat = t(t(mat)/colSums(mat))
  corr=cor(t(mat[unique(union(cand, marks)),]), m="spearman")

  subcorr = corr[marks, !rownames(corr) %in% marks]
  friends = apply(subcorr, 1, which.max)
  mask2 = laply(1:length(friends), function(i) {subcorr[names(friends)[i], friends[i]] >= min_cor})
  if (sum(mask2)>0) {
    new_marks = colnames(subcorr)[friends[mask2]]
    cat("Adding",  sum(mask2), " markers which are correlated with the", length(marks) , "markers found so far\n")
  }



  return (c(marks, new_marks))
}

#' cluster markers
#' @export
tgscm_cluster_markers = function(mat, marks, K)
{
  mat_n = t(t(mat)/colSums(mat))
  score = gstat$niche_norm[marks] - gstat$sz_cor_norm[marks]
  marks_to_clust = tail(marks[order(score)],n=800)

  markcor=cor(t(as.matrix(mat_n[marks_to_clust,])), m="spearman")
  hc = hclust(dist(markcor),"ward.D2")

  return(cutree(hc, K))
}

#' plot statistics about the chosen markers
#' @export
plot_markers = function(mat, marks, gstat, fname, outdir,
                            shades=colorRampPalette(c("navy", "white","brown"))(100),
                            dendrogram.lab.cex = 2,
                            h=4500, w=1800
                             # shades=colorRampPalette(c("chartreuse4", "darkseagreen4", "white", "gray", "blue", "red", "yellow"))(1000)
) {
  if(fname != "interactive") {
    fpath = sprintf("%s/%s", outdir, fname)
    png(fpath, h=h, w=w);
    oldpar = par(cex.lab=3, cex.main=3, cex.axis = 2)
  }


  layout(matrix(c(1,1,1,1,1,2,3,4,5,6), nrow=5),heights=c(1,1,1,1,1))
  oldmar = par(mar = c(6,2,2,6))


  #
  mat_n = t(t(mat)/colSums(mat))
  score = gstat$niche_norm[marks] - gstat$sz_cor_norm[marks]
  marks_to_clust = tail(marks[order(score)],n=800)

  # TODO: this is bad, lowly expressed genes are correlated for no reason.
  #       integrate Gal's code.
  markcor=cor(t(as.matrix(mat_n[marks_to_clust,])), m="spearman")
  hc = hclust(dist(markcor),"ward.D2")

  # 1
  plot(as.dendrogram(hc), horiz=T, nodePar = list(pch=NA,lab.cex = dendrogram.lab.cex))
  if(length(marks) == length(marks_to_clust)) {
    marks = marks[hc$order]
  }

  # 2
  par(mar = c(6,6,6,2))
  image(markcor[hc$order, hc$order], col=shades,zlim=c(-1,1), main = "marker correlations")

  # 3
  plot(gstat$niche_norm, gstat$sz_cor_norm, cex=0.8, pch=19,
       xlab = "normalized niche score", ylab="normalized size correlation");
  points(gstat[marks,"niche_norm"], gstat[marks, "sz_cor_norm"], cex=1.2, pch=19, col="red");
  grid()

  # 4
  plot(log2(gstat$tot), gstat$niche_stat, cex=0.8, pch=19,
       xlab ="log2 total expression", ylab="niche score")
  points(log2(gstat[marks,"tot"]), gstat[marks,"niche_stat"], cex=1.2, col="red", pch=19)
  grid()

  # 5
  plot(log2(gstat$tot), gstat$sz_cor, cex=0.8, pch=19,
       xlab ="log2 total expression", ylab="size correlation")
  points(log2(gstat[marks,"tot"]), gstat[marks,"sz_cor"], cex=1.2, col="red", pch=19)
  grid()

  #6
  f=gstat$ds_log_varmean != 0
  plot(log2(gstat$ds_mean[f]), gstat$ds_log_varmean[f], cex=0.8, pch=19,
       xlab = "log2(mean down-sampled)", ylab="log2(var/mean down-sampled)")
  points(log2(gstat[marks,"ds_mean"]), gstat[marks,"ds_log_varmean"],
         cex=1.2, col="red", pch=19)
  grid()

  if(fname != "interactive") {
    par(oldpar)
    dev.off()
  }
  par(oldmar)
}

