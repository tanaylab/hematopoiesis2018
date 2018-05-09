
# internal utility functions

# rowFunction is a function that works on rows, like rowMeans
# # much faster than tapply
.row_stats_by_factor = function (data, fact, rowFunction = rowMeans) {
	u = as.character(sort(unique(fact)))
	fact[is.na(fact)] = F
	n=length(u)
	centers = matrix(NA,dim(data)[1], n, dimnames = list(rownames(data), u))
	for (i in u) {
		if(sum(fact==i, na.rm=T)>1) {
			centers[,i] = rowFunction(data[,fact==i,drop=F])
		} else {
			centers[,i] = data[,fact==i]
		}
	} # much faster than tapply
	return(centers)
}



.downsamp = function (umis, n, replace = F) {
	# umis = umis[, colSums(umis) >= n]
	# ret = matrix(0, nrow = nrow(umis), ncol = ncol(umis))
	# for(i in 1:ncol(ret)) {
	# 	a = tabulate(sample(rep(1:nrow(umis),times=umis[,i]),replace=replace,size=n))
	# 	ret[1:length(a),i] = a
	# }

	m = nrow(umis)
	.downsamp_one=function(v,n, replace = F){
		a = tabulate(sample(rep(1:length(v),times=v),replace=replace,size=n),nbins=m)
		return (a)
	}
	ret = apply(umis[, colSums(umis) >= n], 2, .downsamp_one, n)
	rownames(ret) = rownames(umis)
	return(ret)
}

# .downsamp_old = function (umis, n, replace = F) {
#
# 	downsamp_one=function(v,n, replace = F){
# 		hist(sample(rep(1:length(v),times=v),replace=replace,size=n),0.5+0:length(v),plot=F)$counts
# 	}
# 	return(apply(umis[, colSums(umis) >= n], 2, .downsamp_one, n))
# }


# from tlsrc/analysis/common/fread.r; one day this will be in a "common" package
# returns rownames as the first column, named row.var
### @importFrom data.table fread
fread <- function(...) data.table::fread(..., data.table=FALSE)

fread_rownames <- function(..., row.var='rowname', set_rownames = F) {
	params <- list(...)
	header <- strsplit(readLines(params[[1]], n=1, warn=FALSE), '\t', fixed=TRUE)[[1]]

	params$header = F; params$skip = 1; params$col.names = c(row.var, header)

	mat = do.call(fread,params)
	if (set_rownames) {
		rownames(mat) = mat[,1]
		mat = mat[,-1]
	}
	return(mat)
}

# save random number generator status, and set seed
.set_seed = function(rseed = 1) {
	if (exists(".Random.seed", .GlobalEnv)) {
		oldseed = .GlobalEnv$.Random.seed
	}
	else {
		oldseed = NULL
	}
	set.seed(rseed)

	return(oldseed)
}

# restor random number generator status
.restore_seed = function(oldseed) {
	if (!is.null(oldseed)) {
		.GlobalEnv$.Random.seed = oldseed
	}
	else {
		rm(".Random.seed", envir = .GlobalEnv)
	}
}

# wrap for opening a plot (png or ps)
.plot_start = function(fn, w, h, device=get_param("plot_device"), res=get_param("plot_ppi")) {
	if (device == "png") {
		png(filename=sub("ps$", "png", fn), width=w, height=h, res=res)
	}
	else if (device == "ps") {
		postscript(file=sub("png$", "ps", fn), width=w/res, height=h/res)
	}
	else {
		stop(sprintf("unknown output device type: %s", device))
	}
}

# translate values vector to integers (for choosing colors by value from a color pallette
.vals_to_n = function(v, n = 1000, zlim = NULL) {

		if (is.null(zlim)) {
		zlim = range(v)
	}
	v = pmin(pmax(v, zlim[1]), zlim[2])
	v = v - min(v)
	1 + floor(v * ((n - 1) / max(v)))
}

# matrix correlation wrap. Parameter selects which function to use, cor or tgs_cor
.cor = function(x, pairwise.complete.obs = F, spearman = F, tidy = F,	threshold = 0)
{
	if (get_param("use_tgs_cor")) {
		return (tgs_cor(x, pairwise.complete.obs, spearman, tidy, threshold))
	}
	else {
		return(cor(x, method=ifelse(spearman, "spearman", "pearson"), use=ifelse(pairwise.complete.obs, "pair", "all")))
	}
}
