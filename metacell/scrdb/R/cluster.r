#' A clustered matrix

#' Representing a clustered matrix, and implementing several basic algorithms for constructing it
#' from a matrix
#'
#' @slot scmat tgScMat. 
#' @slot feat_mat Matrix. 
#' @slot clusts vector. 
#' @slot orphan_clust_id numeric. 
#' @slot clust_fp matrix. 
#' @slot clust_gcov matrix. 
#' @slot alg_params list. 
#' @slot nclust integer. 
#' @slot knn_ordered matrix. k_knn X ncells matrix,  column i contains cell i's k nearest neighbors, ordered fronm closest to farthest.
#' @slot knn_is_balanced logical. 
#' @slot m_knn dgCMatrix. 
#' @slot rseed numeric. 
#'
#' @importClassesFrom Matrix dgCMatrix
#' @exportClass tgScMatClust
#' 
#### @export tgScMatClust
#### @exportClass tgScMatClust
#
#### @include scmat.r

tgScMatClust <- setClass(
  "tgScMatClust",
  slots = c(
    scmat = "tgScMat",
    feat_mat = "Matrix",
    clusts = "vector",
    orphan_clust_id = "numeric",
    clust_fp = "matrix",
    clust_gcov = "matrix",
    alg_params = "list",
    nclust = "integer",
    knn_ordered = "matrix",
    knn_is_balanced = "logical",
    m_knn = "dgCMatrix",
    rseed = "numeric")
)

#' Initialize tgScMatClust object
#'
#' @param scmat 
#' @param features 
#' @param norm_features. Logical - whether to normalize the features
#' @param alg_type 
#' @param rseed 
#' @param factor_by_median Multiple 
#' @param norm_size 
#' @param k_scale_umi 
#'
#' @return
#' @export
#'
#' @examples
.scc_init = function(scmat,
										 features = NULL,
										 norm_features = T,
										 alg_type,
										 rseed = 1,
										 factor_by_median = get_param("scc_factor_by_median"),
										 norm_size = get_param("scc_norm_size"),
										 k_scale_umi = get_param("scc_k_scale_umi"))
{
	scc = new("tgScMatClust")
	scc@scmat = scmat

	scc@feat_mat = create_feature_mat(scmat, features, norm_features = norm_features,
				factor_by_median = factor_by_median,
				norm_size = norm_size, k_scale_umi = k_scale_umi)

	scc@alg_params = list(type=alg_type, rseed = rseed) # this list is extended after clustering

	return(scc)
}


#' Remove bad cells to an 'orphan' cluster and update clusters slots
#' 
#' @param scl 
#'
#' @param bad_cells 
#'
#' @export
scc_move_cells_to_orphan = function(scl, bad_cells)
{
	good_cells = setdiff(colnames(as.matrix(scl@scmat@mat)), bad_cells)
	scl@clusts[good_cells] = as.integer(as.factor(scl@clusts[good_cells]))
	scl@orphan_clust_id = max(scl@clusts[good_cells])+1
	scl@clusts[bad_cells] = scl@orphan_clust_id
	scl = .scc_postprocess(scl)
	scl@nclust = as.integer(scl@orphan_clust_id)
	return(scl)
}

#' @param scc 
#'
#' @return
.scc_postprocess = function(scc) {
	message("will comp fp")
	scc@clust_fp = .calc_clusts_fp(scc)

	message("will comp gcov")
	scc@clust_gcov = .row_stats_by_factor(scc@scmat@mat > 0, fact = scc@clusts, rowFunction = rowMeans)
	
	scc@nclust = length(unique(scc@clusts))

	return(scc)
}

#' Cluster single cell data
#'
#' Use one of the avialable clustering methods
#'
#' @param scmat A tgScMat expression matrix.
#' @param features Either  a matrix of cell features (gene expression, PCA principle componant, etc.),
#' or a list of markers to be extracted from scmat. Clustering will be performed on these features.
#' If left blank, markers will be selected using the default method.
#'
#' @param alg_type One of: knn, kmeans.
#'
#' @return A tgScMatClust object
#'
#' @export
scc_cluster = function(scmat, features = NULL, alg_type = "knn", ...) {

	if(alg_type == "kmeans") {
		.Object = scc_cluster_kmeans(scmat, features, ...)
	} else if(alg_type == "knn") {
		.Object = scc_cluster_knn_graph(scmat, features, ...)
	} else if(alg_type == "hclust") {
		.Object = scc_cluster_hclust(scmat, features, ...)
	# } else if(alg_type == "multinom_em") {
	# 	.Obect = scc_cluster_multinom_em(scmat, features, ...)
	} else {
		stop("algorithm ", alg_type, " is not recognized")
	}
}

#' Cluster single cell data using kmeans
#'
#' @param scmat A tgScMat expression matrix.
#' @param features Either  a matrix of cell features (gene expression, PCA principle componant, etc.),
#' or a list of markers to be extracted from scmat. Clustering will be performed on these features.
#' If left blank, markers will be selected using the default method.
#'
#' @param K number of clusters
#'
#' @param norm_features Apply normalization on features
#' @param k_scale_umi 
#' @param factor_by_median 
#' @param norm_size 
#' @param rseed 
#'
#' @return A tgScMatClust object
#'
#' @export
scc_cluster_kmeans = function(scmat,
															features = NULL,
															K = NA,
															norm_features = T,
															k_scale_umi = get_param("scc_k_scale_umi"),
															factor_by_median = get_param("scc_factor_by_median"),
															norm_size = get_param("scc_norm_size"),
															rseed = 1)
{
	oldseed = .set_seed(rseed)
	scm_cl = .scc_init(
		scmat,
		features,
		norm_features = norm_features,
		alg_type = "kmeans",
		factor_by_median = factor_by_median,
		norm_size = norm_size,
		k_scale_umi = k_scale_umi
	)
	
	if(is.na(K)) {
		K = max(2,ceiling(ncol(scm_cl@feat_mat)/60))
	}
	cat("Running kmeans with K =", K, "...")

	km = TGL.kmeans2(as.matrix(t(scm_cl@feat_mat)), K)
	scm_cl@clusts = km$cluster
	scm_cl@alg_params =  c(scm_cl@alg_params, list (K=K, k_scale_umi = k_scale_umi,
																									factor_by_median = factor_by_median, norm_size = norm_size))

	scm_cl = .scc_postprocess(scm_cl)
	.restore_seed(oldseed)

	return(scm_cl)
}

#' Cluster single cell data using KNN graph
#' 
#' @param scmat A tgScMat expression matrix.
#' @param features Either  a matrix of cell features (gene expression, PCA principle componant, etc.),
#' or a list of markers to be extracted from scmat. Clustering will be performed on these features.
#' If left blank, markers will be selected using the default method.
#' @param norm_features 
#' @param k_scale_umi 
#' @param min_clust_size 
#' @param use_balanced_knn 
#' @param k_knn 
#' @param do_consolidation 
#' @param rseed 
#'
#' @return A tgScMatClust object
#'
#' @export
scc_cluster_knn_graph = function(scmat,
																 features = NULL,
																 norm_features = T,
																 k_scale_umi = get_param("scc_k_scale_umi"),
																 min_clust_size = get_param("scc_min_clust_size"),
																 use_balanced_knn = get_param("scc_use_balanced_knn"),
																 k_knn = get_param("scc_k_knn"),
																 ratio=get_param("scc_clust_size_to_k_knn_ratio"), 
																 do_consolidation = get_param("scc_knn_consolidate"),
																 rseed = 1)
{
	oldseed = .set_seed(rseed)
	scm_cl = .scc_init(scmat, features, norm_features = norm_features, alg_type = "knn", factor_by_median = F,
										 norm_size = F, k_scale_umi = k_scale_umi)

	
	if(is.null(k_knn)) {
		k_knn = floor(max(75, min(500, scm_cl@scmat@ncells / 40)))
	}
	if(is.null(min_clust_size)) {
		min_clust_size = min(floor(ratio * k_knn), 10)
	}

	if (min_clust_size > ratio * k_knn) {
		min_clust_size = floor(ratio * k_knn)
		cat("reducing minimal cluster size to ", min_clust_size)
	}
	cat("Running knn clust with K =", k_knn,"...\n")

	#must filter low complexity cells? of detect by knn?

	cells = scm_cl@scmat@cells
	ncells = length(cells)

	if(ncells < k_knn*2) {
		stop("K ", k_knn, " for clustering is too large given only ", ncells, " cells")
	}
	
	if(use_balanced_knn) {
		scm_cl = .scc_comp_balanced_knn_matrix(scm_cl, k_knn)
	} else {
		scm_cl = .scc_comp_knn_matrix(scm_cl, k_knn)
	}

	clusts = tg_piecewise_graph_cover(scm_cl@m_knn, 
			min_clust_size=min_clust_size, 
			consolidate=do_consolidation)

	clusts1 = as.integer(as.factor(clusts))
	names(clusts1) = names(clusts)
	clusts= clusts1
	message("will reorder")
	clusts = .reorder_knn_clusts(scm_cl@feat_mat, clusts)

	scm_cl@clusts = clusts
	scm_cl@alg_params = c(scm_cl@alg_params, list(k_knn = k_knn, k_scale_umi = k_scale_umi,
				min_clust_size = min_clust_size, do_consolidation = do_consolidation))

	message("will postproc clusts")
	scm_cl = .scc_postprocess(scm_cl)
	.restore_seed(oldseed)

	return(scm_cl)
}

#' bootstrap clustering: sample cells, cluster them (using the supplied scl object clustering parameters), build co-clust matrix and hclust that matrix 
#'
#' @param scl. Initial knn clustering result. 
#' @param scdb. scdb instance, will save co-clust matrix there if supplied.
#' @param k_cut. Number of clusters to extract from the co-clust hclust tree 
#' @param N_boot. Number of boostrap iterations
#' @param boot_ratio. Fraction of cells to sample in each iteration.
#' @param n_procs. Number of processes to use 
#' @param min_clust_size. Pool cells in clusters smaller than min_clust_size to a single outlier cluster.
#' @param verbose. print progress messages.
#' @param rseed 
#'
#' @return
#' @export
#'
#' @examples
scc_bootstrap_clusts = function(scl,
																scdb = NULL,
																k_cut = get_param("scc_bootstrap_k_cut"),
																N_boot = get_param("scc_bootstrap_n_boot"),
																boot_ratio = get_param("scc_bootstrap_ratio"),
																n_procs = get_param("scc_bootstrap_n_procs"),
																min_clust_size = get_param("scc_bootstrap_min_clust_size"),
																verbose = get_param("scc_bootstrap_verbose"),
																rseed=1)
{

	oldseed = .set_seed(rseed)
	if(length(scl@knn_ordered) <= 1) {
		stop("Run bootstrap after initial knn clustering")
	}

	doMC::registerDoMC(n_procs)
	N = ncol(scl@knn_ordered)
	message("boot with ", N, " nodes")
	if(N > 10000) {
		distrib_coclust = F
	} else {
		distrib_coclust = T
	}

	N_boot_each = round(N_boot/n_procs)
	boot_coclust = function(x) {
		tg_piecewise_knn_graph_cover_bootstrap(
			mat=NULL, 
			k_knn = scl@alg_params[["k_knn"]], 
			min_clust_size = scl@alg_params[["min_clust_size"]], 
			boot_ratio = boot_ratio, 
			N_boot = N_boot_each,
			compute_coclust = distrib_coclust,
			full_knn_ord = scl@knn_ordered,
			k_expand=4,
			verbose = verbose)
	}

	res <- plyr::alply(1:n_procs, 1, function(x) boot_coclust(), .parallel=TRUE)
	if(!distrib_coclust) {
		cc = matrix(0,nrow=N, ncol=N)
		cosamp = matrix(0,nrow=N, ncol=N)
		for(co in res[-1]) {
			for(i in 1:length(co$clusts)) {
			   clst = co$clusts[[i]]
			   boot_cells = co$nodes[[i]]
			   isclust_ci = diag(max(clst))[,clst]
			   coclust_ij = t(isclust_ci) %*% isclust_ci
			   if (verbose) { message("will add tot coclus ", paste(dim(coclust_ij),collapse=" ")) }
			   cc[boot_cells, boot_cells] = 
					cc[boot_cells, boot_cells] + 
					coclust_ij
			   if (verbose) { message("will add trials") }
			   cosamp[boot_cells,boot_cells] = 
					cosamp[boot_cells, boot_cells] + 1
			}
		}
	} else {
		cc = res[[1]]$coclust
		cosamp = res[[1]]$num_trials
		for(co in res[-1]) {
			cc = cc + co$coclust
			cosamp = cosamp + co$num_trials
		}
	}
	
	if (sum(cosamp == 0) > 1) {
		stop("Bootstrap sampling is not deep enough - consider increasing the number of bootstrap iterations and/or the fraction of cells sampled in each iteration.")
	}
	cc = cc/cosamp
#	scl = scc_coclust_to_cell_modules(scdb, scl, cc, k_cut, min_clust_size) 
	
#	scl@alg_params = c(scl@alg_params, list(bootstrap = T, bootstrap_N = N_boot, bootstrap_cut = k_cut, min_clust_size = min_clust_size))
	
	.restore_seed(oldseed)
	
#	if (is.null(scdb)) {
#		save(cc, cor_cc, scl, file="coclust.Rda")
#	}
#	else {
#		scdb_save(scdb, cc, file = "coclust.Rda")
#	}
	
	return(cc)
}

#' Update cell modules by a cut tree on a cells' co-clust matrix
#'
#' @param scdb 
#' @param scl 
#' @param cc 
#' @param k_cut 
#' @param min_clust_size 
#'
#' @return
#' @export
#'
#' @examples
scc_coclust_to_cell_modules = function(scdb, scl, cc, k_cut, min_clust_size, method="hclust") 
{
	stopifnot(!is.null(scdb) || (!is.null(scl) && !is.null(cc)))
	
	if (is.null(scl)) {
		sc_cl = scdb_load(scdb, "cell_modules")
		if (is.null(scl)) {
			stop("cell modules not supplied and not found in scdb.")
		}
	}	
	
	if (is.null(cc))
		cc = scdb_load(scdb, "coclust")
		
		if (is.null(cc)) {
			stop("co-clust matrix not supplied and not found in scdb.")
	}

	if (method == "hclust") {
		cor_cc = .cor(cc)
		diag(cor_cc) = 0
		cor_cc[is.na(cor_cc)] = 0
	
		cor_cc2 = .cor(cor_cc)
		cor_cc2[is.na(cor_cc2)] = 0
		hc = hclust(as.dist(1-cor_cc2))
		clusts = cutree(hc, k_cut)
	}
	else if (method == "louvain") {
		g_cc = graph_from_adjacency_matrix(cc, mode="undir", weighted=T, diag=F)
		g_comm = cluster_louvain(g_cc)
		clusts = membership(g_comm)
	}
	else {
		stop(sprintf("Unknown coclust clustering method (%s), supporting hclust or louvain", method))
	}
	
	csz = table(clusts)
	my_size = csz[clusts]
	
	outlier_cell = NULL
	if (any(csz < min_clust_size)) {
		clusts[my_size < min_clust_size] = max(clusts)+1
		outlier_cell = names(clusts)[clusts == max(clusts)][1]
	}
	
	clusts1 = as.integer(as.factor(clusts))
	cl_nms = colnames(scl@scmat@mat) #names(clusts)
	names(clusts1) = cl_nms
	clusts= clusts1
	message("will reorder")
	clusts = .reorder_knn_clusts(scl@feat_mat, clusts)
	names(clusts) = cl_nms
	scl@clusts = clusts
	
	# plot cc matrix
	png(sprintf("%s/bootstrap_coclust_%s%s.png", get_param("outdir"), method, ifelse(method == "hclust", sprintf("_cut%d_minSize%d", k_cut, min_clust_size), "")), width=min(length(clusts), 3000), height=min(length(clusts), 3000))
	cord = order(clusts)
	image(cc[cord, cord], zlim=c(0, 1), col=colorRampPalette(c("white", "darkred", "purple", "black"))(1000), xaxt='n', yaxt='n', xlab="", ylab="")
	clines = cumsum(table(clusts)) / length(clusts)
	
	abline(h=clines, lwd=3)
	abline(v=clines, lwd=3)
	
	at = clines - table(clusts) / (2 * length(clusts))
	if (is.null(outlier_cell)) {
		ccols = 'black'
	}
	else {
		out_cl = clusts[outlier_cell]
		message("outlier clust is: ", out_cl)
		ccols = ifelse(1:max(clusts) == out_cl, 'red', 'black')
	}
	mtext(1:max(clusts), 2, at = at, cex=3, col=ccols)
	mtext(1:max(clusts), 3, at = at, cex=3, col=ccols)
	dev.off()
	
	message("will postproc clusts")
	scl = .scc_postprocess(scl)

	
	return(scl)
}

scc_cluster_hclust = function(scmat, features, K=NA)
{
	#TBA
	stop("ERROR: hclust not supported yet.")
	return(NA)
}

#' Output genes cell modules footprint matrix with metadata on the cell modules
#'
#' @param sc_cl 
#' @param tab_clust_fp_fn 
#' @param min_max_clust_fp_to_report 
#' @param clust_fp_metadata_fields 
#'
#' @return
#' @export
#'
#' @examples
scc_write_clust_fp_tab = function(sc_cl, 
																	outdir = get_param("outdir"),
																	tab_clust_fp_fn = get_param("scc_clust_fp_fn"),
																	min_max_clust_fp_to_report = get_param("scc_min_max_clust_fp_to_report"),
																	clust_fp_metadata_fields = get_param("scc_clust_fp_metadata_fields")) 
{
	f = apply(sc_cl@clust_fp, 1, max) > min_max_clust_fp_to_report
	
	# add #cells per clust and mean #umis (~ cell size)
	out_df = rbind(tapply(colSums(sc_cl@scmat@mat), sc_cl@clusts, mean), table(sc_cl@clusts))
	rownames(out_df) = c('mean_umis', 'n_cells')
	
	# add required breakdown to features
	if (!is.null(clust_fp_metadata_fields)) {
		for (s in clust_fp_metadata_fields) {
			out_df = rbind(table(sc_cl@scmat@cell_metadata[names(sc_cl@clusts), s], sc_cl@clusts), out_df)
		}
	}
	
	# actual clust_fp
	out_df = rbind(out_df, round(log2(sc_cl@clust_fp[f,]), 2))
	
	write.table(out_df, paste0(outdir, "/", tab_clust_fp_fn),	sep = "\t",quote = F)

}

#
#' Extract features mat
#'
#' @param scmat A tgScMat expression matrix.
#' @param features Either  a matrix of cell features (gene expression, PCA principle componant, etc.),
#' or a list of markers to be extracted from scmat. Clustering will be performed on these features.
#' If left blank, markers will be selected using the default method.
#' @param factor_by_median 
#' @param norm_size 
#' @param k_scale_umi 
#' @param norm_features 
#' 
#' @return A matrix with one row per marker gene, one column per cell.
#' Gene expression is normalized according to the parameters.
#'
#' @export
create_feature_mat = function(scmat,
															features = NULL,
															factor_by_median = get_param("scc_factor_by_median"),
															norm_size = get_param("scc_norm_size"),
															k_scale_umi = get_param("scc_k_scale_umi"),
															norm_features = T) {
	
	feat_mat = .extract_feat_mat(scmat, features)
	if (norm_features) {
		feat_mat = .normalize_feat(scmat, feat_mat,
															 factor_by_median = factor_by_median, norm_size = norm_size,
															 k_scale_umi = k_scale_umi)
	}

	return(Matrix(feat_mat))
}

#
# Normalize features matrix
#'
#' @param scmat 
#' @param feat_mat 
#' @param factor_by_median 
#' @param norm_size 
#' @param k_scale_umi 
#'
#' @return
.normalize_feat = function(scmat, feat_mat,
				factor_by_median = get_param("scc_factor_by_median"),
				norm_size = get_param("scc_norm_size"), k_scale_umi = get_param("scc_k_scale_umi"))
{
	umis = feat_mat
	n_median = 1
	if(factor_by_median) {
		n_median = median(colSums(umis))
		# cat("..use median norm ", n_median)
	}
	if(norm_size) {
		# cat(", norm size ")
		umis = n_median * t(t(umis)/(1+colSums(umis)))
	}
	if(scmat@stat_type == "umi") {
		# if(is.na(k_scale_umi)) {
		# 	#find a better heuristic
		# 	k_scale_umi = 7
		# }
		cell_fps = log2(1+k_scale_umi*umis)
	} else {
		if(min(umis) > 0 & max(umis)/min(umis)>100) {
			cell_fps = log2(umis)
		} else {
			cell_fps = umis
		}
	}
	return(cell_fps)
}

#' Extract features matrix
#'
#' @param scmat 
#' @param features 
#'
#' @return If features is null, select marker genes automatically. If features is a vector of genes, extract sub-matrix of those genes. Otherwise expecting features to be a matrix sized k (features) on n (cells), and return it. 
.extract_feat_mat = function(scmat, features = NULL)
{
	#
	# extract expression for markers
	#
	if (is.null(features)) {
		# select features
		cat("Extracting markers using default parameters.\n")
		feat_mat = as.matrix(scmat@mat[select_markers(scmat),])

	} else if (is.null(dim(features)) && length(intersect(features, scmat@genes) > 0)) {

		# got markers
		marks = intersect(features, scmat@genes)
		if (length(marks) != length(features)) {
			warning ("Used only ", length(marks), " out of ", length(features) , " features.")
		}
		feat_mat = as.matrix(scmat@mat[marks,])
	} else if (length(dim(features)) == 2 && dim(features)[2] == scmat@ncells) {

		# got 2D feature matrix
		feat_mat = as.matrix(features)
	} else {
		stop("invalid features used for clustering")
	}
	return(feat_mat)
}

####


#' @export
#'
setMethod(
	"show",
	signature = "tgScMatClust",
	definition =
	 function(object) {
		cat("tgScMatClust, ", length(object@clusts), " cells on ",
				length(unique(object@clusts)), " clusts\n")
	 	cat("Parameters: ")
	 	print(object@alg_params)
		invisible(NULL)
	 }
)

#'
#' Export a clustering results to file.
#'
#' @param file Prefix of file names for outputting.
#'
#' @export
#'
# setGeneric("tgscm_export",
	# function(.Object, file,...) stnadrdGeneric("tgscm_export"))
setMethod(
	"tgscm_export",
	signature = "tgScMatClust",
	definition =
	 function(.Object, file, supress_mat=F, ...) {
		if(!supress_mat) {
		  tgscm_export(.Object@scmat, file)
		}
		write.table(as.matrix(.Object@feat_mat), file=sprintf("%s.feat", file), quote=F, sep="\t")
		write.table(.Object@clusts, file=sprintf("%s.clust", file), quote=F, sep="\t")
		write.table(.Object@clust_fp, file=sprintf("%s.clust_fp", file), quote=F, sep="\t")
		write.table(.Object@alg_params, file=sprintf("%s.alg", file), quote=F, sep="\t")
	 }
)

#' Read a clustering result from file.
#'
#' @param file Name of the file for inputting.
#'
#' @param scmat an instance of the relevant \code{\linkS4class{tgScMat}}, optional.
#'
#' @export
#'
#'
scc_import_clust = function(file, scmat=NULL) {
	if(is.null(scmat)) {
		scmat = scm_import_mat(file)
	}
  .Object = tgScMatClust()
	.Object@scmat = scmat
	.Object@feat_mat = Matrix(as.matrix(read.table(file=sprintf("%s.feat", file), head=T, sep="\t")))
	.Object@clusts = read.table(file=sprintf("%s.clust", file), sep="\t")
	.Object@clusts = data.matrix(.Object@clusts)[,1] # convert to array
	.Object@clust_fp = as.matrix(read.table(file=sprintf("%s.clust_fp", file), sep="\t", check.names=F))
	.Object@alg_params = as.list(read.table(file=sprintf("%s.alg", file), sep="\t", stringsAsFactors = F))
	.Object@nclust = length(unique(.Object@clusts))
	return(.Object)
}


# 
#' Create a footprint on covered genes
#'
#' @param scm_cl 
#' @param ref_func 
#' @param min_total_umi 
#'
#' @return footprint per cluster
#' 
.calc_clusts_fp = function(scm_cl, ref_func = median, min_total_umi = get_param("scc_min_total_gene_count_for_footprint"))
{
	us = scm_cl@scmat@mat
	f_g_cov = rowSums(us) > min_total_umi

	clust_geomean = .row_stats_by_factor(us[f_g_cov,], scm_cl@clusts, function(y) {exp(rowMeans(log(1+y)))-1})
	# clust_geomean =  t(apply(us[f_g_cov,], 1,  #slow!!
	# 		function(x) tapply(x,
	# 			scm_cl@clusts,
	# 			function(y) exp(mean(log(1+y)))-1)))
	clust_meansize = tapply(colSums(us), scm_cl@clusts, mean)
	ideal_cell_size = pmin(1000, median(clust_meansize))
	g_fp = t(ideal_cell_size*t(clust_geomean)/as.vector(clust_meansize))
	#normalize each gene
	fp_reg = 0.1   #0.1 is defined here because 0.1*mean_num_of_cells_in_cluster 
		       #is epxected to be 3-7, which means that we regulairze 
		       #umicount in the cluster by 3-7.
	g_fp_n = (fp_reg+g_fp)/apply(fp_reg+g_fp, 1, ref_func)

	return(g_fp_n)
}

#' Re-order clusters by hclust centers of clusters
#'
#' @param data 
#' @param clusts 
#'
#' @return reordered clusts vector
#' 
.reorder_knn_clusts = function(data, clusts) {
	k = length(unique(clusts))
	# # reorder clusters
	centers = .row_stats_by_factor(data, clusts)
	centers_hc = hclust(dist(cor(centers)), "ward.D2")
	tomap = rep(NA, k)
	tomap[centers_hc$order] = 1:k
	onames = names(clusts)
	clusts = tomap[clusts]
	names(clusts) = onames
	return(clusts)
}


#'
#' Calculate ordered knn matrix. 
#' 
#'  IMPORTANT: calculated matrix might have more rows than expected! always use only the first relevant rows,
#'  if k is small, will calculate more than requested, since the knn object is used in more than 1 place.
#'
#' @param scl 
#' @param k_knn 
#'
#' @return a k_knn X ncells matrix,  column i contains cell i's k nearest neighbors,
#'  ordered from closest to farthest.
#'  
.scc_calc_ordered_knn = function(scl, k_knn = get_param("scc_k_knn"))
{
	message("Calling calc order knn in cluster.r")

	if(length(scl@knn_ordered) > 1 && # knn_ordered exists
		 nrow(scl@knn_ordered) >= k_knn) {
		message("using cached knn ord")
		# no need to calculate again!
		return(scl)
	}
	sim_mat = .calc_sim_mat(scl)
	scl@knn_ordered = tg_simmat_to_knn_ordered(sim_mat, k_knn)
	return(scl)
}


#' Create cell's knn matrix
#'
#' @param scl 
#' @param k_knn 
#'
#' @return
.scc_comp_knn_matrix = function(scl, k_knn = get_param("scc_k_knn"))
{
# calc knn_ordered and save in object
	scl = .scc_calc_ordered_knn(scl, k_knn)
# initialize knn graph adj matrix
	message("will comp knn adjs from knn ord")
	scl@m_knn = tg_knn_ordered_to_knn_adjs(scl@knn_ordered, k_knn)
	message("done comp knn adjs from knn ord")
	scl@knn_is_balanced = F
	return(scl)
}

#' Calculate balanced knn matrix
#'
#' @param scl 
#' @param k_knn 
#' @param k_expand 
#'
#' @return tgScMatClust object
#' 
.scc_comp_balanced_knn_matrix = function(scl, k_knn = get_param("scc_k_knn"), k_expand = get_param("scc_k_expand_for_balanced"))
{
	cells = scl@scmat@cells
	ncells = length(scl@scmat@cells)
	k_knn_e = min(k_knn * k_expand, ncells-1)

# calc knn_ordered and save in object
	scl = .scc_calc_ordered_knn(scl, k_knn_e)
# initialize knn graph adj matrix
	scl@m_knn = tg_knn_ordered_to_knn_balanced_adjs(scl@knn_ordered, k_knn)
	scl@knn_is_balanced = T
	return(scl)
}

#' Calculate features similarity matrix (person correlation)
#'
#' @param scl 
#'
#' @return cor matrix 
#' 
.calc_sim_mat = function(scl) 
{
	message("Will compute cor on feat mat")

	m = as.matrix(.cor(as.matrix(scl@feat_mat)))
	#m = as.matrix(cor(as.matrix(scl@feat_mat)))
	
	return(m)
}

#' Filter bad clusters. Either clusters with out high enough fold-change genes, and/or clusters with low self-confusion neighbours
#'
#' @param scl 
#' @param min_max_fold_change 
#' @param confu 
#' @param k_knn_for_confu 
#' @param min_self_confusion 
#'
#' @return tgScMatClust object with updated good clusters
#' @export
#'
#' @examples
scc_filter_bad_clusters = function(scl, 
				min_max_fold_change= get_param("scc_min_max_clust_fp_to_filter"), 
				confu = NULL,
				k_knn_for_confu = get_param("scc_k_knn_for_confu"),
				min_self_confusion = get_param("scc_min_self_confusion"))
{
	if(is.null(min_max_fold_change) & is.null(min_self_confusion)) {
		stop("at least one threshold (fold change or self confusion) need to be determined in scc_filter_bad_clusters")
	}
	if(!is.null(min_self_confusion) & is.null(confu)) {
		if(is.null(k_knn_for_confu)) {
			stop("in scc_filter_bad_clusters, either provide a confusion matrix or the k parameter for computing it, or avoid setting the confusion threshold")
		}
		scl = .scc_comp_balanced_knn_matrix(scl, k_knn_for_confu, get_param("scc_k_expand_for_balanced"))
		mknnio_ij = scl@m_knn
		isclust_ci = diag(max(scl@clusts))[,scl@clusts]
		confu = isclust_ci %*% (mknnio_ij>0) %*% t(isclust_ci)
	}
	max_fc = apply(scl@clust_fp, 2, max)

	bad_clusts = c()
	if(!is.null(min_max_fold_change)) {
		bad_clusts = which(max_fc < min_max_fold_change)
		message("removing ", length(bad_clusts), " clusters due to low fold change in their gene expression")
	}
	if(!is.null(min_self_confusion)) {
		self_confu = diag(confu)/rowSums(confu)
		bad_clusts1 = which(self_confu < min_self_confusion)
		bad_clusts = union(bad_clusts, bad_clusts1)
		message("removing ", length(bad_clusts1), " clusters due to low self confusion, and in total ditching ", length(bad_clusts), " clusters")
	}
	
	if(length(bad_clusts) == 0) {
		return(scl)
	}
	bad_cells = names(scl@clusts)[which(scl@clusts %in% bad_clusts)]

	good_cells = setdiff(colnames(scl@scmat@mat), bad_cells)

	scl@scmat = scm_sub_mat(scl@scmat, cells=good_cells)
	scl@feat_mat = scl@feat_mat[,good_cells]
	scl@clusts = as.integer(as.factor(scl@clusts[good_cells]))
	names(scl@clusts) = good_cells
        scl@clusts = .reorder_knn_clusts(scl@feat_mat, scl@clusts)
	scl = .scc_postprocess(scl)
	scl@nclust = scl@nclust - length(bad_clusts)
        if(length(scl@knn_ordered) > 1) {
		scl@knn_ordered = scl@knn_ordered[good_cells,]
	}
	return(scl)
}



#' TODO: go over and remove/generelize specific params
#'
#' @param scl 
#' @param k_knn 
#' @param min_genemod_fold 
#' @param min_confu_score 
#' @param plot_dir 
#'
#' @return
#'
#' @examples
scc_orphan_unlinked_cells = function(
	scl,
	k_knn = get_param("scc_k_knn"),
	min_genemod_fold = 1.5,
	min_confu_score = get_param("scc_min_confusion_score"),
	plot_dir = get_param("outdir")
)
{
	
	if (is.null(k_knn)) {
		k_knn = min(max(scl@scmat@ncells / 40, 75, 500))
	}
	
	#compute confusion
	scl = .scc_comp_balanced_knn_matrix(scl, k_knn, 10)
	mknnio_ij = scl@m_knn
	isclust_ci = diag(max(scl@clusts))[,scl@clusts]
	confu_ic = as.matrix((mknnio_ij>0) %*% t(isclust_ci))
	confu_ci = as.matrix(isclust_ci %*% (mknnio_ij>0) )
	
#what should we filter? i with low ci
	n = nrow(confu_ic)
	confubal_ic = confu_ic * t(confu_ci)
	confu_w = confubal_ic[n*(scl@clusts-1) + 1:n]
	confu_tot = rowSums(confubal_ic)
	#score must be normalize since each cell has a different number of outgoing edges
	confu_score = confu_w/(confu_tot+k_knn)

	T_mark = 2
	ismark_gc = (scl@clust_fp > T_mark)	
	rna = as.matrix(scl@scmat@mat[rownames(scl@clust_fp),])
	rna = t(t(rna)/colSums(rna))*1000
	lumi_gi = log2(1+7*rna)
	
	totmark_ic = t(lumi_gi) %*% ismark_gc
	totmark_ic = t(t(totmark_ic)/apply(totmark_ic,2,median))
	fpenrich = totmark_ic[n*(scl@clusts-1) + 1:n]

	#also compute score for each e
	f_to_filt = fpenrich < min_genemod_fold & confu_score < min_confu_score

	png(sprintf("%s/orphan_params.png", plot_dir), w=800, h=800)
	plot(log2(0.05+fpenrich), confu_score, pch=19, cex=0.6)
	points(log2(0.05+fpenrich[f_to_filt]), confu_score[f_to_filt], pch=19, cex=0.6, col="red")
	dev.off()
	png(sprintf("%s/orphan_on_clusts.png", plot_dir), w=1200, h=600)
	orph_stat = table(scl@clusts, f_to_filt)
	barplot(t(orph_stat/rowSums(orph_stat)),las=2)
	dev.off()

	nms_to_filt = colnames(lumi_gi)[f_to_filt]
	scl = scc_move_cells_to_orphan(scl, bad_cells = nms_to_filt)
	
	return(scl)
}

#' Find cells that are outliers in the current clusters
#'
#' @param scl 
#' @param T_lfc 
#' @param min_outlier_u 
#'
#' @return vector of outlier cells
#' @export
#'
#' @examples
scc_find_clustering_outliers = function(scl, 
																		T_lfc = get_param("scc_clust_outliers_t_lfc"), 
																		min_outlier_u = get_param("scc_clust_outliers_min_gene_count"))
{
	# 1 encoding of the clusters, rows= clusters, cols = cells
	isclust_ci = diag(max(scl@clusts))[,scl@clusts]
	
	#we compute the expected number of umi per gene per cell given clust
	u_gi = scl@scmat@mat

	u_i = colSums(as.matrix(u_gi))

	ishigh_g = apply(u_gi, 1, max)>= min_outlier_u

	u_gi = u_gi[ishigh_g,]

	u_gc = as.matrix(u_gi %*% t(isclust_ci))
	u_c = colSums(u_gc)

	p_gc = t(t(u_gc) / u_c)
	p_gi = p_gc %*% isclust_ci
	exp_gi = t(t(p_gi) * u_i)

	#we estimate a premissive sd on that
	lfc_gi = log2((1+u_gi)/(1+exp_gi))

#	sd_gi = pmax(sqrt(ceiling(exp_gi)),1)
#	z_gi = (u_gi - exp_gi)/sd_gi

	maxlfc_g = apply(lfc_gi, 1, max)
	maxlfc_i = apply(lfc_gi, 2, max)

	if(sum(maxlfc_g > T_lfc) > 1 & sum(maxlfc_i > T_lfc) > 2) {
		outu_gi = log2(1 + as.matrix(u_gi[maxlfc_g > T_lfc, maxlfc_i > T_lfc]))

	#reporting the outlier gene / cell matrix

		hc1 = hclust(dist(cor(outu_gi)), "ward.D2")
		hc2 = hclust(dist(cor(t(outu_gi))), "ward.D2")
		png(sprintf("%s/outlier_mat.png", get_param("outdir")), w=2000, h=2000)
		shades = colorRampPalette(c("white", "blue", "red", "yellow", "black"))(1000)
		image(t(outu_gi[hc2$order, hc1$order]), col=shades, xaxt='n', yaxt='n')
		mtext(rownames(outu_gi)[hc2$order], at=seq(0,1,l=length(hc2$order)), las=2, side=2, cex=0.8)
		mtext(rownames(outu_gi)[hc2$order], at=seq(0,1,l=length(hc2$order)), las=2, side=4, cex=0.8)
		dev.off()
	}

	return(names(which(maxlfc_i > T_lfc)))
}

#' Clean ambient noise on clusters. TODO: move more hard-coded to params?
#'
#' @param scl 
#' @param orig_mat 
#' @param epsilon 
#' @param batch_attr 
#' @param T_zclean 
#'
#' @return
#'
.scc_clean_ambient_on_clusts = function(scl,
																				orig_mat,
																				epsilon = get_param("scm_amb_noise_epsilon"),
																				batch_attr = get_param("scm_batch_meta_attr"),
																				rpt_base_dir = sprintf("%s/ambient", get_param("outdir")),
																				T_zclean = get_param("scc_filter_ambient_Z_max"))
{
	message("will clean ambient on clusts")
	batch_factor = as.integer(as.factor(scl@scmat@cell_metadata[names(scl@clusts),batch_attr]))

	isclust_ci = diag(max(scl@clusts))[,scl@clusts]    # 1 encoding of the clusters, rows= clusters, cols = cells
	if(max(batch_factor)==1) {
		isbatch_bi = matrix(rep(1,length(batch_factor)),nrow=1)
	} else {
		isbatch_bi = diag(max(batch_factor))[,batch_factor]    # 1 encoding of the clusters, rows= clusters, cols =nodes
	}
	#we break clusts on batches
	n_cb = table(scl@clusts, batch_factor)
	n_b = colSums(n_cb)
	n_c = as.vector(table(scl@clusts))
	n = sum(n_c)

	#we compute total U_{g,b} - number of umis per gene per batch
	u_gi = orig_mat[,names(scl@clusts)]

#	u_gb = t(apply(u_gi, 1, function(x) tapply(x, batch_factor, sum)))
	u_gb = u_gi %*% t(isbatch_bi)
	amb_gb = t(t(u_gb) * (epsilon/n_b))

#f_g just to save time and
	f_g = apply(amb_gb, 1, max) > 0.02

	u_gi = u_gi[f_g,]
	u_gb = u_gb[f_g,]
	amb_gb = amb_gb[f_g,]

	Tnotamb_gb = ceiling(ceiling(amb_gb)+3*sqrt(ceiling(amb_gb)))
	Tnotamb_gi = Tnotamb_gb[,batch_factor]
	umaxamb_gi = pmin(as.matrix(u_gi), as.matrix(Tnotamb_gi))
#	Tnotamb_f = u_gi > Tnotamb_gi

	#we compute amb_{g,c} and obs_{g,c} as the expected number of
	#ambient umi's per cluster, and the observed number
	#the observed number is scissoring high umi count cells to avoid
	#outlier take over
	amb_gc = as.matrix(amb_gb) %*% t(as.matrix(n_cb))

	obs_gc = umaxamb_gi %*% t(isclust_ci)
#	obs_gc = t(apply(umaxamb_gi, 1, function(x) tapply(x, scl@clusts, sum)))

	d_gc = obs_gc - amb_gc
	z_gc = d_gc / sqrt(amb_gc)

#	save(Tnotamb_gb, Tnotamb_gi, umaxamb_gi, amb_gc, obs_gc, d_gc, z_gc, "dbg_amb.Rda")

#filtered genes are those for which:
# clusters on z_gc>T_Z capture more than 1-k*epsilon of the umi
# at least 1/4 of the cells are in cells with Z<2
	u_gc = u_gi %*% t(isclust_ci)
#	u_gc =  t(apply(u_gi, 1, function(x) tapply(x, scl@clusts, sum)))
	upos_gc = u_gc * (z_gc < T_zclean)
	ambumitot_g = rowSums(upos_gc)/rowSums(u_gc)
	ambcelltot_g = (z_gc < T_zclean) %*% as.matrix(n_c)

	tofilt_g = (ambcelltot_g > (n/8)) & (ambumitot_g < epsilon*2)

# for these case we clean umis below Tnotamb_gi in clusters with Z<1
	z_gi = z_gc[,scl@clusts]
	z_gi[!tofilt_g,] = 8
	cleanu_gi = u_gi
	cleanu_gi[z_gi<T_zclean & u_gi < Tnotamb_gi] = 0
	full_gi = orig_mat[,names(scl@clusts)]
	denoise_gi = full_gi
	denoise_gi[f_g,] = cleanu_gi

	all_emp = c()
	denoise_nms = names(which((rowSums(full_gi) - rowSums(denoise_gi))>200))
	for(nm in denoise_nms) {
		f = n_c>10 & z_gc[nm,]<T_zclean
		f_out = n_c>10 & z_gc[nm,]>T_zclean
		p_back = median((3+u_gc[nm,f])/n_c[f])
		empiric = p_back/mean(orig_mat[nm,])
		all_emp = c(all_emp, empiric)
		all_gc= tapply(full_gi[nm,], scl@clusts, sum)
		min_signal = min((all_gc[f_out]+3)/n_c[f_out])
		d = density(log2((3+all_gc)/n_c))
		if(dir.exists(rpt_base_dir)) {
			png(sprintf("%s/%s.png", rpt_base_dir, nm), w=600,h=600)
			plot(d$x, d$y, type="b",lwd=3, col=ifelse(d$x<log2(min_signal),"black", "red"), pch=19, cex=0.5)
			abline(v=log2(p_back), lwd=3)
			dev.off()
		}
	}
	if(length(all_emp) > 0) {
		names(all_emp) = denoise_nms
	}
	save(denoise_gi,
		u_gi,
		cleanu_gi,
		z_gi,
		tofilt_g,
		ambcelltot_g,
		ambumitot_g,
		upos_gc,
		u_gc, z_gc, d_gc,
		amb_gc, obs_gc,
		all_emp,
		Tnotamb_gb, Tnotamb_gi, umaxamb_gi, orig_mat,
		u_gb,
		file=sprintf("%s/dbg_celanamb.Rda", rpt_base_dir))

	write.table(data.frame(eps=all_emp, 
				n = rowSums(orig_mat)[names(all_emp)], 
				clean_n =(rowSums(full_gi) - rowSums(denoise_gi))[names(all_emp)]), 
				sprintf("%s/empiric_eps.txt", rpt_base_dir),
				quote=F, sep="\t")
	return(denoise_gi)
}

#' Title
#'
#' @param scl 
#' @param orig_mat 
#' @param batch_attr 
#' @param n_reg 
#'
#' @return observed and expected number of umis per gene and batch
#'
.scc_find_batchy_genes_on_clusts = function(scl,
																						orig_mat,
																						batch_attr = get_param("scm_batch_meta_attr"),
																						n_reg = 10
)
{
	message("will search for batchy genes on clusters")
	batch_factor = as.integer(as.factor(scl@scmat@cell_metadata[names(scl@clusts),batch_attr]))

	isclust_ci = diag(max(scl@clusts))[,scl@clusts]    # 1 encoding of the clusters, rows= clusters, cols = cells
	isbatch_bi = diag(max(batch_factor))[,batch_factor]    # 1 encoding of the batch, rows= batch, cols =nodes

	#number of cells in cluster,batch
	n_cb = table(scl@clusts, batch_factor)
	n_b = colSums(n_cb)
	n_c = as.vector(table(scl@clusts))
	n = sum(n_c)

	u_gi = orig_mat[,names(scl@clusts)]

	#we compute total U_{g,b} - number of umis per gene per batch
	n_gb = u_gi %*% t(isbatch_bi)

	n_gc = as.matrix(u_gi %*% t(isclust_ci))
	
	#now we need to find the expected number of umi per gene and batch e_gb\	
	e_gb = n_gc %*% (as.matrix(n_cb)/n_c)
	
	return(list(ratio=(n_gb+n_reg)/(e_gb+n_reg), e=e_gb, o=n_gb))
}

scc_get_fpcor_genes = function(scl, 
		filt_lateral_genes_supervised, 
		cor_thresh,
		filt_focus_genes_supervised = NULL,
		fig_fn = "supervised_lateral.png")
{
	fp = log2(scl@clust_fp)

	g_lateral = intersect(rownames(fp), filt_lateral_genes_supervised)
	g_focus = intersect(rownames(fp), filt_focus_genes_supervised)

	all_cor = rep(0, nrow(fp))
	for(nm in g_lateral) {
		y = fp[nm,]
		y_cor = apply(fp, 1, function(x) { cor(y, x, m="spearman") })
		all_cor = pmax(y_cor, all_cor)
	}
	for(nm in g_focus) {
		y = fp[nm,]
		y_cor = apply(fp, 1, function(x) { cor(y, x, m="spearman") })
		all_cor = pmin(ifelse(y_cor > cor_thresh, 0, 1),  all_cor)
	}

	to_filt = rownames(fp)[which(all_cor>cor_thresh)]
	png(fig_fn, w=200+12*length(to_filt), h=200+12*length(to_filt))
	cr = cor(t(fp[to_filt,]))
	hc = hclust(dist(cr), "ward.D2")
	shades = colorRampPalette(c("blue", "white", "yellow"))(200)
	image(cr[hc$order, hc$order], zlim=c(-1,1), col=shades, xaxt='n', yaxt='n')
	mtext(to_filt[hc$order], at = seq(0,1,l=length(to_filt)), side =1, las=1, cex=0.8)
	mtext(to_filt[hc$order], at = seq(0,1,l=length(to_filt)), side =2, las=1, cex=0.8)
	dev.off()

	return(rownames(fp)[which(all_cor>cor_thresh)])
}
