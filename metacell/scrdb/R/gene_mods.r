
#' Gene modules (groups of correlated genes)
#
#' Supporting algorithms for detection of modules and some stats on them
#'
#' @slot scmat tgScMat. 
#' @slot foc_genes vector. selected genes
#' @slot feat_mat matrix. 
#' @slot gene_cor matrix. 
#' @slot rnd_gene_cor matrix. 
#' @slot hc_order vector. 
#' @slot gmods vector. 
#' @slot indep_genes vector. 
#' @slot gmod_mean matrix.
#' @slot gmod_fp matrix. 
#' @slot alg_params list. 
#' @slot nclust integer. 
#' @slot knn_ordered matrix. k_knn X ncells matrix,  column i contains cell i's k nearest neighbors, ordered from closest to farthest.
#' @slot m_knn matrix. 
#' @slot rseed numeric. 
#'
#' @importClassesFrom Matrix dgCMatrix
#'
#### @export tgScMatClust
#### @exportClass tgScMatClust
#
#### @include scmat.r
tgScGeneMods <- setClass(
  "tgScGeneMods",
  slots = c(
    scmat = "tgScMat",
    foc_genes = "vector",
    feat_mat = "matrix",
    gene_cor = "matrix",
    rnd_gene_cor = "matrix",
    hc_order = "vector",
    gmods = "vector",
    indep_genes = "vector",
    gmod_mean = "matrix",
    gmod_fp = "matrix",
    alg_params = "list",
    nclust = "integer",
    knn_ordered = "matrix",
    m_knn = "matrix",
    rseed = "numeric")
)

## Houskeeping staf 
####
#' @export
#'
setMethod(
	"show",
	signature = "tgScGeneMods",
	definition =
	 function(object) {
		cat("tgScGeneMods, ", length(object@gmods), " genes on ",
				length(unique(object@gmods)), " mods\n")
		invisible(NULL)
	 }
)

#'
#' Export a gene mod object to files
#'
#' @param file Prefix of file names for outputting.
#'
#' @export
#'
setMethod(
	"tgscm_export",
	signature = "tgScGeneMods",
	definition =
	 function(.Object, file, supress_mat=F, ...) {
		if(!supress_mat) {
		  tgscm_export(.Object@scmat, file)
		}
		stop("export gmods not yet implemented")
	 }
)

#' Read a gmod object from file
#'
#' @param file Name of the file for inputting.
#'
#' @param scmat an instance of the relevant \code{\linkS4class{tgScMat}}, optional.
#'
#' @export
#'
#'
scgm_import = function(file, scmat=NULL) {
	if(is.null(scmat)) {
		scmat = scm_import_mat(file)
	}
	stop("import gmods not yet implemented (not a bad idea to implement)")
	return(.Object)
}

### ================

scgm_init_cclust_fp = function(scl, foc_genes = NULL, min_fold=1) 
{
	if(is.null(foc_genes)) {
		foc_genes = names(which(
				apply(abs(log2(scl@clust_fp)), 1, max) > min_fold))
	}
	if(length(foc_genes) < 10) {
		stop("trying to generate gene modules with less only ", length(foc_genes), " (i.e. <10) candidates - not recommended")
	}

	scgm = new("tgScGeneMods")
	scgm@scmat = scl@scmat
	scgm@foc_genes = foc_genes
	scgm@feat_mat = log2(scl@clust_fp[foc_genes,])
	return(scgm)
}

scgm_init_cclust_resid = function(scl, foc_genes = NULL, min_tot = 200) 
{
	scgm = new("tgScGeneMods")
	scgm@scmat = scl@scmat
	if(is.null(foc_genes)) {
		foc_genes = names(which(rowSums(scl@scmat@mat) > min_tot))
	}
	scgm@foc_genes = foc_genes

	scgm@feat_mat = scc_calc_clust_resid_exp(scl, foc_genes)

	return(scgm)
}

scgm_build_cor_metric = function(scgm, cor_method="pearson")
{
	scgm@gene_cor = .cor(t(scgm@feat_mat), spearman=cor_method == "spearman") #this is it?
	diag(scgm@gene_cor) = 0
	
  #permutating the matrix and recompute?
	rmat = t(apply(t(scgm@feat_mat),2,sample))
	scgm@rnd_gene_cor = .cor(t(rmat), spearman=cor_method == "spearman") #this is it?
	diag(scgm@rnd_gene_cor) = 0
	return(scgm)
}

.scgm_postprocess = function(scgm) 
{
	message("will comp fp")
	scgm@gmod_mean = .scgm_calc_gmod_mean(scgm)
	scgm@gmod_fp = .scgm_calc_gmod_fp(scgm)
	scgm@nclust = length(unique(scgm@gmods))
	return(scgm)
}

#' @export
scgm_from_lists = function(scmat, gmods, cor_method="pearson") 
{
	scgm = new("tgScGeneMods")
	scgm@scmat = scmat

	genes = names(gmods)
	if(length(setdiff(genes, rownames(scmat@mat)))>0) {
		stop("trying to initialize gmod from a list with mismatched gene name, e.g. ", head(setdiff(genes, rownames(scmat@mat)),5),")")
	}
	scgm@feat_mat = scmat@mat[genes,]

	scgm = scgm_build_cor_metric(scgm, cor_method)

	scgm@gmods = gmods

	scgm = .scgm_postprocess(scgm)

	return(scgm)

}

#' Given a cell cluster solution, this will cluster genes given either their
#' cluster enrichment per cluster, or given the deviation of their expression
#' from the cluster mean predictions (residual expression)
#'
#' @param scl A tgScMatClust cell cluster object
#'
#' @param foc_genes - genes for analysis (i.e. only genes with enough stat)
#'
#' @param feat_type - "clust_fp" is using cluster enrichment, "resid" for residuals
#' @param use_hclust - if this is false, use knn_clust
#' @param k_clusts number of clusters if using hclust
#' @param k_knn number of clusters if using knn clustering
#'
#' @return A tgScGeneMod object
#'
#' @export
scgm_from_clusts = function(scl, feat_type = "clust_fp", 
				foc_genes=NULL, min_fold=1,
				use_hclust = F, k_clusts = NA, 
				k_knn = NA, 
				min_tot = 200,
				report_dir = ".",
				cor_method = "pearson",
				rseed = 1)
{
	oldseed = .set_seed(rseed)

	if(use_hclust & is.na(k_clusts)) {
		stop("must define number of clust when using hclust in gmod construction")
	}
	if(!use_hclust & is.na(k_knn)) {
		stop("must define k_knn graph clustering for gmod construction")
	}

	if(feat_type == "clust_fp") {
		scgm = scgm_init_cclust_fp(scl, foc_genes, min_fold)
	} else {
		scgm = scgm_init_cclust_resid(scl, foc_genes, min_tot) 
	}

	scgm = scgm_build_cor_metric(scgm, cor_method)

	if(use_hclust) {
		scgm = .scgm_hclust(scgm, k_clusts, report_dir)
	} else {
		scgm = .scgm_knn_clust(scgm, k_knn, report_dir, rseed)
	}

	scgm = .scgm_postprocess(scgm)

	.restore_seed(oldseed)

	return(scgm)
}

.scgm_hclust = function(scgm, k_clusts, report_dir)
{
#cluster on similarity between correlation with all others
#outliers detection
	g_best_cor = apply(scgm@gene_cor, 1, max)
	T_cor = apply(scgm@rnd_gene_cor, 1, max)+0.1
	good_g = names(which(g_best_cor > T_cor))
#clustering
	scgm@indep_genes= setdiff(rownames(scgm@gene_cor), good_g)
	hc = hclust(dist(scgm@gene_cor[good_g, good_g]), "ward.D2")
	scgm@hc_order = hc$order
	scgm@gmods = cutree(hc, k_clusts)
	return(scgm)
}

.scgm_knn_clust = function(scgm, k_knn, report_dir, rseed)
{
	oldseed = .set_seed(rseed)

#call the man
	min_clust_size = 1

	cat("Running gmod knn clust with K =", k_knn,"...\n")

	genes = rownames(scgm@feat_mat)
	ngenes = length(genes)

	if(ngenes < k_knn*2) {
		stop("K ", k_knn, " for k_knn is too large given only ", ngenes, " feature genes")
	}

	scgm = .scgm_calc_balanced_knn_matrix(scgm, k_knn)

	m_knn = scgm@m_knn
	m_knn = m_knn / k_knn

	clusts = .tg_knn_clust(m_knn, min_clust_size=min_clust_size, consolidate=1)

	clusts1 = as.integer(as.factor(clusts))
	names(clusts1) = names(clusts)
	clusts= clusts1
	message("will reorder")
	clusts = .scgm_reorder_knn_clusts(scgm@feat_mat, clusts)

	scgm@gmods = clusts

	message("will postproc clusts")
	scgm = .scgm_postprocess(scgm)
	.restore_seed(oldseed)

	return(scgm)
}

# calculate mean module expression in cells
.scgm_calc_gmod_mean = function(scgm)
{
	us = as.matrix(scgm@scmat@mat[names(scgm@gmods),])
	
	u_i = colSums(us)
	us_n = t(t(us)/u_i) * mean(u_i)
	
	gmod_mean = .row_stats_by_factor(t(us_n), scgm@gmods, rowSums)

	return(gmod_mean)
}

# create footprints of modules in cells
.scgm_calc_gmod_fp = function(scgm, ref_func = median, fp_reg=0.1)
{
	gmod_fp = t(t(scgm@gmod_mean+fp_reg)/(apply(scgm@gmod_mean, 2, ref_func)+fp_reg))

	return(gmod_fp)
}

.scgm_reorder_knn_gmods = function(feat_mat, gmods) 
{
	n_c = length(unique(clusts))
	# # reorder clusters
	centers = .row_stats_by_factor(t(feat_mat), gmods)
	centers_hc = hclust(dist(cor(centers)), "ward.D2")
	tomap = rep(NA, n_c)
	tomap[centers_hc$order] = 1:n_c
	onames = names(gmods)
	gmods = tomap[gmods]
	names(gmods) = onames
	return(gmods)
}

# return a k_knn X ncells matrix,  column i contains cell i's k nearest neighbors,
# ordered fronm closest to farthest.
# IMPORTANT: calculated matrix might have more rows than expected! always use only the first relevant rows,
# if k is small, will calculate more than requested, since the knn object is used in more than 1 place.
.scgm_calc_knn_ordered = function(scgm, k_knn = -1)
{
	message("Calling calc order knn in gmod.r")

	if(length(scgm@knn_ordered) > 1 && # knn_ordered exists
		 nrow(scgm@knn_ordered) >= k_knn) {
		# no need to calculate again!
		return(scgm)
	}

	ngenes = nrow(scgm@feat_mat)
	minimal_k = sqrt(ngenes)
	if (k_knn == -1)  {
		k_knn = max(30, minimal_k)
	}

	if(k_knn < minimal_k) {
		k_knn = minimal_k
	}

	sim_mat =  scgm@gene_cor

	m_knn_ordered = matrix(NA, nrow = k_knn, ncol = ngenes)

	# ranks = rev(1 - (0:k_knn+1)/(k_knn + 2))
	for(gene in 1:ngenes) {
		knn = head(order(sim_mat[gene,], decreasing = T),n=k_knn+1)
		m_knn_ordered[,gene] = knn[-1] # remove self
		# m_knn[knn, cell] = ranks
	}
	scgm@knn_ordered = m_knn_ordered
	return(scgm)
}

.scgm_calc_balanced_knn_matrix = function(scgm, k_knn, k_expand=2)
{
	ngenes = nrow(scgm@feat_mat)
	genes = rownames(scgm@feat_mat)

	k_knn10 = min(k_knn*10, scgm@scmat@ncells-1)

	scgm = .scgm_calc_ordered_knn(scgm, k_knn10)

	m_knn = matrix(k_knn10, nrow = ngenes,ncol = ngenes, dimnames = list(genes, genes))
	for(gene in 1:ngenes) {
		m_knn[scl@knn_ordered[1:k_knn10, gene], gene] = 1:k_knn10
	}
	M = k_knn*k_knn*k_expand
	m_knn_io = pmax(-m_knn * t(m_knn) + M,0)

	A = nrow(m_knn_io)-k_knn*3	#no more than k_nn*3 contribs per row, or incoming neighbots
	m_knn_io_rows = t(apply(m_knn_io, 1, function(x) pmax(rank(x)-A,0)))
	A = nrow(m_knn_io)-k_knn	#no more than k_nn contribs per column, or outgoing neighbors
	m_knn_io = apply(m_knn_io_rows, 2, function(x) pmax(rank(x)-A,0))

	scgm@m_knn = m_knn_io

	return(scgm)
}

scc_calc_clust_resid_exp = function(scl, foc_genes = NULL, min_g=200)
{
	clusts = scl@clusts
	mat = as.matrix(scl@scmat@mat)[,names(clusts)]
	u_i = colSums(mat)
	mat = t(t(mat)/u_i)

	if(is.null(foc_genes)) {
		foc_genes = names(which(1000*rowSums(mat) > min_g))
	}

	u_gi = mat[foc_genes,]

	isclust_ci = diag(max(scl@clusts))[,scl@clusts]    # 1 encoding of the clusters, rows= clusters, cols = cells

	n_c = as.vector(table(scl@clusts))
	n = sum(n_c)

#we compute total U_{g,b} - number of umis per gene per batch
	e_gc = t(t(u_gi %*% t(isclust_ci))/n_c)
	e_gi = e_gc[,clusts]

	diff_gi = u_gi - e_gi

	return(diff_gi)
}

#' Title
#'
#' @param gm 
#' @param base_dir 
#' @param w 
#' @param h 
#' @param h_fp 
#' @param cex_gene 
#' @param filter_mods 
#' @param colspec 
#' @param mark_mods 
#'
#' @return
#' @export
#'
#' @examples
scgm_plot_matrices = function(gm, base_dir, name="", w=1000, h=1000, h_fp= h, cex_gene=0.6, filter_mods=NA, mark_mods=F, colspec = get_param("divergent_cols"), zlim=NULL, zlim_q=1e-4, res=180)
{

	dir.create(base_dir, showWarnings = F, recursive = T)
	
	.plot_start(sprintf("%s/%sgmod_cor%s.png", base_dir, name, ifelse(is.na(filter_mods), "", sprintf("_%d_mods", length(filter_mods)))), w=w,h=h, res=res)

	good_g = names(gm@gmods)[gm@hc_order]
	if (length(filter_mods) > 1 || !is.na(filter_mods)) {
		good_g = good_g[gm@gmods[good_g] %in% filter_mods]	
	}
	
	shades = colorRampPalette(colspec)(1000)
	
	if (is.null(zlim)) {
		zlim = c(-1, 1) * max(abs(quantile(gm@gene_cor, c(zlim_q, 1-zlim_q))))
	}
	image(pmin(pmax(gm@gene_cor[good_g, good_g], zlim[1]), zlim[2]),
				zlim=zlim,
		xaxt='n', 
		yaxt='n',
		col=shades, 
		main=sprintf("%.2f - %.2f", zlim[1], zlim[2]))

	mtext(good_g, at=seq(0,1,l=length(good_g)), line=1, las=2, side=2, cex=cex_gene)
	mtext(good_g, at=seq(0,1,l=length(good_g)), line=1, las=2, side=1, cex=cex_gene)

	mod_rle = rle(gm@gmods[good_g])
	mlines = (cumsum(mod_rle$lengths)-1) / length(good_g)
	mcent = (cumsum(mod_rle$lengths) - mod_rle$lengths / 2) / length(good_g)
	
	if (mark_mods) {
		abline(h=mlines, lwd=0.5)
		abline(v=mlines, lwd=0.5)
		
		mtext(mod_rle$values, at=mcent, line=1, las=2, side=4, cex=cex_gene*20)
		mtext(mod_rle$values, at=mcent, line=1, las=2, side=3, cex=cex_gene*20)
	} 
	else {
		mtext(good_g, at=seq(0,1,l=length(good_g)), line=1, las=2, side=4, cex=cex_gene)
		mtext(good_g, at=seq(0,1,l=length(good_g)), line=1, las=2, side=3, cex=cex_gene)
		
	}
	dev.off()

	zlim = c(-1, 1) * max(abs(quantile(gm@feat_mat[good_g, ], c(zlim_q, 1-zlim_q))))
	
	.plot_start(sprintf("%s/%sgmod_fp%s.png", base_dir, name, ifelse(is.na(filter_mods), "", sprintf("_%d_mods", length(filter_mods)))), w=w,h=h_fp, res=res)
	
	image(t(pmax(pmin(gm@feat_mat[good_g,], zlim[2]), zlim[1])), zlim=zlim,
					col=shades, yaxt='n', xaxt='n')
	mtext(good_g, at=seq(0,1,l=length(good_g)), line=1, las=2, side=2, cex=cex_gene)
	#mtext(good_g, at=seq(0,1,l=length(good_g)), line=1, las=2, side=4, cex=cex_gene)
	
	abline(h=mlines, lwd=0.5)

	mtext(colnames(gm@feat_mat), at=seq(0, 1, length=ncol(gm@feat_mat)), line=1, las=2, side=1, cex=cex_gene*5)
	mtext(mod_rle$values, at=mcent, line=1, las=2, side=4, cex=cex_gene*5)
	
	dev.off()
	
}

#' Plot total gmod expression per cell on the 2d projection
#'
#' @param sc_2d 
#' @param gm 
#' @param output_dir 
#' @param name 
#' @param gm_labels 
#' @param colspec 
#' @param nshades 
#' @param width 
#' @param height 
#' @param cex 
#' @param lwd 
#' @param zlim 
#'
#' @return
#' @export
#'
#' @examples
scgm_plot_gmod_on_2d = function(sc_2d, gm, output_dir=sprintf("%s/gmod_2d", get_param("outdir")), name="", gm_labels=NULL, colspec = c('white', RColorBrewer::brewer.pal(9,"YlOrRd")), nshades=1000, width=800, height=800, cex=0.5, lwd=0.5, zlim=NULL, zlim_q=1e-4, show_top_n_genes=30) 
{
	dir.create(output_dir, showWarnings = F, recursive=T)
	
	shades = colorRampPalette(colspec)(nshades)

	if (is.null(zlim)) {
		zlim = c(-1, 1) * max(abs(quantile(log2(gm@gmod_mean+1), c(zlim_q, 1-zlim_q))))
	}	

	xlim = range(sc_2d@x)
	if (show_top_n_genes > 0) {
		xlim[1] = xlim[1] - diff(xlim)/5
	}
	
	for (i in 1:gm@nclust) {
		.plot_start(
			sprintf(
				"%s/%sgm%d%s.png",
				output_dir,
				name,
				i,
				ifelse(is.null(gm_labels), "", paste0("_", gm_labels[i]))
			),
			w = width,
			h = height
		)
		
		
		plot(
			sc_2d@x,
			sc_2d@y,
			type = 'p',
			pch = 21,
			col = 'darkgrey',
			bg = shades[.vals_to_n(log2(gm@gmod_mean[, i] + 1), zlim = zlim)],
			lwd = lwd,
			cex = cex,
			xaxt = 'n',
			yaxt = 'n',
			xlab = "",
			ylab = "",
			xlim = xlim,
			main = paste(i, ifelse(is.null(gm_labels), "", gm_labels[i]), zlim[1], zlim[2])
		)

		if (show_top_n_genes > 0) {
			cand_genes = names(gm@gmods)[gm@gmods == i]
			
			if (length(cand_genes) > 1) {
				top_genes = names(head(sort(apply(
					gm@feat_mat[cand_genes,], 1, max), decreasing = T), n = min(length(cand_genes), show_top_n_genes)))
				
				labels = sprintf("%s\t%.2f", top_genes, apply(gm@feat_mat[top_genes, ], 1, max))
			}
			else {
				labels = sprintf("%s\t%.2f", cand_genes, max(gm@feat_mat[cand_genes, ]))
			}
			legend("topleft", legend = labels, bty='n', seg.len = 0)
		}
		
		dev.off()
	}
} 
