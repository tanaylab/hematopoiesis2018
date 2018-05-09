#' scrdb cell sorter
#'
#' given a matrix and a bunch of marker sets, this function split cells by marker expression, and write separate matrices for downstream analysis
#'
#' @param mat 
#' @param markers_fn
#' @param k_scale_umi 
#'
#' @return
#' @export
#'
#' @examples
sc_proj_mat_on_marks = function(mat, markers_fn, 	factor_by_median = T, 
																norm_size = T, k_scale_umi = 7)
{
	markers = read.table(markers_fn, header = T, stringsAsFactors = F)
	
	istype_tg = table(markers$subtype, markers$gene)
	
	good_genes = intersect(rownames(mat@mat), colnames(istype_tg))
	if (length(good_genes) < 3) {
		stop(
			"found only ",
			length(good_genes),
			" marker genes represented in the matrix, terminating"
		)
	}
	
	istype_tg = as.matrix(istype_tg[, good_genes])
	
	umis = as.matrix(mat@mat)[good_genes, ]

	n_median = 1
	if (factor_by_median) {
		n_median = median(colSums(umis))
	}
	if (norm_size) {
		umis = n_median * t(t(umis)/(1+colSums(umis)))
	}
	umis_n = log2(1 + k_scale_umi * umis)
	
	istype_tg = t(t(istype_tg) / colSums(istype_tg))
	
	totmark_it = t(umis_n) %*% t(istype_tg)
	
	totmark_it
}

#' cluster mat projected on markers genes
#'
#' @param mat 
#' @param markers_fn 
#' @param n_mark_clusts 
#' @param analysis_dir 
#' @param alg_type 
#' @param min_clust_size 
#' @param tab_clust_fp_fn 
#' @param clust_fp_metadata_fields 
#'
#' @return
#'
#' @export
sc_marker_clusts = function(mat,
														markers_fn,
														factor_by_median = T,
														norm_size = T,
														k_scale_umi = 7,
														n_mark_clusts = 120,
														analysis_dir = getwd(),
														alg_type = "kmeans",
														min_clust_size=20, 
														tab_clust_fp_fn = "clust_fp.txt",
														clust_fp_metadata_fields = c("Patient", "Sample.Name", "Cell.type"),
														min_max_clust_fp_to_report = get_param("scc_min_max_clust_fp_to_report"),
														ordered_subtypes = c("ab_T_cell", "regulatory_CD4_T_cell", "CD4_T_cells", "CD8_T_cells", "gd_T_cell", "NK_cell", "cytotoxic_B_cells", "IgA_B_cell", "IgM_B_cell", "IgG_B_cell", "plasma_cells", "DC", "macrophage", "granulocytes", "erythrocyte", "osteoclast"), type_col=list("CD3"="#e41a1c", "CD45"="#377eb8", "CD3&CD45"="#4daf4a"), totmark_cl=NULL)
{
	totmark_it = sc_proj_mat_on_marks(mat, markers_fn, factor_by_median = factor_by_median, norm_size = norm_size, k_scale_umi = k_scale_umi)
	
	if (is.null(totmark_cl)) {
		if (alg_type == "knn") {
			totmark_cl = scc_cluster_knn_graph (mat, t(totmark_it), norm_features = F, k_knn = n_mark_clusts, min_clust_size = min_clust_size)
		}
		else if (alg_type == "kmeans") {
			totmark_cl = scc_cluster_kmeans (mat, t(totmark_it), norm_features = F, K=n_mark_clusts)
			#totmark_cl = scc_cluster_kmeans (mat, features=select_markers(mat, gstat = NULL, type = "any",mark.min_var_mean = 0.2, mark.sz_cor_norm_max = -Inf, mark.niche_T = Inf, mark.min_tot = NA), norm_features = T, K=n_mark_clusts)
		}
	}
	
	clusts = totmark_cl@clusts
	
	centers = t(.row_stats_by_factor(t(totmark_it), clusts, rowFunction=rowMeans))
	
	centers_norm = t(t(centers) / colMeans(centers))
	
	if (!is.null(ordered_subtypes)) {
		ordered_subtypes = intersect(ordered_subtypes, colnames(centers_norm))
		subtypes = 1:ncol(centers_norm)
		names(subtypes) = colnames(centers_norm)
		hc_marks = subtypes[ordered_subtypes]
		
		hc_cl = order(apply(centers_norm[, hc_marks], 1, which.max) + 1e-3 * apply(centers_norm[,hc_marks], 1, max))
	}
	else {
		hc_cl = hclust(dist(cor(t(centers_norm))), "ward.D2")$order
		hc_marks = hclust(dist(t(centers_norm)), "ward.D2")$order
	}
	
	
	fp_shades = colorRampPalette(RColorBrewer::brewer.pal(n = 9, "YlOrRd"))(1000)

	cn_m = melt(pmin(centers_norm, 10), value.name='enr')
	
	
	cn_m$cl = factor(cn_m$Var1, levels = hc_cl)
	cn_m$mark = factor(cn_m$Var2, levels = ordered_subtypes)
	
	tot_cl = table(clusts)
	tot_m = melt(tot_cl, value.name='count')
	tot_m$total = sprintf("avg %.1f", mean(tot_cl))
	tot_m$cl = factor(tot_m$clusts, levels=hc_cl)
	

	dist_to_center = sqrt(rowSums((totmark_it - centers[clusts,]) ** 2))
	m_ess = tapply(dist_to_center, clusts, mean)
	ess_m = data.frame(clust=names(m_ess), ess=m_ess, mean_m=sprintf("avg %.2f", mean(dist_to_center)))
	ess_m$cl = factor(ess_m$clust, levels=hc_cl)

		
	type_cl = table(totmark_cl@clusts, totmark_cl@scmat@cell_metadata[names(totmark_cl@clusts), 'Cell.type'])
	type_cl_n = type_cl / rowSums(type_cl)
	
	type_m = melt(type_cl_n, value.name='frac')
	type_m$cl = factor(type_m$Var1, levels=hc_cl)
	
	#png(sprintf("%s/subtype_clusts_%s.png", analysis_dir, alg_type),	w = 1500,	h = 1000)
	
	p_hm <- ggplot(cn_m, aes(x=factor(cl), y=factor(mark))) + geom_tile(aes(fill = enr)) + scale_fill_gradientn(colours=RColorBrewer::brewer.pal(n=9, "YlOrRd")) +  xlab('') + ylab('') + theme(axis.text.x=element_text(size=5, angle=90, vjust=0.5), axis.text.y=element_text(size=8), legend.text=element_text(size=8), legend.title=element_blank(), plot.margin=margin(0, 3, 1, 3))
	
	p_tot_barplot <- ggplot(data=tot_m, aes(x=factor(cl), y=count, fill=total)) + geom_bar(stat="identity") + xlab('') + ylab('') + theme(axis.text.x=element_text(size=5, angle=90, vjust=0.5), axis.text.y=element_text(size=8), legend.text=element_text(size=8), legend.title=element_text(size=8), plot.margin=margin(0, 3, 1, 3)) + scale_fill_discrete(name="#cells")
	
	p_ess_barplot <- ggplot(data=ess_m, aes(x=factor(cl), y=ess, fill=mean_m)) + geom_bar(stat="identity") + xlab('') + ylab('') + theme(axis.text.x=element_text(size=5, angle=90, vjust=0.5), axis.text.y=element_text(size=8), legend.text=element_text(size=8), legend.title=element_text(size=8), plot.margin=margin(0, 3, 1, 3)) + scale_fill_discrete(name="withinss")
	
	p_type_barplot <- ggplot(data=type_m, aes(x=factor(cl), y=frac, fill=Var2)) + geom_bar(stat="identity") + xlab('') + ylab('') + theme(axis.text.x=element_text(size=5, angle=90, vjust=0.5), axis.text.y=element_text(size=8), legend.text=element_text(size=8), legend.title=element_blank(), plot.margin=margin(0, 3, 1, 3))
	
	plot_grid(p_hm, p_tot_barplot, p_ess_barplot, p_type_barplot, nrow=4, ncol=1, align='v', rel_heights=c(4,1,1,1))
	
	ggsave(sprintf("%s/subtype_clusts_%s_byMed%.1s_bySize%.1s.png", analysis_dir, alg_type, factor_by_median, norm_size), width=24, height=16, units="cm")
	
	
	max_center_dist = tapply(dist_to_center, clusts, max)
	
	#message(sprintf("withinss = %.2f", mean(sqrt(rowSums((totmark_it - centers[permute(clusts),]) ** 2)))))
	
	centers_data = cbind(max_center_dist, centers)
	colnames(centers_data)[1] = 'max_center_dist'
	write.table(
		centers_data[hc_cl, ],
		sprintf("%s/subtype_centers_%s.txt", analysis_dir, alg_type),
		quote = F,
		sep = "\t"
	)
	
	if (!is.na(tab_clust_fp_fn)) {
		f = apply(totmark_cl@clust_fp, 1, max) > min_max_clust_fp_to_report
		
		if (length(clust_fp_metadata_fields) > 1 || !is.na(clust_fp_metadata_fields)) {
			for (s in clust_fp_metadata_fields) {
				write.table(table(totmark_cl@scmat@cell_metadata[names(totmark_cl@clusts), s], totmark_cl@clusts), sprintf("%s/%s_%s.%s", analysis_dir, alg_type, tab_clust_fp_fn, s), sep = "\t", quote = F)	
			}
		}
		write.table(
			round(totmark_cl@clust_fp[f, ], 2),
			sprintf("%s/%s_%s", analysis_dir, alg_type, tab_clust_fp_fn),
			sep = "\t",
			quote = F
		)
	}
	
	return(totmark_cl)
}

#' assign cells to types by given projected markers cluster centers
#'
#' @param mat 
#' @param markers_fn gene to type table, project umi table onto types. If null, use the colnames of ref_centers_fn as markers and calc feature matrix
#' @param ref_centers_fn 
#' @param ref_centers_type_assign_fn 
#' @param factor_by_median 
#' @param norm_size 
#' @param k_scale_umi 
#' @param fig_pref 
#' @param out_base_dir 
#'
#' @return
#'
#' @export

sc_marker_split = function(mat,
													 markers_fn,
													 ref_centers_fn,
													 ref_centers_type_assign_fn,
													 factor_by_median = T,
													 norm_size = T,
													 k_scale_umi = 7,
													 fig_pref = NULL,
													 out_base_dir = NULL)
{
	
		
	ref_centers = read.table(
		ref_centers_fn,
		sep = "\t",
		header = T,
		stringsAsFactors = F,
		check.names = F
	)
	ref_centers = ref_centers[order(as.numeric(rownames(ref_centers))), ]
	
	ref_centers_types = read.table(
		ref_centers_type_assign_fn,
		sep = "\t",
		header = T,
		stringsAsFactors = F
	)
	ref_centers_types = ref_centers_types[order(as.numeric(rownames(ref_centers_types))), 'type']
	
	
	cl_max_dist = ref_centers[, 1]
	ref_centers = ref_centers[, -1]

	if (is.null(markers_fn)) {
		totmark_it = t(as.matrix(create_feature_mat(mat, colnames(ref_centers), norm_features = T,
												 factor_by_median = factor_by_median,
												 norm_size = norm_size, k_scale_umi = k_scale_umi)))
	}
	else {
		totmark_it = sc_proj_mat_on_marks(mat, markers_fn, factor_by_median = factor_by_median, norm_size = norm_size, k_scale_umi = k_scale_umi)
	}
	
	good_g = intersect(colnames(totmark_it), colnames(ref_centers))
	totmark_it = totmark_it[, good_g]
	ref_centers = ref_centers[, good_g]
	
	dists = as.matrix(pdist(totmark_it, ref_centers))
	clusts = apply(dists, 1, which.min)
	dist_to_center = sqrt(rowSums((totmark_it - ref_centers[clusts,]) ** 2))
	valid_cells = dist_to_center <= cl_max_dist[clusts]
	
	clusts[!valid_cells] = max(clusts) + 1
	
	ref_centers_types = c(ref_centers_types, 'None')
	
	mat@cell_metadata$assigned_type = ref_centers_types[clusts]
	mat@cell_metadata$vsort_cl = clusts
	
	if (!is.null(out_base_dir)) {
		ind = split(1:nrow(mat@cell_metadata), mat@cell_metadata$assigned_type)
		lapply(names(ind), function(st) {
			message(sprintf("writing %s cells...", st))
			dir.create(sprintf("%s/%s", out_base_dir, st), showWarnings = F)
			write.table(as.matrix(mat@mat[, ind[[st]]]), sprintf("%s/%s/%s_sc_mat.txt", out_base_dir, st, st), quote=F, sep="\t")
			write.table(mat@cell_metadata[ind[[st]],], sprintf("%s/%s/%s_sc_mat_md.txt", out_base_dir, st, st), quote=F, sep="\t")
		})
	}
	
	if (!is.null(fig_pref)) {
		tnames = unique(ref_centers_types[-length(ref_centers_types)])
		tcols = RColorBrewer::brewer.pal(n = length(tnames), 'Set1')
		names(tcols) = tnames
		
		cell_comp = table(c(ref_centers_types[clusts], tnames)) - 1
		
		png(sprintf("%s_cell_type_comp.png", fig_pref),
				w = 800,
				h = 600)
		barplot(cell_comp / sum(cell_comp),
						ylab = "Fraction",
						col = tcols[names(cell_comp)])
		dev.off()
		
		
		png(sprintf("%s_cells_per_ref_cl.png", fig_pref),
				w = 800,
				h = 400)
		tab = table(c(clusts, 1:nrow(ref_centers))) - 1
		
		barplot(tab,
						ylab = "#cells",
						col = tcols[ref_centers_types],
						ylim = c(0, 1.1 * max(tab)))
		legend(
			"topleft",
			legend = names(tcols),
			fill = tcols,
			ncol = length(tcols),
			bty = 'n',
			cex = 0.9
		)
		dev.off()
		
	}
	
	mat
}

#' Plots markers x clusters, supporting addition of additional similiar markers. 
#'
#' @param sc_cl 
#' @param markers_fn 
#' @param type 
#' @param add_genes_by_cor_above 
#' @param add_n_similar_genes 
#' @param k_scale_umi 
#' @param max_fp 
#' @param width 
#' @param height 
#' @param res 
#'
#' @return
#' @export
#'
#' @examples
sc_plot_sort_markers_on_clusts <-
	function(sc_cl,
					 markers_fn,
					 remove_clustering_outliers=T,
					 type = "clust_fp",
					 min_max_clust_fp_to_report = get_param("scc_min_max_clust_fp_to_report"),
					 add_genes_by_cor_above = NA,
					 add_n_similar_genes = 0,
					 k_scale_umi = 7, 
					 max_fp = 2,
					 width = 2500, 
					 height = 2500, 
					 res = 220,
					 cl_ord=NULL)
	{
		if (remove_clustering_outliers) {
			outl = find_clustering_outliers(sc_cl)
			sc_cl = scc_move_cells_to_orphan(sc_cl, outl)
		}
		
		if (type == "clust_fp") {
			feat_mat = log2(sc_cl@clust_fp)
		}
		else if (type == "norm_umi") {
			feat_mat = sc_cl@scmat@mat
			feat_mat = feat_mat[rowSums(feat_mat) > 0,]
			feat_mat_n = t(t(feat_mat) / (1 + colSums(feat_mat)))
			feat_mat = as.matrix(log2(1 + k_scale_umi * feat_mat_n))
		}
		
		markers = read.table(markers_fn, header = T, stringsAsFactors = F)
		m_genes = unique(markers$gene)
		non_m_genes = setdiff(rownames(feat_mat), m_genes)
		
		add_genes = c()
		if (!is.na(add_genes_by_cor_above) | add_n_similar_genes > 0) {
			cc2 = cor(t(feat_mat[m_genes,]), t(feat_mat[non_m_genes,]))
			
			
			if (add_n_similar_genes == 0) {
				add_genes = unique(unlist(apply(cc2, 1, function(v) {
					colnames(cc2)[which(v >= add_genes_by_cor_above)]
				})))
				m_genes = c(m_genes, add_genes)
			}
			else {
				for (i in 1:add_n_similar_genes) {
					curr_add = unique(colnames(cc2)[apply(cc2, 1, which.max)[apply(cc2, 1, max, na.rm =
																																						T) > add_genes_by_cor_above]])
					message(sprintf("iter %d: adding %d genes", i, length(curr_add)))
					m_genes = c(m_genes, curr_add)
					cc2 = cc2[, setdiff(colnames(cc2), curr_add)]
					add_genes = c(add_genes, curr_add)
				}
			}
		}
		
		cc = cor(cor(t(as.matrix(feat_mat[m_genes,]))))
		hc_g = hclust(dist(cc), method = "ward.D2")
		
		tcol = RColorBrewer::brewer.pal(n = length(unique(markers$type)), 'Set1')
		names(tcol) = unique(markers$type)
		
		subt_ann = table(markers$gene, markers$subtype)
		subt_df = data.frame(row.names = rownames(subt_ann))
		
		for (i in 1:ncol(subt_ann)) {
			subt_df[, i] = subt_ann[, i]
		}
		colnames(subt_df) = colnames(subt_ann)
		
		t_subt = unique(markers[, c('type', 'subtype')])
		rownames(t_subt) = t_subt$subtype
		t_subt$color = unlist(tcol[t_subt$type])
		
		subt_df = subt_df[, t_subt[order(t_subt$type), 'subtype']]
		
		subt_contig_col = list()
		for (subt in colnames(subt_df)) {
			subt_contig_col[[subt]] = list('0' = "white", '1' = t_subt[subt, 'color'])
		}
		
		#png(sprintf("%s/marker_cor_hm_add_top%d_or_above%.2f.png", get_param("outdir"), add_n_similar_genes, add_genes_by_cor_above), w=1024, h=960, res=600, pointsize=2)
		png(
			sprintf(
				"%s/marker_cor_hm_%s_add_top%d_or_above%.2f.png",
				get_param("outdir"),
				type,
				add_n_similar_genes,
				add_genes_by_cor_above
			),
			w = width,
			h = height,
			res = res,
			pointsize = 3
		)
		pheatmap(
			cc[hc_g$order, hc_g$order],
			cluster_rows = F,
			cluster_col = F,
			annotation_col = subt_df,
			annotation_colors = subt_contig_col,
			annotation_legend = F,
			fontsize = 3
		)
		dev.off()
		
		png(
			sprintf(
				"%s/marker_lfp_hm_%s_add_top%d_or_above%.2f.png",
				get_param("outdir"),
				type,
				add_n_similar_genes,
				add_genes_by_cor_above
			),
			w = width,
			h = height,
			res = res,
			pointsize = 4
		)
		lfp = as.matrix(feat_mat[m_genes, ])
		
		rhc = hclust(dist(cor(t(lfp))), method="ward.D2")
		if (is.null(cl_ord)) {
			chc = hclust(dist(cor(lfp)), method="ward.D2")
			cl_ord = chc$order
		}
		pheatmap(
			pmin(pmax(lfp[rhc$order, cl_ord], -max_fp), max_fp), 
			cluster_rows = F,
			cluster_col = F,
			annotation_row = subt_df,
			annotation_colors = subt_contig_col,
			annotation_legend = F,
			fontsize = 2,
			fontsize_row = 4
		)
		dev.off()
		
		# plot clusters stats
		png(sprintf(
			"%s/marker_cl_stat_%s_add_top%d_or_above%.2f.png",
			get_param("outdir"),
			type,
			add_n_similar_genes,
			add_genes_by_cor_above
		),
		w = width,
		h = width/2,
		res = res,
		)
		
		bd_field = get_param("scm_cells_breakdown_field")
		
		if (!is.null(bd_field)) {
			layout(matrix(1:2, 2, 1))
		}
		par(mar=c(1,4,1,1))
		mcl = tapply(colSums(as.matrix(sc_cl@scmat@mat)), sc_cl@clusts, mean)
		barplot(mcl[cl_ord], col='black', xaxt='n', ylab='mean umi')
		
		if (!is.null(bd_field)) {
			tcl = table(sc_cl@clusts, sc_cl@scmat@cell_metadata[, bd_field])
			tcl = tcl[cl_ord, ]
			tcl_n = tcl / rowSums(tcl)
			
			tcl_cols = 'black'
			
			tcols = get_param("scm_cell_col_dict")
			
			if (!is.null(tcols)) {
				tcl_cols = unlist(tcols[colnames(tcl_n)])
			}
			barplot(t(tcl_n), col=tcl_cols, border=tcl_cols, xaxt='n', ylim=c(-0.2, 1), ylab="fraction")
			if (!is.null(tcols)) {
				legend("bottomright", legend=colnames(tcl), fill=unlist(tcols[colnames(tcl)]), ncol=ncol(tcl), bty='n')
			}
		}
		
		dev.off()								
		
		# report mean cor of newly found genes to old markers (stratify by marker type)
		g2type = unique(markers[, c('gene', 'type')])
		rownames(g2type) = g2type$gene
		cc_new = cor(t(as.matrix(feat_mat[add_genes,])), t(as.matrix(feat_mat[unique(markers$gene),])))
		new_type_sim = apply(cc_new, 1, function(v) { tapply(v, g2type[colnames(cc_new), 'type'], mean) })
		write.table(t(new_type_sim), sprintf("%s/marker_cor_%s_add_top%d_or_above%.2f_found_marks.txt", get_param("outdir"), type, add_n_similar_genes, add_genes_by_cor_above), quote=F, sep="\t")
		
		# report features table clusters' centers 
		clusts = sc_cl@clusts
		fmat = sc_cl@feat_mat
		
		centers = .row_stats_by_factor(fmat, clusts, rowFunction=rowMeans)
		dist_to_center = sqrt(colSums((fmat - centers[,clusts]) ** 2))
		max_center_dist = tapply(dist_to_center, clusts, max)
		
		centers_data = cbind(max_center_dist, t(centers))
		colnames(centers_data)[1] = 'max_center_dist'
		write.table(
			centers_data[cl_ord, ],
			sprintf("%s/sort_markers_%s_add_top%d_or_above%.2f_centers.txt", get_param("outdir"), type, add_n_similar_genes, add_genes_by_cor_above),
			quote = F,
			sep = "\t"
		)
		
		
		# output ordered clust_fp table
		f = apply(sc_cl@clust_fp, 1, max) > min_max_clust_fp_to_report
		write.table(round(sc_cl@clust_fp[f, cl_ord], 2), sprintf("%s/sort_markers_%s_add_top%d_or_above%.2f_clust_fp.txt", get_param("outdir"), type, add_n_similar_genes, add_genes_by_cor_above), quote=F, sep="\t")
		
		lfp[rhc$order, cl_ord]
}

