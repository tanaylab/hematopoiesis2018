#' scrdb pipeline
#'
#' pipeline for reading mars-seq, cleaning, finding cell and gene modules, 
#' and plotting results.
#'

#' sc_pipe_build_mat
#'
#' reading mars-seq into a single umi matrix
#'
#' @param scdb initialized scdb object
#' @param index_fn 
#' @param base_dir 
#' @param batch_meta_attr 
#' @param mars_batches 
#' @param sample_n_batches 
#' @param outdir 
#' @param outdir_add_timestamp 
#' @param min_umi_n 
#' @param min_umis_init 
#'
#' @return returns invisibly a tgScMat
#'
#' @export
sc_pipe_build_mat = function(scdb,
														 index_fn = get_param("scm_batch_table_fn"),
														 base_dir = get_param("scm_mars_base_dir"),
														 batch_meta_attr = get_param("scm_batch_meta_attr"),
														 mars_batches = get_param("scm_mars_batches"),
														 sample_n_batches = get_param("scm_n_batches_to_sample"),
														 outdir = get_param("outdir"),
														 outdir_add_timestamp = F,
														 min_umi_n = get_param("scm_min_cell_umis"),
														 min_umis_init = min(get_param("scm_min_cell_umis_init"), min_umi_n),
														 dataset_type = get_param("scm_dataset_type"),
														 tenex_mat_ifn = get_param("scm_10x_mat_ifn"),
														 tenex_batches = get_param("scm_10x_batches"))
{
	sc_mat = scdb_load(scdb, "raw_mat")
	
	if (is.null(sc_mat)) {
		if (outdir_add_timestamp) {
			outdir = paste0(outdir, "/", format(Sys.time(), "%Y.%m.%d_%H.%M"))
		}
		
		if (!dir.exists(outdir)) {
			dir.create(outdir)
		}
		
		# read batch index
		index = read.delim(index_fn, stringsAsFactors = F, header = T)
		
		# use selected batches
		use_batches = mars_batches

		if (dataset_type == "10x") {
			use_batches = tenex_batches
		}
		if (!is.null(use_batches)) {
			bat = intersect(use_batches, index[, batch_meta_attr])
			index = index[index[, batch_meta_attr] %in% bat, ]
		}
		
		# sample batches
		if (!is.null(sample_n_batches)) {
			sample_n_batches = min(sample_n_batches, nrow(index))
			index = index[sample(1:nrow(index), sample_n_batches, replace = F), ]
		}
		
		# read all mars UMIs
		if (dataset_type == "mars") {
			sc_mat = scm_read_scmat_mars(
				base_dir = base_dir,
				mars_batches = index[, batch_meta_attr],
				batch_meta = index,
				min_umis_n = min_umis_init
			)
		}
		else if (dataset_type == "10x") {
			sc_mat = scm_read_scmat_10x(
				sparse_matrix_ifn = tenex_mat_ifn,
				use_batches = index[, batch_meta_attr],
				min_umis_n = min_umis_init
			)
		}
		
		# add FACS indices (if those exist, see scm_fac_idx_... params)
		if (!is.null(get_param("scm_facs_idx_base_dir"))) {
			sc_mat = scm_add_facs_indices_from_tables(sc_mat)
		}
		
		scdb_save(scdb, sc_mat, "raw_mat")
	}
	
	invisible(sc_mat)
}


#' sc_pipe_mat_to_clean_mat
#'
#' cleaning the umi table
#'
#' @param outdir 
#' @param outdir_add_timestamp 
#' @param clust_knn 
#' @param min_umi_n 
#' @param min_umis_init 
#' @param mark.min_var_mean 
#' @param mark.sz_cor_norm_max 
#' @param mark.niche_T 
#' @param mark.min_tot 
#' @param glob_blacklist_terms 
#' @param mark_blacklist_terms regular expression terms to search blacklisted genes
#' @param mark_blacklist_genes exact names of blacklist genes 
#' @param outlier_gene_top1 
#' @param outlier_gene_on_count 
#' @param remove_outliers_before_clustering 
#' @param min_clust_size 
#' @param filt_amb_on_clusts 
#' @param filt_outliers_on_clusts 
#' @param scdb 
#' @param raw_sc_mat 
#' @param batch_meta_attr 
#' @param amb_epsilon 
#'
#' @return returns invisibly a tgScMat
#'
#' @export
sc_pipe_mat_to_clean_mat = function(scdb, raw_sc_mat,
												 outdir = get_param("outdir"),
												 outdir_add_timestamp = F,
												 batch_meta_attr = get_param("scm_batch_meta_attr"),
												 clust_knn = get_param("scc_k_knn"),
												 min_umi_n = get_param("scm_min_cell_umis"),
												 min_umis_init = min(get_param("scm_min_cell_umis_init"), min_umi_n),
												 mark.min_var_mean = get_param("scc_marker_min_var_mean"),
												 mark.sz_cor_norm_max = get_param("scc_marker_sz_cor_norm_max"),
												 mark.niche_T = get_param("scc_marker_niche_T"),
												 mark.min_tot = get_param("scc_marker_min_tot"),
												 glob_blacklist_terms = get_param("scm_gene_blacklist_names"),
												 mark_blacklist_terms = get_param("scc_marker_blacklist_names"),
												 mark_blacklist_genes = get_param("scc_marker_blacklist_genes"),
												 outlier_gene_top1 = get_param("scm_outlier_max_gene_top1"),
												 outlier_gene_on_count = get_param("scm_outlier_gene_on_count"),
												 remove_outliers_before_clustering = get_param("scm_remove_outliers"),
												 min_clust_size = get_param("scc_min_clust_size"),
												 filt_amb_on_clusts = get_param("scm_filter_amb_noise_on_clusts"),
												 filt_outliers_on_clusts = get_param("scm_filter_outliers_on_clusts"),
												 amb_epsilon = get_param("scm_amb_noise_epsilon"),
												 filt_lateral_supervised_genes = get_param("scm_filt_lateral_supervised_genes"),
												 filt_lateral_supervised_thresh = get_param("scm_filt_lateral_supervised_thresh"),
												 filt_focus_supervised_genes = get_param("scm_filt_focus_supervised_genes")
												 )
{

	sc_mat = scdb_load(scdb, "clean_mat")
	
	if (is.null(sc_mat)) {
		
		if (is.null(raw_sc_mat)) {
			raw_sc_mat = scdb_load(scdb, "raw_mat")
			if (is.null(raw_sc_mat)) {
				stop("raw umi mat not supplied and not found in scdb object")
			}
		}
		
		sc_mat = raw_sc_mat
		
		if (outdir_add_timestamp) {
			outdir = paste0(outdir, "/", format(Sys.time(), "%Y.%m.%d_%H.%M"))
		}
		
		if (!dir.exists(outdir)) {
			dir.create(outdir)
		}
		
		if (length(glob_blacklist_terms) > 0) {
			# blacklist for markers
			glob_blist = c()
			allg = rownames(sc_mat@mat)
			for (bad in glob_blacklist_terms) {
				glob_blist = union(glob_blist, grep(bad, allg, v = T))
			}
			if (length(glob_blist) > 0) {
				sc_mat = scm_sub_mat(sc_mat,
														 genes = setdiff(rownames(sc_mat@mat), glob_blist))
			}
		}
		
		# remove cross well contamination
		sc_mat = scm_remove_ambient_by_epsilon(sc_mat, amb_epsilon, batch_meta_attr, min_umi_n = min_umi_n)
		# calc gstat
		gstat = umi_gene_stat(sc_mat)
		
		# remove_outliers (extract to scmat?)
		if (remove_outliers_before_clustering) {
			outlier_genes = rownames(gstat[gstat$ds_top1 > outlier_gene_top1
																		 &
																		 	gstat$is_on_count < outlier_gene_on_count, ])
			if (length(outlier_genes) > 1) {
				outlier_cells = names(which(colSums(as.matrix(sc_mat@mat)[outlier_genes, ]) >
																			outlier_gene_top1))
				message(
					"removing ",
					length(outlier_cells),
					" outlier cells and ",
					length(outlier_genes),
					" outlier genes"
				)
				if (length(outlier_cells) > 0) {
					sc_mat = scm_sub_mat(
						sc_mat,
						cells = setdiff(colnames(sc_mat@mat), outlier_cells),
						genes = setdiff(rownames(sc_mat@mat), outlier_genes)
					)
				}
				gstat = gstat[setdiff(rownames(gstat), outlier_genes), ]
			}
		}
		
		if (filt_amb_on_clusts | filt_outliers_on_clusts) {
			# find markers
			marks = select_markers(
				sc_mat,
				gstat = gstat,
				type = "any",
				mark.min_var_mean = mark.min_var_mean,
				mark.sz_cor_norm_max = mark.sz_cor_norm_max,
				mark.niche_T = mark.niche_T,
				mark.min_tot = mark.min_tot
			)
			
			# blacklist for markers
			blist = c()
			allg = rownames(gstat)
			for (bad in mark_blacklist_terms) {
				blist = union(blist, grep(bad, allg, v = T))
			}
			blist = union(blist, mark_blacklist_genes)
			
			marks = setdiff(marks, blist)
			if (length(marks) < 10) {
				stop(
					"found ",
					length(marks),
					" markers only, consider changing config to allow more genes qualifying as informative"
				)
			}
			
			plot_markers(
				mat = sc_mat@mat,
				marks =  marks,
				gstat = gstat,
				fname = "markers.png",
				outdir = outdir
			)
		
			# cluster as an additional tool for cleaning
			sc_cl = scc_cluster_knn_graph(sc_mat,
																		sc_mat@mat[marks, ],
																		k_knn = clust_knn,
																		min_clust_size = min_clust_size)
			
			if (filt_amb_on_clusts) {
				rpt_base_dir <- sprintf("%s/ambient", outdir)
				if (!dir.exists(rpt_base_dir)) {
					dir.create(rpt_base_dir)
				}
				deamb_mat = .scc_clean_ambient_on_clusts(sc_cl,
																								 raw_sc_mat@mat[sc_mat@genes, ],
																								 amb_epsilon,
																								 batch_meta_attr,
																								 rpt_base_dir = rpt_base_dir)
				sc_mat = tgScMat(deamb_mat, "umi", cell_metadata = sc_mat@cell_metadata[colnames(deamb_mat), ])
			}
			if(!is.null(filt_lateral_supervised_genes)) {
				message("attempt to filter lateral genes")
				genes_lateral = scc_get_fpcor_genes(sc_cl, filt_lateral_supervised_genes, filt_lateral_supervised_thresh, filt_focus_supervised_genes)
				fn = sprintf("%s/%s_%s.txt", scdb@basedir, scdb@fn_prefix, "lateral_genes")
				write.table(genes_lateral, fn)
			}
			
			if (filt_outliers_on_clusts) {
				cell_outliers_nms = scc_find_clustering_outliers(sc_cl)
				message("remove ", length(cell_outliers_nms), " outlier cells")
				sc_mat = scm_sub_mat(sc_mat, cells = setdiff(colnames(sc_mat@mat), cell_outliers_nms))
				reclust = T
			}
			
		}
		
		scdb_save(scdb, sc_mat, "clean_mat")
		
	}
	
	invisible(sc_mat)
}


#' sc_pipe_mat_to_cmods
#'
#' cluster input sc_mat
#'
#' @param outdir 
#' @param outdir_add_timestamp 
#' @param clust_knn 
#' @param mark.min_var_mean 
#' @param mark.sz_cor_norm_max 
#' @param mark.niche_T 
#' @param mark.min_tot 
#' @param mark_blacklist_terms
#' @param mark_blacklist_genes 
#' @param filt_outliers_on_clusts 
#' @param min_clust_size 
#' @param tab_clust_fp_fn 
#' @param clust_fp_metadata_fields 
#' @param scdb 
#' @param sc_mat 
#' @param use_bootstapping 
#' @param k_cut 
#' @param N_boot 
#' @param boot_ratio 
#' @param n_proces 
#' @param boot_in_clust_size 
#'
#' @return invisible tgScMatClust clustering object
#' @export
#'
#' @examples
sc_pipe_mat_to_cmods = function(scdb,
																sc_mat,
																outdir = get_param("outdir"),
																outdir_add_timestamp = F,
																use_bootstapping = get_param("scc_bootstrap"),
																k_cut = get_param("scc_bootstrap_k_cut"),
																N_boot = get_param("scc_bootstrap_n_boot"),
																boot_ratio = get_param("scc_bootstrap_ratio"),
																n_procs = get_param("scc_bootstrap_n_procs"),
																boot_min_clust_size = get_param("scc_bootstrap_min_clust_size"),
																clust_knn = get_param("scc_k_knn"),
																mark.min_var_mean = get_param("scc_marker_min_var_mean"),
																mark.sz_cor_norm_max = get_param("scc_marker_sz_cor_norm_max"),
																mark.niche_T = get_param("scc_marker_niche_T"),
																mark.min_tot = get_param("scc_marker_min_tot"),
																mark_blacklist_terms = get_param("scc_marker_blacklist_names"),
																mark_blacklist_genes = get_param("scc_marker_blacklist_genes"),
																markers_for_blacklist_cells = get_param("scc_marker_for_blacklist_cells"),
																thresh_for_blacklist_cells = get_param("scc_thresh_for_blacklist_cells"),
																filt_outliers_on_clusts = get_param("scm_filter_outliers_on_clusts"),
																min_clust_size = get_param("scc_min_clust_size"),
																
																tab_clust_fp_fn = get_param("scc_clust_fp_fn"),
																min_max_clust_fp_to_report = get_param("scc_min_max_clust_fp_to_report"),
																clust_fp_metadata_fields = get_param("scc_clust_fp_metadata_fields"))
{
	sc_cl = scdb_load(scdb, "cell_modules")
	
	if (is.null(sc_cl)) {
		if (is.null(sc_mat)) {
			sc_mat = scdb_load(scdb, "clean_mat")
			if (is.null(sc_mat)) {
				stop("clean umi mat not supplied and not found in scdb object")
			}
		}
		
		if (outdir_add_timestamp) {
			outdir = paste0(outdir, "/", format(Sys.time(), "%Y.%m.%d_%H.%M"))
		}
		
		if (!dir.exists(outdir)) {
			dir.create(outdir)
		}

		gstat = umi_gene_stat(sc_mat)
		
		# markers
		marks = select_markers(
			sc_mat,
			gstat = gstat,
			type = "any",
			mark.min_var_mean = mark.min_var_mean,
			mark.sz_cor_norm_max = mark.sz_cor_norm_max,
			mark.niche_T = mark.niche_T,
			mark.min_tot = mark.min_tot
		)
		
		# blacklist for markers
		blist = c()
		allg = rownames(gstat)
		for (bad in mark_blacklist_terms) {
			blist = union(blist, grep(bad, allg, v = T))
		}
		blist = union(blist, mark_blacklist_genes)

#		lat_fn = sprintf("%s/%s_%s.txt", scdb@basedir, scdb@fn_prefix, "lateral_genes")
#		lateral_blist = read.table(lat_fn)
#		blist = union(blist, lateral_blist$x)
		
		message(sprintf("Out of %d genes to blacklist as markers, %d intersect with candidate markers and removed", length(blist), length(intersect(blist, marks))))
		marks = setdiff(marks, blist)
		if (length(marks) < 10) {
			stop(
				"found ",
				length(marks),
				" markers only, consider changing config to allow more genes qualifying as informative"
			)
		}
		plot_markers(
			mat = sc_mat@mat,
			marks =  marks,
			gstat = gstat,
			fname = "markers.png",
			outdir = outdir
		)
		
		# cluster
		sc_cl = scc_cluster_knn_graph(sc_mat,
																	sc_mat@mat[marks,],
																	k_knn = clust_knn,
																	min_clust_size = min_clust_size)
		
		reclust = F
		cell_outliers_nms = c()
		
		if (filt_outliers_on_clusts) {
			cell_outliers_nms = scc_find_clustering_outliers(sc_cl)
			message("remove ", length(cell_outliers_nms), " outlier cells")
			if(!is.null(markers_for_blacklist_cells)) {
				thresh = thresh_for_blacklist_cells
				markers_for_blacklist_cells = intersect(markers_for_blacklist_cells, rownames(sc_cl@scmat@mat))
				if(length(markers_for_blacklist_cells)>1) {
					cell_blist = names(which(apply(sc_cl@scmat@mat[markers_for_blacklist_cells,], 2, max) > thresh))
					message("remove ", length(cell_blist), "cells due to markers blisted")
				} else {
					cell_blist=  names(which(sc_cl@scmat@mat[markers_for_blacklist_cells,]>thresh))
				}
				cell_outliers_nms = union(cell_outliers_nms, cell_blist)
			}
			sc_mat = scm_sub_mat(sc_mat, cells = setdiff(colnames(sc_mat@mat), cell_outliers_nms))
			reclust = T
		}
		sc_cl_1st = sc_cl
		if (reclust) {
			sc_cl = scc_cluster_knn_graph(sc_mat,
																		sc_mat@mat[marks,],
																		k_knn = clust_knn,
																		min_clust_size = min_clust_size)
		
		}
		if (use_bootstapping) {
			message("starting bootstrap clustering...")
			sc_cl = scc_bootstrap_clusts(sc_cl, scdb, k_cut, N_boot, boot_ratio, n_procs, boot_min_clust_size)
		}	
		
		if (!is.na(tab_clust_fp_fn)) {
			scc_write_clust_fp_tab(sc_cl, 
														 outdir = outdir,
														 tab_clust_fp_fn = tab_clust_fp_fn,
														 min_max_clust_fp_to_report = min_max_clust_fp_to_report,
														 clust_fp_metadata_fields = clust_fp_metadata_fields)
		}
		
		scdb_save(scdb, sc_cl, "cell_modules")
		
	}
	
	invisible(sc_cl)
}

#' Find correletaed gene modules on cell modules
#'
#' @param scdb 
#' @param sc_cl 
#' @param feat_type either clust_fp (correlate genes in clust_fp matrix) or residual (correlate genes in residual matrix after cleaning cluster specific expression)
#' @param outdir 
#' @param foc_genes list of genes to use (null by default, and then finding genes by parameters below)
#' @param min_abs_log2_fold if feat_type is clust_fp, select genes with abs log2 fold change above this value
#' @param min_tot if feat_type is residual, select genes with total umi coutn above this value
#' @param rseed 
#'
#' @return invisible tgScMatClust clustering object
#' @export
#'
#' @examples
sc_pipe_cmods_to_gmods = function(scdb,
																  sc_cl,
																	feat_type = "clust_fp",
																	outdir = get_param("outdir"), 
																	n_mods = get_param("scgm_n_gmods"),
																	foc_genes = get_param("scgm_foc_genes"),
																	min_abs_log2_fold = get_param("scgm_min_abs_log2_fold"),
																	min_tot = get_param("scgm_min_tot_umis"),
																	rseed = 1) 
{
	obj_name = paste0("gene_modules_by_", feat_type)
	sc_gmod = scdb_load(scdb, obj_name)
	
	if (is.null(sc_gmod)) {
		if (is.null(sc_cl)) {
			sc_cl = scdb_load(scdb, "cell_modules")
			if (is.null(sc_cl)) {
				stop("cell modules not supplied and not found in scdb object")
			}
		}
		
		sc_gmod = scgm_from_clusts(
			sc_cl,
			feat_type,
			foc_genes,
			min_fold = min_abs_log2_fold,
			use_hclust = T,
			k_clusts = n_mods, 
			min_tot = min_tot,
			report_dir = outdir,
			rseed = rseed
		)
		
		write.table(data.frame(row.names=names(sc_gmod@gmods), gmod=sc_gmod@gmods), sprintf("%s/%s_gmods.txt", outdir, obj_name), quote=F, sep="\t")
		scdb_save(scdb, sc_gmod, obj_name)
		
	}
	
	invisible(sc_gmod)

}

#' Generate new cell modules by cut tree on a cells' co-clust matrix (does not save the new cmods object in scdb)
#'
#' @param scdb 
#' @param sc_cl 
#' @param k_cut 
#' @param min_clust_size 
#'
#' @return new tgScMatClust object
#' @export
#'
#' @examples
sc_pipe_coclust_to_cell_modules = function(scdb, sc_cl, 
																					 k_cut = get_param("scc_bootstrap_k_cut"), 
																					 min_clust_size = get_param("scc_bootstrap_min_clust_size"),
																					 method=get_param("scc_bootstrap_clust_method"),
																					 outdir = get_param("outdir"),
																					 tab_clust_fp_fn = get_param("scc_clust_fp_fn"),
																					 min_max_clust_fp_to_report = get_param("scc_min_max_clust_fp_to_report"),
																					 clust_fp_metadata_fields = get_param("scc_clust_fp_metadata_fields"))
{
	if (is.null(sc_cl)) {
		sc_cl = scdb_load(scdb, "cell_modules")
		if (is.null(sc_cl)) {
			stop("cell modules not supplied and not found in scdb object")
		}
	}
	scl = scc_coclust_to_cell_modules(scdb, sc_cl, cc=NULL, k_cut, min_clust_size, method=method)
	
	scl@alg_params[['bootstrap_cut']] = k_cut
	scl@alg_params[['min_clust_size']] = min_clust_size

	if (!is.na(tab_clust_fp_fn)) {
		scc_write_clust_fp_tab(scl, 
													 outdir = outdir,
													 tab_clust_fp_fn = tab_clust_fp_fn,
													 min_max_clust_fp_to_report = min_max_clust_fp_to_report,
													 clust_fp_metadata_fields = clust_fp_metadata_fields)
	}	
	
	return (scl) 
}

#' sc_pipe_plots
#'
#' @param scdb. initialized scdb object
#' @param sc_cl. cell modules object to use for plotting. Ignored if plot object is available from scdb, search it in scdb if sc_cl = NULL
#' @param outdir 
#' @param mark_blacklist_terms 
#' @param mark_blacklist_genes
#' @param mark_colors_fn 
#' @param fig_width_confu 
#' @param fig_mat_per_clust_genes 
#' @param fig_mat_gene_min_fold 
#' @param fig_mat_gene_min_cov 
#' @param fig_mat_width 
#' @param fig_mat_height 
#' @param fig_mat_text_cex 
#' @param fig_mat_smooth_n 
#' @param fig_2d_height 
#' @param fig_2d_width 
#' @param K_2dproj 
#' @param K_cells_2dproj 
#' @param force_max_deg_2dproj 
#' @param T_edge_2dproj 
#' @param T_edge_asym 
#' @param restrict_edges_by_fpcor 
#' @param expand_K_connect_2dproj 
#' @param fold_shades 
#' @param focus_tfs_fn 
#' @param plot_genes_2d 
#' @param meta_field_nm 
#' @param highlight_cells_by_meta_field_nms 
#' @param show_2d_color_markers 
#' @param mat_focus_genes_fn 
#'
#' @return returns invisibly a tgScMatClust2D plotting object.
#'
#' @export
sc_pipe_plots = function(scdb,
												 sc_cl,
												 outdir = get_param("outdir"),
												 mark_blacklist_terms = get_param("scc_marker_blacklist_names"),
												 mark_blacklist_genes = get_param("scc_marker_blacklist_genes"),
												 mark_colors_fn = get_param("scp_mark_colors_fn"),
												 fig_width_confu = get_param("scp_width_confu"),
												 fig_mat_per_clust_genes = get_param("scp_ord_per_clust_genes"),
												 fig_mat_gene_min_fold = get_param("scp_ord_gene_min_fold"),
												 fig_mat_gene_min_cov = get_param("scp_ord_gene_min_cov"),
												 fig_mat_width = get_param("scp_mat_width"),
												 fig_mat_height = get_param("scp_mat_height"),
												 fig_mat_text_cex = get_param("scp_mat_text_cex"),
												 fig_mat_smooth_n = get_param("scp_mat_smooth_n"),
												 fig_2d_height = get_param("scp_2d_height"),
												 fig_2d_width = get_param("scp_2d_width"),
												 K_2dproj = get_param("scp_K_2dproj"),
												 K_cells_2dproj = get_param("scp_K_cells_2dproj"),
												 force_max_deg_2dproj = get_param("scp_force_max_deg_2dproj"),
												 T_edge_2dproj = get_param("scp_T_edge_2dproj"),
												 T_edge_asym = get_param("scp_T_edge_asym"),
												 restrict_edges_by_fpcor = get_param("scp_restrict_edges_by_fpcor"),
												 expand_K_connect_2dproj = get_param("scp_expand_K_connect_2dproj"),
												 fold_shades = colorRampPalette(get_param("divergent_cols"))(1000),
												 focus_tfs_fn = get_param("scp_focus_tfs_fn"),
												 plot_genes_2d = get_param("scp_plot_genes_2d"),
												 meta_field_nm = get_param("scp_metadata_fig_factor"),
												 highlight_cells_by_meta_field_nms = get_param("scp_highlight_cells_by_meta_field_nms"),
												 mat_focus_genes_fn = get_param("scp_mat_focus_genes_fn"),
												 show_2d_color_markers = get_param("scp_show_2d_color_markers")
) {
	sc_2d = scdb_load(scdb, "plot")
	
	if (is.null(sc_2d)) {
		if (is.null(sc_cl)) {
			sc_cl = scdb_load(scdb, "cell_modules")
			if (is.null(sc_cl)) {
				stop("cell modules not supplied and not found in scdb object")
			}
		}
		
		#  plot  #
		
		sc_2d = scp_compute_clust_knn_graph(
			scp_init_plot(sc_cl),
			force_max_deg = force_max_deg_2dproj,
			K = K_2dproj,
			T_edge = T_edge_2dproj,
			T_edge_asym = T_edge_asym,
			subselect_by_fpcor = restrict_edges_by_fpcor,
			k_expand_inout_factor = expand_K_connect_2dproj
		)
	
		blist = c()
		allg = sc_cl@scmat@genes
		for (bad in mark_blacklist_terms) {
			blist = union(blist, grep(bad, allg, v = T))
		}
		sc_2d@blist_genes = union(blist, mark_blacklist_genes)
		
		# assign colors to clusters and cells
		marker_colors = NULL
		if (!is.null(mark_colors_fn)) {
			mcol = read.table(
				mark_colors_fn,
				sep = "\t",
				h = T,
				stringsAsFactors = F
			)
			marker_colors = data.frame(color = mcol$color)
			rownames(marker_colors) = mcol$gene
			marker_colors$group = mcol$group
			if ("priority" %in% colnames(mcol)) {
				marker_colors$priority = mcol$priority
			}
			else {
				marker_colors$priority = 1
			}
		}
		sc_2d = scp_set_clust_cols(sc_2d, marker_colors = marker_colors)
		
		mat_genes = list(marks = NULL)
		sc_2d = scp_plot_confusion_mat(sc_2d,
																	 outdir = outdir,
																	 height = fig_width_confu,
																	 width = fig_width_confu)
		
		if (!is.null(mat_focus_genes_fn)) {
			ct_marks = read.table(mat_focus_genes_fn,
														header = T,
														stringsAsFactors = F)
			mat_genes[['sort_marks']] = unique(ct_marks$gene)
		}

		for (s in names(mat_genes)) {
			sc_2d = scp_plot_mat(
				sc_2d,
				png_fn = sprintf("%s/%s_all_mat.png", outdir, s),
				per_clust_genes = fig_mat_per_clust_genes,
				gene_min_fold = fig_mat_gene_min_fold,
				gene_min_cov = fig_mat_gene_min_cov,
				width = fig_mat_width,
				height = fig_mat_height,
				text_cex = fig_mat_text_cex,
				smooth_n = fig_mat_smooth_n,
				blacklist = sc_2d@blist_genes,
				remove_blacklist = F,
				genes = mat_genes[[s]]
			)
			sc_2d = scp_plot_mat(
				sc_2d,
				plot_cells = F,
				png_fn = sprintf("%s/%s_clust_mat.png", outdir, s),
				fp_shades = fold_shades,
				per_clust_genes = fig_mat_per_clust_genes,
				gene_min_fold = fig_mat_gene_min_fold,
				gene_min_cov = fig_mat_gene_min_cov,
				width = fig_mat_width,
				height = fig_mat_height,
				text_cex = fig_mat_text_cex,
				blacklist = sc_2d@blist_genes,
				remove_blacklist = F,
				genes = mat_genes[[s]]
			)
		}
		
		# plot metadata
		if (is.null(meta_field_nm)) {
			meta_field_nm = colnames(sc_2d@scl@scmat@cell_metadata)[1]
		}
		sc_2d = scp_plot_metadata_factor(
			sc_2d,
			clust_breakdown = T,
			heatmap = T,
			meta_field_nm = meta_field_nm,
			width = fig_mat_width,
			height = fig_mat_width / 2,
			fname = paste0(outdir, "/metadata.png")
		)
		
		# projections
		sc_2d = scp_plot_clust_2d(
			sc_2d,
			fname = paste0(outdir, "/all_2d.png"),
			height = fig_2d_height,
			width = fig_2d_width,
			plot_markers = show_2d_color_markers,
			K_for_cells = K_cells_2dproj
		)
		sc_2d = scp_plot_clust_2d(
			sc_2d,
			fname = paste0(outdir, "/graph_2d.png"),
			height = fig_2d_height,
			width = fig_2d_width,
			plot_edges = T,
			plot_markers = show_2d_color_markers,
			K_for_cells = K_cells_2dproj
		)
		
		write.table(
			cbind(
				data.frame(
					id = names(sc_2d@scl@clusts),
					x = sc_2d@x,
					y = sc_2d@y,
					clust = sc_2d@scl@clusts,
					color = sc_2d@c_colors
				),
				sc_2d@scl@scmat@cell_metadata
			),
			paste0(outdir, "/cells_xy_metadata.txt"),
			sep = "\t",
			quote = F
		)
		
		# color cells by metadata fields
		if (!is.null(highlight_cells_by_meta_field_nms)) {
			for (meta_field in highlight_cells_by_meta_field_nms) {
				message(sprintf("2d projection by %s", meta_field))
				for (single_plot in c(F, T)) {
					sc_2d = scp_plot_clust_2d_by_meta_field(
						sc_2d,
						meta_field,
						odir = sprintf("%s/all_2d_by_%s", outdir, meta_field),
						K_for_cells = K_cells_2dproj,
						bg_col = "grey90",
						fg_col = NULL,
						panel_size = 250 * ifelse(single_plot, 1, 3),
						cex = 0.5,
						single_plot = single_plot
					)
				}
			}
		}
		
		# TF plots
		if (!is.null(focus_tfs_fn)) {
			tfs = read.table(focus_tfs_fn, h = T)
			tfs = rownames(tfs)
			sc_2d = scp_plot_mat(
				sc_2d,
				png_fn = paste0(outdir, "/tf_all_mat.png"),
				genes_pool = tfs,
				per_clust_genes = fig_mat_per_clust_genes,
				gene_min_fold = fig_mat_gene_min_fold,
				gene_min_cov = fig_mat_gene_min_cov,
				width = fig_mat_width,
				height = fig_mat_height,
				text_cex = fig_mat_text_cex,
				smooth_n = fig_mat_smooth_n,
				blacklist = sc_2d@blist_genes
			)
			
			sc_2d = scp_plot_mat(
				sc_2d,
				plot_cells = F,
				genes_pool = tfs,
				png_fn = paste0(outdir, "/tf_clust_mat.png"),
				fp_shades = fold_shades,
				per_clust_genes = fig_mat_per_clust_genes,
				gene_min_fold = fig_mat_gene_min_fold,
				gene_min_cov = fig_mat_gene_min_cov,
				width = fig_mat_width,
				height = fig_mat_height,
				text_cex = fig_mat_text_cex,
				blacklist = sc_2d@blist_genes
			)
			
			cl_ord = sc_2d@cl_ord
			
			covgenes = names(rowSums(sc_cl@scmat@mat) > 80)
			tfs = intersect(tfs, covgenes)
			tfs = intersect(tfs, rownames(sc_cl@clust_fp))
			mask_folds = sc_cl@clust_fp[tfs, ] * (sc_cl@clust_gcov[tfs, ] > 0.5)
			spec_tfs = names(which(apply(mask_folds, 1, max) > 2))
			
			tf_cor = cor(log2(t(sc_cl@clust_fp[spec_tfs, ])))
			
			tf_hc = hclust(dist(tf_cor), "ward.D2")
			tf_max_at = apply(sc_cl@clust_fp[spec_tfs, cl_ord], 1, which.max)
			tf_hc = as.hclust(reorder(as.dendrogram(tf_hc), tf_max_at, agglo.FUN =
																	mean))
			png(paste0(outdir, "/tfs_hclust.png"),
					w = 1500,
					h = 1500)
			plot(as.phylo(tf_hc), type = "fan")
			dev.off()
			png(paste0(outdir, "/tfs_cor.png"),
					w = 2000,
					h = 2000)
			diag(tf_cor) = NA
			tfs_shades = colorRampPalette(
				c(
					"blue",
					"blue",
					"blue",
					"lightblue",
					"white",
					"red",
					"orange",
					"yellow",
					"black"
				)
			)(1000)
			
			image(
				tf_cor[tf_hc$order, tf_hc$order],
				col = tfs_shades,
				zlim = c(-1, 1),
				xaxt = 'n',
				yaxt = 'n'
			)
			mtext(
				spec_tfs[tf_hc$order],
				at = seq(0, 1, l = length(spec_tfs)),
				las = 2,
				side = 2,
				cex = 0.8
			)
			mtext(
				spec_tfs[tf_hc$order],
				at = seq(0, 1, l = length(spec_tfs)),
				las = 2,
				side = 4,
				cex = 0.8
			)
			mtext(
				spec_tfs[tf_hc$order],
				at = seq(0, 1, l = length(spec_tfs)),
				las = 2,
				side = 1,
				cex = 0.8
			)
			dev.off()
		}
		
		if (plot_genes_2d) {
			message("plotting genes 2d")
			if (!exists("spec_tfs")) {
				spec_tfs = rownames(sc_cl@feat_mat)
				if (!is.null(marker_colors)) {
					spec_tfs = unique(c(spec_tfs, rownames(marker_colors)))
				}
			}
			dir.create(paste0(outdir, "/genes/"))
			for (nm in spec_tfs) {
				cat("plotting ", nm, "\n")
				scp_plot_gene_2d(
					sc_2d,
					gene_nm = nm,
					base_dir = paste0(outdir, "/genes/"),
					w = 1200,
					h = 1200,
					reg_factor = 10
				)
			}
		}

		# plot gene modules
		pref = "^gene_modules_by_"
		gmod_objs = scdb_get_existing_obj_names(scdb, pattern=pref)
		
		for(gmod_obj in gmod_objs) {
			gmod_obj_nm = gsub(pref, "", gmod_obj)
			
			message(sprintf("generating gene modules plots for %s", gmod_obj_nm))
			
			gm = scdb_load(scdb, gmod_obj)
			
			scgm_plot_matrices(
				gm,
				base_dir = sprintf("%s/gmods/%s", outdir, gmod_obj_nm),
				name = paste0(gmod_obj_nm, "_"),
				w = max(800, length(gm@gmods)*8),
				h = max(800, length(gm@gmods)*8),
				mark_mods = T,
				cex_gene = 0.2, 
				res = min(max(220, length(gm@gmods)/5), 500)
			)
			
			scgm_plot_gmod_on_2d(
				sc_2d,
				gm,
				output_dir = sprintf("%s/gmods/%s/gmod_2d", outdir, gmod_obj_nm),
				name = paste0(gmod_obj_nm, "_"),
				width = 800,
				height = 800
			) 
				
			# plot breakdown of cells to modules by feature (e.g. patient, sample etc.)
			if (!is.null(get_param("scp_cell_mod_tab_by"))) {
				plot_cells_breakdown_to_mods_by_feature(sc_2d)
			}
		}
		
		scdb_save(scdb, sc_2d, "plot")
	}
	
	invisible(sc_2d)
	
}


#' sc_pipe_init_sort
#' 
#' Project sc_mat on given cell type markers, cluster projected mat. Should be followed by manual
#' assignment of cell types to cluster centers
#'
#' @param sc_mat 
#' @param markers_fn table with columns: subtype,gene,type
#' @param max_umi_count 
#' @param alg_type 
#' @param min_clust_size 
#' @param n_mark_clusts 
#' @param ordered_subtypes 
#' @param sc_cl 
#'
#' @return kmeans object of the projected mat
#' @export
#'
#' @examples
sc_pipe_init_sort = function(sc_mat,
														 markers_fn,
														 factor_by_median = get_param("scc_factor_by_median"),
														 norm_size = get_param("scc_norm_size"),
														 k_scale_umi = get_param("scc_k_scale_umi"),
														 alg_type = "knn",
														 min_clust_size = 20,
														 n_mark_clusts = 80,
														 ordered_subtypes = c(
														 	"ab_T_cell",
														 	"regulatory_CD4_T_cell",
														 	"CD4_T_cells",
														 	"CD8_T_cells",
														 	"gd_T_cell",
														 	"NK_cell",
														 	"cytotoxic_B_cells",
														 	"IgA_B_cell",
														 	"IgM_B_cell",
														 	"IgG_B_cell",
														 	"plasma_cells",
														 	"DC",
														 	"macrophage",
														 	"granulocytes",
														 	"erythrocyte",
														 	"osteoclast"
														 ),
														 sc_cl = NULL)
{
	sc_marker_clusts(
		sc_mat,
		markers_fn = markers_fn,
		factor_by_median = factor_by_median,
		norm_size = norm_size,
		k_scale_umi = k_scale_umi,
		n_mark_clusts = n_mark_clusts,
		alg_type=alg_type,
		min_clust_size = min_clust_size,
		analysis_dir = getwd(),
		ordered_subtypes = ordered_subtypes,
		totmark_cl=sc_cl
	)
	
}

#' sc_pipe_sort
#'
#' Assign cells to types by finding the nearest center per cell
#' @param sc_mat
#' @param sort_types_markers_fn
#' @param sort_ref_centers_pref_fn
#' @param out_base_dir
#' @param out_base_dir_add_timestamp
#' @param sc_mat_tables_fn_pref
#'
#' @return returns invisibly a the original tgScMat with assigned types in @cell_metadata$assignd_type
#' @export
#'
#' @examples
sc_pipe_sort = function(scdb,
												sc_mat=NULL,
												sort_types_markers_fn,
												sort_ref_centers_pref_fn,
												factor_by_median = get_param("scc_factor_by_median"),
												norm_size = get_param("scc_norm_size"),
												k_scale_umi = get_param("scc_k_scale_umi"),
												out_base_dir,
												out_base_dir_add_timestamp = F
)
{
	if (is.null(sc_mat)) {
		sc_mat = scdb_load(scdb, "clean_mat")
	}
	
	stopifnot(!is.null(sc_mat))
	
	if (out_base_dir_add_timestamp) {
		out_base_dir = paste0(out_base_dir, "/", format(Sys.time(), "%Y.%m.%d_%H.%M"))
	}
	
	if (!dir.exists(out_base_dir)) {
		dir.create(out_base_dir)
	}
	
	# in-silico sort to types by reference type markers and cluster centers
	message(
		sprintf(
			"in-silico sorting cells based on cell type markers: %s and pre-defined clusters %s",
			sort_types_markers_fn,
			sort_ref_centers_pref_fn
		)
	)
	
	sc_mat = sc_marker_split(
		sc_mat,
		sort_types_markers_fn,
		sprintf("%s.txt", sort_ref_centers_pref_fn),
		sprintf("%s_assign.txt", sort_ref_centers_pref_fn),
		factor_by_median = factor_by_median,
		norm_size = norm_size,
		k_scale_umi = k_scale_umi,
		fig_pref = sprintf("%s/sort", out_base_dir),
		out_base_dir
	)
	
	invisible(sc_mat)
}
