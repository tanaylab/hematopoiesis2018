require("KernSmooth")
require("reshape2")
require("RANN")
require("plyr")
require("plotrix")
require("gplots")

library("compositions")

sc_pipeline = function(basedir, outdir, batches, index_fn, umi_dir,
	batch_meta_attr = "Amp_batch_ID", amb_epsilon = 0.03, mark_blacklist_terms = c(),
	T_edge_2dproj = 0.05, force_max_deg_2dproj = 8, read_and_clean_only = F, cells = NULL) {

	dir.create(basedir)
	dir.create(outdir)

	scdb = scdb_init(basedir = basedir)

	sc_raw_mat = sc_pipe_build_mat(scdb,
        index_fn = index_fn,
        batch_meta_attr = batch_meta_attr,
        base_dir= umi_dir,
        mars_batches = batches,
        outdir = outdir)

	sc_clean_mat = sc_pipe_mat_to_clean_mat(scdb, sc_raw_mat,
        batch_meta_attr = batch_meta_attr,
        filt_amb_on_clusts = T,
        filt_outliers_on_clusts = T,
        remove_outliers_before_clustering = F,
        mark_blacklist_terms = mark_blacklist_terms,
        amb_epsilon = amb_epsilon)
	
	if (!read_and_clean_only) {
		if (!is.null(cells)) {
			sc_clean_mat@cells = cells; sc_clean_mat@ncells = length(cells)
			sc_clean_mat@mat = sc_clean_mat@mat[,cells]
		}
		sc_cl = sc_pipe_mat_to_cmods(scdb, sc_clean_mat,
	        mark_blacklist_terms = mark_blacklist_terms,
	        outdir = outdir,
	        filt_outliers_on_clusts = T,
	        clust_fp_metadata_fields = NA,
	        tab_clust_fp_fn = NA)

		sc_2d = sc_pipe_plots(scdb, sc_cl,
	      	outdir = outdir,
	        mark_blacklist_terms = mark_blacklist_terms,
	        T_edge_2dproj = T_edge_2dproj,
	        force_max_deg_2dproj = force_max_deg_2dproj)
	}
	scdb
}

perform_bootstrap = function(scdb, n, outdir, mark_blacklist_terms, boot_ratio = 0.7, T_edge_2dproj = 0.05, force_max_deg_2dproj = 8, 
	n_procs = 16, analysis_only = F) {
	sc_cl = sc_pipe_mat_to_cmods(scdb)
	filename = paste0(scdb@basedir, "/scrdb_data_cocluster.Rda")
	if (!file.exists(filename)) {
	   cc = scc_bootstrap_clusts(sc_cl,
        	scdb = scdb,
		k_cut = get_param("scc_bootstrap_k_cut"),
        	boot_ratio = boot_ratio,
		n_procs = n_procs,
        	min_clust_size = 20,
        	verbose = T,
        	rseed=1)

	   save(cc, file = filename)
	} else { load(filename)}
	if (!analysis_only) {
		bs_cl = scc_coclust_to_cell_modules(scdb, sc_cl, cc, n, 20, method="hclust")
		file.remove(paste0(scdb@basedir, "/scrdb_data_plot.RData"))
		sc_2d = sc_pipe_plots(scdb, bs_cl,
        		outdir = outdir,
	        	mark_blacklist_terms = mark_blacklist_terms,
	        	T_edge_2dproj = T_edge_2dproj,
	        	force_max_deg_2dproj = force_max_deg_2dproj)
	}
}

plot_bootstrapping = function(scdb) {

}

score_on_gene_prograns = function(umis, gene_map, thresh = 0.5) {

	genes = intersect(rownames(umis), names(gene_map))
	umis = umis[genes,]; gene_map = gene_map[genes]
	cell_modules = apply(umis, 2, tapply, gene_map, sum)
	cell_modules = cell_modules[ -nrow(cell_modules),]
	below_med = cell_modules <= apply(cell_modules,1,quantile, thresh)
	cell_bg = apply(matrix(rownames(cell_modules)),1, function(x) 
		{ifelse(sum(umis[,below_med[x,]]) > 0, sum(cell_modules[x, below_med[x,]]) / sum(umis[,below_med[x,]]), 0)})
	names(cell_bg) = rownames(cell_modules)
	cell_uc = colSums(umis)
	zscore = apply(matrix(rownames(cell_modules)),1,
               function(x) { (cell_modules[x,] - cell_bg[x] * cell_uc) / sqrt(cell_bg[x] * (1 - cell_bg[x]) * cell_uc)})

	zscore
}

xy_scatter_genes = function(x,y, bad_genes = c(), disp_genes = c(), cols = c("navyblue", "chocolate2"), text = F, fc = 1, reg = 10, lwd = 1) {
        good_genes = setdiff(names(x), bad_genes)
	lim = log2(c(reg, reg + max(c(x,y))))
        plot(log2(x[good_genes] + reg), log2(y[good_genes] + reg), pch = 20, cex = 2, col = cols[1 + good_genes %in% disp_genes],
        xlim = lim, ylim = lim, axes = F, xlab = "", ylab = "")
        if (text) { text(log2(x[disp_genes] + reg), log2(y[disp_genes] + reg), disp_genes)}
        abline(coef = c(fc,1), lty = 2); abline(coef = c(-fc,1), lty = 2)
        axis(1); axis(2)
}

choose_genes_from_clust = function(sc_cl, good_clusts = colnames(sc_cl@clust_fp), 
	nms_per_clust = 5, nms_thresh = 5, max_num = Inf, bad_genes = c(), must_haves = c(),
	ord = "none") {

	lfp = log2(sc_cl@clust_fp[,good_clusts])
	nms = unique(as.vector(apply(lfp,2,function(x){ head(names(sort(x, decreasing = T)),nms_per_clust)})))
	nms = setdiff(nms, c(bad_genes, names(which(apply(lfp[nms,],1,max) < nms_thresh))))
	nms = union(must_haves, head(names(sort(apply(lfp[nms, ], 1, max),T)), max_num - length(must_haves)))
	if (ord == "hc") {
		nms = nms[ hclust(dist(cor(t(lfp[nms,]))), "ward.D2")$order]
	} else if (ord == "max.col") {
		nms = nms[ order(max.col(lfp[nms,]), rowSums(as.matrix(sc_cl@scmat@mat[nms,])))]
	}
	nms
}

sc_to_bulk = function(sc_cl, comb, bad_genes = c(), cells = names(sc_cl@clusts), min_comb = 0, choose_genes = T) {
	umis = as.matrix(sc_cl@scmat@mat)
	umis_n = sweep(umis,2,colSums(umis),"/") * 1000

	MAP = as.numeric(factor(comb[cells])); names(MAP) = cells
	if (choose_genes) {
		genes = setdiff(scr_chi_square_diff_genes(umis[,cells], MAP = MAP[cells], fdr = T, pval = 1e-3), bad_genes)
	} else {
		genes = rownames(umis)
	}
	m = t(apply(umis_n[genes, cells],1,tapply,comb[cells],sum))
	sizes = table(comb[cells])
	good_clusts = names(which(sizes >= min_comb))
	m = m[,good_clusts]; sizes = sizes[good_clusts]
	m = sweep(m,2,as.vector(sizes),"/") * min(sizes)
	m
}

create_batch_matrix = function(sc_cl, comb, fname, batch_iden = "Amp_batch_ID", cord = colnames(sc_cl@clust_fp),
	batch_shades = colorRampPalette(c("white", "navyblue"))(1000)) {
	cells = names(comb)
	cell_stats = sc_cl@scmat@cell_metadata[cells,]
	batch_meta = unique(cbind(as.vector(cell_stats[,batch_iden]), comb))
	B = batch_meta[,1]
	batch_dist = table(cell_stats[,batch_iden], sc_cl@clusts)
	batch_dist = batch_dist[B,]
	batch_dist = sweep(batch_dist,2,colSums(batch_dist),"/")
	bord = order(batch_meta[,2])
	bcls = cumsum(table(batch_meta[,2])) / nrow(batch_meta)
	png(fname, height = 2000, width = 3000)
	par(mar = c(3,20,0,10))
	image(t(batch_dist[ bord, cord]), col = batch_shades, axes = F, zlim = c(0, 1))
        mtext(batch_meta[bord,1], side = 4, las = 2, at = (1 - seq_len(nrow(batch_meta))) / (1 - nrow(batch_meta)))
	mtext(names(bcls), side = 2, las = 2, at = rowMeans(cbind(c(0,bcls[-length(bcls)]), bcls)))
	mtext(cord, side = 1, las = 1, at = (1 - seq_along(cord)) / (1 - length(cord)))
	#abline(h = bcls, lwd = 3)
	dev.off()
}

plot_sc_heatmap = function(sc_cl, nms, good_clusts = colnames(sc_cl@clust_fp), cells = NULL, fname = NULL, annotate = F, mar = rep(0,4),
	genes_shades = colorRampPalette(c("white", "orange", "tomato","mediumorchid4", "midnightblue"))(1000), draw_cls = T) {

	if (is.null(cells)) {
		cells = names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% good_clusts))
	}
	umis = as.matrix(sc_cl@scmat@mat)
	umis_n = sweep(umis,2,colSums(umis),"/") * 1000
	foc = log(1 + 7 * umis_n)
	cls = cumsum(table(factor(sc_cl@clusts[cells], levels = good_clusts))) / length(cells)	
	if (!is.null(fname)) { 
		png(fname, height = 2000, width = 3300)
		par(mar = mar)
	}
	cell_ord = cells[order(factor(sc_cl@clusts[cells], levels = good_clusts), sample(length(cells)))]
	IM = foc[nms, cell_ord]
	image(t(IM), col = genes_shades, axes = F)
	zlim = quantile(IM, c(0,1))
	message("zlim: ", zlim[1], " - ", zlim[2])
	if (annotate) {
		mtext(nms, side = 2, las = 2, at = (1 - seq_along(nms)) / (1 - length(nms)))
		mtext(names(cls), side = 1, las = 2, at = rowMeans(cbind(c(0,cls[-length(cls)]), cls)))
	}
	if (draw_cls) {
		abline(v = cls, lty = 2)
	}
        if (!is.null(fname)) {
		dev.off()
	}
	cell_ord
}

reposition_cc = function(sc_2d, margin = 500, coords = NULL) {
        xlim = quantile(sc_2d@x, c(0,1)) + c(-margin, margin)
        ylim = quantile(sc_2d@y, c(0,1)) + c(-margin, margin)
	P = .graph_con_comp(sc_2d@clust_graph)
	cols = rainbow(max(P))
	if (is.null(coords)) {
		plot(sc_2d@x, sc_2d@y, pch = 21, bg = cols[P[sc_2d@scl@clusts]], xlim = xlim, ylim = ylim)
        	coords <- unlist(locator(2, type="l"))
	}
	nearest_cc = P[sc_2d@scl@clusts[names(which.min( rowSums(cbind(sc_2d@x - coords[1], sc_2d@y - coords[3]) ^ 2)))]]
        new_2d = sc_2d
	new_2d@x = new_2d@x + (P[new_2d@scl@clusts] == nearest_cc) * (coords[2] - coords[1])
        new_2d@y = new_2d@y + (P[new_2d@scl@clusts] == nearest_cc) * (coords[4] - coords[3])
	new_2d@x_cl = tapply(new_2d@x, new_2d@scl@clusts, mean)
	new_2d@y_cl = tapply(new_2d@y, new_2d@scl@clusts, mean)
        list(sc_2d = new_2d, coords = coords)
}

scr_write_models_file = function(cl, filename = "models.txt") {
	umis = as.matrix(cl@scmat@mat)
        genes = rownames(cl@clust_fp)
        write.table(data.frame(log2(cl@clust_fp), umicount = rowSums(umis[genes,])), col.names = NA, quote = F, sep = "\t", file = filename)
}

proj_sc_on_clusts = function(umis, new_cells, old_cells, clusts, markers, knn_for_cor = 10, K = 10) {

	umis_n = sweep(umis,2,colSums(umis),"/") * 1000
	ranks = rev(1 - (0:knn_for_cor+1)/(knn_for_cor + 2))
	d_old = t(log2(1+7*umis_n))[old_cells ,markers]
	d_new = t(log2(1+7*umis_n))[new_cells ,markers]
	dcor = as.matrix(cor(t(d_new), t(d_old)))

	clusts = clusts[old_cells]
	clusts1 = diag(max(clusts))[,clusts]    # 1 encoding of the clusters, rows= clusters, cols =nodes
	csizes = rowSums(clusts1); csizes[csizes==0] = 1 # add 1 to avoid division by 0
	clusts1 = clusts1/csizes

	m_knn = t(apply(dcor, 1, function(x) ranks[ rank(-x)])) #1 - pmin(rank(-x)/knn_for_cor,1) ))
	m_knn[ is.na(m_knn)] = 0

	t_knn = t(apply(t(dcor), 1, function(x) ranks[ rank(-x)])) #1 - pmin(rank(-x)/knn_for_cor,1) ))
	t_knn[ is.na(t_knn)] = 0

	forward_votes = t(clusts1 %*% t(m_knn)) *1000
	backward_votes = tcrossprod(t(t_knn), clusts1) * 1000
	votes = backward_votes * forward_votes + 1e-10 * forward_votes  # + 1e-10 * forward_votes is for nodes that no one wants
	clusts = apply(votes,1,which.max)
	names(clusts) = rownames(d_new)
	clusts
}

gate_facs = function(facs, xlab, ylab, gate = TRUE, colgate = NULL, roof = 1e4, rect = FALSE, log = "", polygon = NULL, n = 1000) {

  if (is.null(polygon)) {
    disp_facs = facs[gate, ]
    if (roof < dim(disp_facs)[1]) {
      ind = sample(dim(disp_facs)[1], roof)
    } else {
      ind = 1:dim(disp_facs)[1]
    }

    x = disp_facs[ind, xlab]
    y = disp_facs[ind, ylab]

    if (!is.null(colgate)) {
      colgate = colgate[gate, ]
      colgate = colgate[ind, ]
      colors = unbinary(apply(matrix(1:dim(colgate)[1]), MARGIN = 1, FUN = function(x){paste(colgate[x,] * 1, collapse = "")}))
      #smoothScatter(x, y)
      #points(x, y, col = colors, pch=21, cex = 0.4)
      plot(x, y, pch = 21, cex = 0.4, log = log, col = colors + 1)
    } else {
      plot(x, y, pch = 21, cex = 0.4, log = log)
    }

    coords <- locator(n, type="l") # add lines
   C = unlist(coords)
    if (!rect) {
      n = length(coords$x)
      c = rep(0, n*4)
      dim(c) = c(n*2,2)
      c[1:n,] = C
      c[(n+1):(n*2),] = C
    } else {
      n = 4
      c = matrix(0,n*2,2)
      c[1,] = c(C[1], C[3])
      c[2,] = c(C[1], C[4])
      c[3,] = c(C[2], C[4])
      c[4,] = c(C[2], C[3])
      c[5:8,] = c[1:4,]
    }
    lines(c, col = "red")
  } else {
    c = polygon
    n = nrow(c) / 2
  }

  x = facs[, xlab]
  y = facs[, ylab]

  all_pos = apply(matrix(1:n), MARGIN = 1, FUN =
                    function(i){sign( (c[i + 1,1]- c[i,1])*(y-c[i,2]) - (c[i + 1,2]-c[i,2])*(x-c[i,1]) )})
  w = unlist(all_pos)
  within = apply(w < 0, 1, prod) == 1

  return(list(gate = within, polygon = c))

}

scr_chi_square_diff_genes = function(umis, MAP = NULL, g1 = NULL, g2 = NULL, pval, fdr = F) {

  if (is.null(MAP)) {
    MAP = c(rep(1,length(g1)), rep(2, length(g2)))
    names(MAP) = c(g1, g2)
  }
  cells = names(MAP)
  umis = umis[,cells]
  uniform_a = rowSums(umis)/sum(umis)
  exp_count = matrix(uniform_a, ncol = 1) %*% matrix(colSums(umis),1) # exp_counts per cell
  dimnames(exp_count)  = dimnames(umis)
  ex = t(daply(.data= data.frame(cbind(V1 = MAP, t(exp_count)), check.names = F), .(V1), colSums))[-1,]
  obs = t(daply(.data= data.frame(cbind(V1 = MAP, t(umis)), check.names = F), .(V1), colSums))[-1,]

  x2 = rowSums(((obs-ex)^2 )/ex ) # chi^2 with df = ncluster-1

  if (!fdr) {
    sig_genes = x2 > qchisq(1-pval,df= length(unique(MAP)) - 1)
  } else {
    pvals = p.adjust(1 - pchisq(x2, df = length(unique(MAP)) - 1), "fdr")
    sig_genes = pvals < pval
  }
  sig_genes[ is.na(sig_genes)] = F
  return (names(sig_genes)[sig_genes])

}


scr_create_gene_modules_table = function(umis, clusts, clusters, K = 20, Z = 2) {

  cluster_genes = matrix(FALSE, nrow = nrow(umis), ncol = length(clusters))
  colnames(cluster_genes) = clusters
  rownames(cluster_genes) = rownames(umis)

  idx_n = t(t(umis) / colSums(umis))
  idx_n = idx_n[, clusts %in% clusters]

  m = aggregate(x = t(idx_n), by = list(clusts[clusts %in% clusters]), FUN = mean)

  M = t(m[,-1])
  colnames(M) = m[,1]
  means = rowMeans(M)

  for (c in clusters) {
    cluster_genes[,as.character(c)] = M[, as.character(c)] > Z * means & rowSums(umis[, clusts == c]) > K
  }

  return (cluster_genes)
}


scr_find_gene_modules = function(umis, nmodules = 30, min_gene_umis = 5, min_cell_umis = 500,
                                 min_var_mean = 1.2, min_cor_within_module = 0.05, rseed = 1) {
  set.seed(rseed)
  u = apply(umis[,colSums(umis)>min_cell_umis],2,.downsamp_one,min_cell_umis)
  rownames(u) = rownames(umis)
  vm = apply(u, 1, var) / rowMeans(u)
  genes = rowSums(u)>min_gene_umis & vm >= min_var_mean
  u = u[genes,]
  u = log2(u+1)
  #   u = u-rowMeans(u)
  corr = cor(t(u))
  hc = hclust(as.dist(1-corr), method = "ward.D2")
  cls = cutree(hc, nmodules)
  #   cor_means = laply(1:nmodules, function(x) {c = corr[cls==x,cls==x]; mean(c[lower.tri(c)])})
  modules = list()
  for(k in 1:nmodules) {
    c = corr[cls==k,cls==k];
    if (mean(c[lower.tri(c)]) > min_cor_within_module) {
      modules[[length(modules)+1]] = names(cls[cls==k])
    }
  }
  return(list(modules = modules, cls = cls, corr = corr))
}


scr_find_outline = function(clusts_2d, reg = 0.7, bw=50, cells = names(clusts_2d@x)) {
	#cells_graph_x = clusts_2d@x[ !is.nan(clusts_2d@x)]
        #cells_graph_y = clusts_2d@y[ !is.nan(clusts_2d@y)]
	xl = c(min(clusts_2d@x[cells], na.rm=T), max(clusts_2d@x, na.rm=T))
        yl = c(min(clusts_2d@y[cells], na.rm=T), max(clusts_2d@y, na.rm=T))

        rngx = xl[2]-xl[1]
        rngy = yl[2]-yl[1]
        pt_reg = merge(seq(xl[1],xl[2],length.out=bw)
                        ,seq(yl[1],yl[2],length.out=bw))
	pt_reg = matrix(unlist(pt_reg),ncol=2)
        back <- bkde2D(rbind(pt_reg, cbind(clusts_2d@x[cells], clusts_2d@y[cells])),
                        b=c(rngx/bw, rngy/bw),
                        gridsize=c(500, 500))

        outline = matrix(0, nrow = 500, ncol = 500)
        Z = back$fhat > quantile(back$fhat, reg)
        outline[ Z[,c(500, 1:499)] != Z] = 1
	outline = melt(outline)

	outline[,1] = back$x1[ outline[,1]]
        outline[,2] = back$x2[ outline[,2]]
        outline = outline[ outline[,3] == 1,1:2]
        neigh =nn2(outline,k = 20)

        ord = c(); i = 1;
        while(length(ord) != nrow(outline)){
                ord = c(ord,i);
                i = neigh$nn.idx[i,which(!(neigh$nn.idx[i,] %in% ord))[1]]
        }
        ord = ord[ !is.na(ord)]
        outline[ord,]
}

plot_gene_2d = function(scl2d,
        gene_nm, w, h, base_dir = NULL, rna_mat = NULL,
        top_color_q = 0.95,
        cex=1.5,
        n_reg_grid=50,
        reg_factor=1,
        low_umi=0,
        mid_umi=1,
        cont_levels = 0,
        bw_bins = 50,
        positive_psize = 0.5,
        negative_psize = 0.5,
        min_rna_to_color = 0,
        outline = NULL,
        return_values = F,
        rna_shades = colorRampPalette(c("white", "white", "lightgray", "darkorange1", "darkgoldenrod1", "darkgoldenrod4")),
        pt_shades = colorRampPalette(c("white", "white", "lightgray", "darkgray", "orange", "burlywood1", "chocolate4")),
        bg_shades = colorRampPalette(c("white", "white", "lightgray", "darkgray", "cyan", "blue1", "blue3")),
        bg_cells = c())
{
        graph_x = scl2d@x
        graph_y = scl2d@y
        if(is.null(rna_mat)) {
                rna = scl2d@scl@scmat@mat[gene_nm,names(graph_x)]
                rna_tot = colSums(as.matrix(scl2d@scl@scmat@mat)[,names(graph_x)])
        } else {
                rna = rna_mat[gene_nm,names(graph_x)]
                rna_tot = colSums(rna_mat[,names(graph_x)])
        }
        med_n = median(rna_tot)
        rna_base = floor(med_n*rna/rna_tot)
        rna = rna_base + sapply(med_n*rna/rna_tot-rna_base, function(p) rbinom(1,1,p))
        rna[is.nan(rna)] = 0

        xl = c(min(graph_x,na.rm=T), max(graph_x, na.rm=T))
        yl = c(min(graph_y,na.rm=T), max(graph_y, na.rm=T))
        rngx = xl[2]-xl[1]
        rngy = yl[2]-yl[1]
        pt_reg = merge(seq(xl[1],xl[2],length.out=n_reg_grid)
                        ,seq(yl[1],yl[2],length.out=n_reg_grid))
        pt_reg = matrix(unlist(pt_reg),ncol=2)
        epsilon=0
        if(!is.null(base_dir)) {
                fnm = gene_nm
                if(nchar(fnm) > 15) {
                        fnm = substr(gene_nm, 1, 15)
                }
                fnm = sub("/", "_", gene_nm)
                fn = sprintf("%s/%s.png", base_dir, fnm)
                png(fn, heigh=h, w=w);
                par(mar=c(0,0,0,0))
        }
        rna_x = unlist(apply(cbind(graph_x,rna),1,function(x) rep(x[1],each=x[2]*reg_factor)))
        rna_y = unlist(apply(cbind(graph_y,rna),1,function(x) rep(x[1],each=x[2]*reg_factor)))
        cells = sapply(strsplit(names(rna_x), "\\."), "[", 1)
        fore <- bkde2D(rbind(pt_reg, cbind(rna_x, rna_y)[!(cells %in% bg_cells),]),
                b=c(rngx/bw_bins, rngy/bw_bins),
                gridsize=c(500, 500))
        if (length(bg_cells) > 0) {
           bg_rna = cbind(rna_x, rna_y)[cells %in% bg_cells,]
        } else {
           bg_rna = cbind(graph_x, graph_y)
        }
        back <- bkde2D(rbind(pt_reg, bg_rna),
                b=c(rngx/bw_bins, rngy/bw_bins),
                gridsize=c(500, 500))
        fore$fhat = fore$fhat * sum(back$fhat)/sum(fore$fhat)
#               smoothScatter(rna_x, rna_y, colramp=rna_col, xlim=xl, ylim=yl)
        lrs = log2(fore$fhat/back$fhat)
        if (length(bg_cells) > 0) {
           zlim = max(abs(quantile(lrs))); zlim = c(-zlim, zlim)
        } else {
                if(median(lrs,na.rm=T)>-1) { #background regul is too dominant
                        lrs = lrs - median(lrs,na.rm=T) - 1
                }
                zlim = c(-1,max(4, max(lrs)))
        }
        message("plot ", gene_nm, " tot ", sum(rna))

#        if (is.null(outline)) {
           xlim = quantile(back$x1, c(0,1))
           ylim = quantile(back$x2, c(0,1))
#        } else {
#           xlim = quantile(outline[,1], c(0,1))
#           ylim = quantile(outline[,2], c(0,1))
#        }
        image(x=back$x1, y=back$x2, z=lrs,zlim=zlim, col=rna_shades(1000),xaxt='n', yaxt='n', xlim = xlim, ylim = ylim)
        if(cont_levels > 0) {
                contour(x=back$x1, y=back$x2, z=lrs, nlevels=cont_levels, lwd=3,add=T, drawlabels=F)
        }

        low_umi_g = ifelse(low_umi == -1, floor(median(rna)), low_umi)
        high_rna_t = as.numeric(quantile(rna, top_color_q))
#               pt_cols = pt_shades[ifelse(rna <= low_umi_g, 1, ifelse(rna <= mid_umi, 2, ifelse(rna <= high_rna_t, 3,4)))]
        pt_val = round(1+999/length(rna) * rank(pmax(rna,min_rna_to_color), ties.method="min"))
        pt_cols = ifelse(names(rna) %in% bg_cells, bg_shades(1000)[pt_val], pt_shades(1000)[pt_val])

#                       cex=ifelse(pt_cols=="white", 0.5, 1),
        points(graph_x,
                graph_y, pch=21,
                bg=pt_cols,
                cex=ifelse(rna>0, positive_psize, negative_psize),
                col=ifelse(pt_cols=="white", "black", "black"),
                lwd=0.5)
        if (is.null(outline)) {
                grid(lwd=3,col="black",lty=1)
        } else {
          points(outline[,1], outline[,2], type = "l", lwd = 4)
        }
        if(!is.null(base_dir)) {
                dev.off()
        }
        if (return_values) {
           return (list(rna_cols = pt_cols, back = back, lrs = lrs))
        }
}

proj_cells_on_clust_graph = function(sc_2d, clusts, blur = 0, omit_clust=-1, nn.idx, main = T, use_top_k = -1)
{
        if(use_top_k == -1) {
                use_top_k = ncol(nn.idx)
        }

        x_cl = sc_2d@x_cl
        y_cl = sc_2d@y_cl

        blurx = blur*(max(x_cl) - min(x_cl))
        blury = blur*(max(y_cl) - min(y_cl))

        omit_x = NA
        omit_y = NA
        omit_clust = as.character(omit_clust)
        clusts = as.character(clusts)
        if(omit_clust != -1) {
                omit_x = x_cl[omit_clust]
                omit_y = y_cl[omit_clust]
                x_cl[omit_clust] = NA
                y_cl[omit_clust] = NA
        }

        px = apply(nn.idx[,1:use_top_k], 1,
                                function(x) ifelse(sum(clusts[x]!=omit_clust)>0,
                                                mean(x_cl[clusts[x]], na.rm=T),
                                                omit_x))
        py = apply(nn.idx[,1:use_top_k], 1,
                                function(x) ifelse(sum(clusts[x]!=omit_clust)>0,
                                                mean(y_cl[clusts[x]], na.rm=T),
                                                omit_y))


        message("Blur x ", blurx, " y ", blury)
	px = px + rnorm(mean=0, sd=blurx, n=length(px))
        py = py + rnorm(mean=0, sd=blury, n=length(py))
        list(x = px, y = py)
}

proj_ds_on_graph = function(sc_2d, K = 50, umis, fn = "new_umis.png", bg_cells = NULL, coords = NULL, knn_for_cor = 100,
	     outline = NULL, bw = 50, cex = 1.5, reg = 10, bg_reg = 1, lwd = 1, 
	     clust_shades = colorRampPalette(c("gray88", "orange", "red"))(101), markers = rownames(sc_2d@scl@feat_mat))
{

	sc_cl = sc_2d@scl
	old_umis = as.matrix(sc_cl@scmat@mat)
	# create gmod_fp for new ds
	cells = names(sc_cl@clusts)
	clust_names = names(sc_2d@x_cl)
	if (is.null(coords)) {

		# compute knn
		d_old = t(log2(1+7*old_umis[markers, cells]))
		d_new = t(log2(1+7*umis[markers,]))
                dcor = as.matrix(cor(t(d_new), t(d_old)))
		m_knn = t(apply(dcor, 1, function(x) 1 - pmin(rank(-x)/knn_for_cor,1) ))
		nn.idx = t(apply(m_knn, 1, function(x) as.numeric(factor(names(tail(sort(x),n=K)), levels = cells))))
		ass_clust = matrix(sc_cl@clusts[cells[nn.idx]], ncol = K, dimnames = dimnames(nn.idx))
       		clust = as.numeric(apply(ass_clust,1,function(x){ clust_names[ which.max(table(factor(x, levels = clust_names)))]}))
       		names(clust) = colnames(umis)
	# compute coordinates
		coords = proj_cells_on_clust_graph(sc_2d, sc_cl@clusts, 0.02, nn.idx = nn.idx, main = F)	
	}

	x = sc_2d@x
      	y = sc_2d@y
	x_cl = sc_2d@x_cl #nodeRenderInfo(scr_clust_render)$nodeX
	y_cl = sc_2d@y_cl #nodeRenderInfo(scr_clust_render)$nodeY

	if (!is.null(bg_cells)) {
	  wfg = !(colnames(umis) %in% bg_cells)
	  fg_coords = list(x = coords$x[wfg], y = coords$y[wfg])
	  bg_coords = list(x = coords$x[!wfg], y = coords$y[!wfg])
	  fg_clust = clust[ wfg]; bg_clust = clust[ !wfg]
	} else {
	  wfg = T
	  fg_coords = coords
	  bg_coords = list(x = sc_2d@x, y = sc_2d@y)
	  fg_clust = clust; bg_clust = sc_cl@clusts
	}

	xl = c(min(x, na.rm = T), max(x, na.rm = T))
   	yl = c(min(y, na.rm = T), max(y, na.rm = T))
        rngx = xl[2]-xl[1]
	rngy = yl[2]-yl[1]
        pt_reg = merge(seq(xl[1],xl[2],length.out=bw)
                        ,seq(yl[1],yl[2],length.out=bw))
        pt_reg = matrix(unlist(pt_reg),ncol=2)
	fore <- bkde2D(rbind(pt_reg, 
			cbind(rep(bg_coords$x, each = bg_reg), rep(bg_coords$y, each = bg_reg)),
	     	 	cbind(rep(fg_coords$x, each = reg), rep(fg_coords$y, each = reg))),
                        b=c(rngx/bw, rngy/bw),
                        gridsize=c(500, 500))
        back <- bkde2D(rbind(pt_reg, 
	     		cbind(rep(bg_coords$x, each = bg_reg), rep(bg_coords$y, each = bg_reg))),
                        b=c(rngx/bw, rngy/bw),
                        gridsize=c(500, 500))
        fore$fhat = fore$fhat * sum(back$fhat)/sum(fore$fhat)
        lrs = log2(fore$fhat/back$fhat)
	lrs[ is.nan(lrs)] = 0
        if(median(lrs, na.rm = T)>-1) {
       #                 lrs = lrs - median(lrs, na.rm = T) - 1
        }
        lrs = pmax(pmin(lrs,4,  na.rm = T),-1,  na.rm = T)
        
	gridx = findInterval(coords$x, seq(min(x, na.rm=T), max(x, na.rm=T), length.out=500), rightmost.closed=T)
	gridy = findInterval(coords$y, seq(min(y, na.rm=T), max(y, na.rm=T), length.out=500), rightmost.closed=T)
	lrs_range = (lrs - min(lrs)) / (max(lrs) - min(lrs)) * 100

	png(fn, height = 2000, width = 2000)
	plot(x, y, axes = F, xlab = "", ylab = "", type = "n")
	points(coords$x, coords$y, pch = 20, cex = ifelse(wfg, cex,0), lwd = 2,
	    col = clust_shades[ 1 + lrs_range[ cbind(gridx, gridy)]])
	if (!is.null(outline)) {points(outline[,1], outline[,2], type = "l", lwd = 4)}
	dev.off()
	return(list( coords = coords, clust = clust))
}

plot_virtual_facs = function(plate_facs, xlab, ylab, fname, filter = T, gates = list()) {
	png(fname, height=1500, width=1500)
	short_facs = plate_facs[ filter & plate_facs[,xlab] > 0 & plate_facs[,ylab] > 0 & !is.na(plate_facs[,xlab]) & !is.na(plate_facs[,ylab]),]
	k = bkde2D(cbind(log10(short_facs[,xlab]), log10(short_facs[,ylab])), bandwidth = 0.1)
	plot(log10(short_facs[,xlab]), log10(short_facs[,ylab]), pch = 20, cex = 1.5, col = "gray", log = "", xlim = c(1,5), axes = F, xlab = "", ylab = "")
	axis(1, at = 1:5, labels = 10^(1:5))
	axis(2, at = 1:5, labels = 10^(1:5))
	contour(k$x1, k$x2, k$fhat, nlevels = 20, drawlabels=F, lwd = 5, add = T)
	lapply(gates, function(x) lines(log10(x), lwd = 5, col = "red"))
	dev.off()
}

detect_batchy_genes = function(sc_mat, group_by, marks = NULL, batch_meta_attr = "Amp_batch_ID", cell_per_batch = 384) {
	if (is.null(marks)) {	
		gstat = umi_gene_stat(sc_mat)
		marks = select_markers(sc_mat,gstat = gstat,
                  type="any",
                  mark.min_var_mean = 0.2,
                  mark.sz_cor_norm_max=-Inf,
                  mark.niche_T = Inf)
	}
	umis = as.matrix(sc_mat@mat)
	umis_n = sweep(umis,2,colSums(umis),"/") * 1000
	cell_stats = sc_mat@cell_metadata
	m = t(apply(umis_n[marks,] > 1,1,tapply, cell_stats[, batch_meta_attr], sum))
	m = sweep(m,2, table(cell_stats[, batch_meta_attr]), "/") * cell_per_batch
	batch_temp = unique(cbind(as.vector(cell_stats[,c(batch_meta_attr, group_by)]),1))
	batch2area = apply(cbind(batch_temp[,group_by], "") ,1, paste0,collapse="."); names(batch2area) = batch_temp[,batch_meta_attr]
	mm = t(apply(m,1,tapply,batch2area,mean)); m2 = log2((m + 10) / (mm[,batch2area] + 10))
	m2 = log2((m + 10) / (apply(mm,1,mean) + 10))
	m2
}
