
#' @title Plotting object for a clustering solution
#'
#' @description Plotting object for a clustering solution
#'
#' Representing a 2d layout for a cluster model, used for plotting.
#'
#' \code{tgScMatClustLayout2D} represents a 2D layout for a cluster model.
#'
#' @slot scl tgScMatClust.
#' @slot cl_ord vector.
#' @slot c_ord vector.
#' @slot cl_colors vector.
#' @slot c_colors vector.
#' @slot mark_colors vector.
#' @slot mark_groups vector
#' @slot current_subset vector.
#' @slot x_cl vector.
#' @slot y_cl vector.
#' @slot clust_graph matrix.
#' @slot x vector.
#' @slot y vector.
#' @slot cell_2d_sd vector.
#' @slot blist_genes vector.
#'
#### @export tgScMatClustLayout2D
#### @exportClass tgScMatClustLayout2D
tgScMatClustLayout2D <- setClass(
        "tgScMatClustLayout2D",
	slots = c(scl = "tgScMatClust",
		  mat_k_confusion = "table",
		  mat_bal_confusion = "matrix",
		  mknnio_ij = "matrix",
		  k_confusion = "numeric",
		  cl_ord = "vector",
		  c_ord = "vector",
		  cl_colors = "vector",
		  c_colors = "vector",
		  mark_colors = "vector",
		  mark_groups = "vector",
		  current_subset = "vector",
		  x_cl = "vector",
		  y_cl = "vector",
		  clust_graph = "matrix",
		  x = "vector",
		  y = "vector",
		  cell_2d_sd = "vector",
		  blist_genes = "vector")
)

#' @export
setMethod(
	"initialize",
	signature = "tgScMatClustLayout2D",
	definition =
	  function(.Object, scl=NULL) {
	   if(!is.null(scl)) {
		   .Object@scl = scl
		   .Object@k_confusion = -1
	   }
	   return(.Object)
	 }
	)

#' Init a single cell plotting object
#'
#' Init a single cell plotting object, tgScMatClust2D
#'
#' @param scc A tgScMatClust object, created using scc_cluster().
#'
#' @return A tgScMatClustLayout2D object, used in all of the plotting functions.
#'
#' @export
scp_init_plot = function(scc) {
	return(tgScMatClustLayout2D(scc))
}

#' Compute cells order either by confusion matrix or by selected marker genes
#'
#' @param scm2d 
#' @param K 
#' @param use_confu 
#' @param per_clust_genes 
#' @param gene_min_fold 
#' @param gene_min_cov 
#' @param ignore_blist_genes 
#'
#' @return
#' @export
#'
#' @examples
scp_compute_cell_order = function(scm2d, K=-1, use_confu=F, 
																	per_clust_genes= get_param("scp_ord_per_clust_genes"),
																	gene_min_fold = get_param("scp_ord_gene_min_fold"), 
																	gene_min_cov = get_param("scp_ord_gene_min_cov"), 
																	ignore_blist_genes = get_param("scp_ord_ignore_blist_genes"))
{
	if(length(scm2d@cl_ord) > 1) {
		message("using cached cell order")
		return(scm2d)
	}
	if(use_confu) {
		scm2d = scp_compute_cell_order_confu(scm2d, K)
	} else {
		if (ignore_blist_genes) {
			ignore_genes = scm2d@blist_genes
		}
		else {
			ignore_genes = NULL
		}
		scm2d = scp_compute_cell_order_markers(scm2d, per_clust_genes=per_clust_genes, gene_min_fold=gene_min_fold, gene_min_cov = gene_min_cov, ignore_genes = ignore_genes)
	}
	return(scm2d)
}

#' Order cells by confusion matrix
#'
#' @param scm2d 
#' @param K 
#'
#' @return
#' @export
#'
#' @examples
scp_compute_cell_order_confu = function(scm2d, K=-1)
{
	message("compute cell order using confu")
	if(K == -1) {
		scm2d = .scp_compute_k_confusion(scm2d)
		K = scm2d@k_confusion
	}
#currently we support ordering using correlation on the confusion matrix
	scm2d = scp_compute_balanced_clust_confusion(scm2d, 
					K, 
					k_expand_inout_factor=10)

	confu = scm2d@mat_bal_confusion
	# c_confu = cor(confu)
  	hc = hclust(as.dist(log(1/(confu+1e-4))),"ward.D2")
	scm2d@cl_ord = hc$order
#will be nice to hclust within each cluster given confusion matrix
	new_rank = rep(NA, length(hc$order))
	new_rank[hc$order] = (1:length(hc$order))
	scm2d@c_ord = order(new_rank[scm2d@scl@clusts])
	return(scm2d)
}

#' Order cells by hclust on marker genes
#'
#' @param scm2d 
#' @param per_clust_genes 
#' @param gene_min_fold 
#' @param gene_min_cov 
#' @param ignore_genes 
#'
#' @return
#' @export
#'
#' @examples
scp_compute_cell_order_markers = function(scm2d,
		per_clust_genes=5,
		gene_min_fold = 2.5, gene_min_cov = 0.5, ignore_genes = NULL)
{
	message("compute cell order using hc on marks")
#'
	scl = scm2d@scl
	gene_folds = scm2d@scl@clust_fp

	genes_pool = rownames(gene_folds)
	unmask_folds = gene_folds[genes_pool,]
	mask_folds = gene_folds[genes_pool,]*(scl@clust_gcov[genes_pool,]>gene_min_cov)
	mask_folds[mask_folds == 0] = 1

	good_marks = unique(as.vector(unlist(
			apply(mask_folds,
				2,
				function(x)  {
				   names(head(sort(-x[x>gene_min_fold]),n=per_clust_genes)) })
		     )))

	good_marks = setdiff(good_marks, ignore_genes)
	
	if(is.null(good_marks) | length(good_marks) < 2) {
		return(scp_compute_cell_order_confu(scm2d))
	}
	feat = log2(unmask_folds[good_marks,])

	hc = hclust(dist(cor(feat)), "ward.D2")

	g_ncover = apply(feat > 1, 1, sum)
	
	main_mark = names(g_ncover)[which.max(g_ncover)]
	f = feat[main_mark,] < 0.25
	g_ncover = apply(feat[,f]>1, 1, sum)
	second_mark = names(g_ncover)[which.max(g_ncover)]
	message("reorder on ", main_mark, " vs ", second_mark)
	d = reorder(as.dendrogram(hc), 
				feat[main_mark,]-feat[second_mark,], 
				agglo.FUN=mean)
	hc2 = as.hclust(d)
	scm2d@cl_ord = hc2$order
#will be nice to hclust within each cluster given confusion matrix
	new_rank = rep(NA, length(hc2$order))
	new_rank[hc2$order] = (1:length(hc2$order))
	scm2d@c_ord = order(new_rank[scm2d@scl@clusts])
	return(scm2d)
}
	
.scp_compute_k_confusion = function(scm2d)
{
	if(scm2d@k_confusion == -1) {
		#think about better heuristics here
		scm2d@k_confusion = 20
	}
	return(scm2d)
}

#' @export
scp_compute_balanced_clust_confusion = function(scm2d, K = -1, k_expand_inout_factor=10)
{
	if(length(scm2d@mat_bal_confusion)>1 & (K == -1 | K == scm2d@k_confusion)) {	#exists
		return(scm2d)
	}
	message("computing balanced confusion matrix")
	if(K == -1) {
		scm2d = .scp_compute_k_confusion(scm2d)
		K = scm2d@k_confusion
	}

	scl= scm2d@scl
	scl = .scc_comp_balanced_knn_matrix(scl, K, k_expand_inout_factor)
	message("gen full knn matrix")
	mknnio_ij = as.matrix(scl@m_knn)
#column j contain the K-nn (somehow balanced) neighbors for cell j
#row i may have anything between 0 to a large number of cells that consider i as a neighbor

	message("gen clust _ci")
	isclust_ci = diag(max(scl@clusts))[,scl@clusts]

	message("gen confu cc")
	confu_cc =isclust_ci %*% (mknnio_ij>0) %*% t(isclust_ci)
	message("gen confu cc")

	rownames(confu_cc) = 1:nrow(confu_cc)
	colnames(confu_cc) = 1:nrow(confu_cc)

	scm2d@mat_bal_confusion = confu_cc
	scm2d@mknnio_ij = mknnio_ij
	scm2d@k_confusion = as.integer(K)

	message("out of computing balanced confu")
	return(scm2d)
}

#' @export
scp_compute_confusion = function(scm2d, K = -1, sim_mat=NULL)
{
	if(length(scm2d@mat_k_confusion)>1 & (K == -1 | K == scm2d@k_confusion)) {	#exists
		return(scm2d)
	}
	message("computing confusion matrix")
	if(K == -1) {
		scm2d = .scp_compute_k_confusion(scm2d)
		K = scm2d@k_confusion
	}
	scm_cl = scm2d@scl
	if(is.null(sim_mat)) {
		sim_mat = as.matrix(cor(as.matrix(scm_cl@feat_mat)))
	}
	clust = scm_cl@clusts

	confu = .knn_on_clusts(sim_mat, clust,
					rownames(sim_mat),
					rownames(sim_mat), K)
	# confu_self = .knn_on_clusts(sim_mat, clust,
	# 			rownames(sim_mat),
	# 			rownames(sim_mat), round(max(K/3,5)))

	diag(confu) = diag(confu) - tabulate(clust)
	# diag(confu_self) = diag(confu_self) - tabulate(clust)

	# print(self_score)
	rownames(confu) = 1:nrow(confu)
	colnames(confu) = 1:ncol(confu)
	scm2d@mat_k_confusion = confu
	# scm2d@mat_k_confusion_self = confu_self
	scm2d@k_confusion = as.integer(K)
	return(scm2d)
}

#' @export
scp_set_clust_cols = function(scm2d, clust_cols=NULL,
		marker_colors = NULL, Tmark=2, Tu=5, cell_color_by_marks = F)
{
	if(!is.null(clust_cols)) {
		scm2d@cl_colors = clust_cols
		scm2d@c_colors = clust_cols[scm2d@scl@clusts]
	} else if(!is.null(marker_colors)) {
		scm2d@cl_colors = rep("white", ncol(scm2d@scl@clust_fp))
		good_marks = intersect(rownames(marker_colors),
						rownames(scm2d@scl@clust_fp))
		good_colors = marker_colors[good_marks, ]
		
		if(length(good_marks) > 1) {
			scm2d@mark_colors = as.character(good_colors$color)
			names(scm2d@mark_colors) = good_marks
			
			scm2d@mark_groups = as.character(good_colors$group)
			names(scm2d@mark_groups) = as.character(good_colors$color)

			cl_colors = rep(NA, length(scm2d@cl_colors))
			c_colors = rep(NA, length(scm2d@scl@scmat@ncells))
			
			for (p in 1:max(good_colors$priority)) {
				mark_fp = scm2d@scl@clust_fp[good_marks[good_colors$priority == p],]
				mark_fp = rbind(mark_fp, mark_fp) # quick and dirty - make sure mark_fp has > 1 row
				
				f = is.na(cl_colors) & apply(mark_fp, 2, max) > Tmark
				cl_colors[f] = as.character(good_colors[good_colors$priority == p, 'color'])[apply(t(mark_fp[, f]), 1, which.max)]

				mark_u = scm2d@scl@scmat@mat[good_marks[good_colors$priority == p],]
				mark_u = rbind(mark_u, mark_u)
				f = is.na(c_colors) & apply(mark_u, 2, max) > Tu
				c_colors[f] = as.character(good_colors[good_colors$priority == p, 'color'])[apply(t(mark_u[, f]), 1, which.max)]
			}
			cl_colors[is.na(cl_colors)] = "white"
			scm2d@cl_colors = cl_colors
			
			if(cell_color_by_marks) {
				c_colors[is.na(c_colors)] = "white"
				scm2d@c_colors = c_colors
			}
			else {
				scm2d@c_colors = scm2d@cl_colors[scm2d@scl@clusts]
			}
		} else if(length(good_marks) == 1) {
			f = scm2d@scl@clust_fp[good_marks[1],]>Tmark
			scm2d@cl_colors[f] = as.character(marker_colors[good_marks[1], 'color'])
			scm2d@c_colors = scm2d@cl_colors[scm2d@scl@clusts]
		}
	} else {
		clust_cols = colorRampPalette(c("white", "lightgray", "darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue", " cyan"))(max(scm2d@scl@clusts))
		scm2d = scp_compute_cell_order(scm2d)
		cl_reord = rep(NA, length(scm2d@cl_ord))
		cl_reord[scm2d@cl_ord] = (1:length(scm2d@cl_ord))
		scm2d@cl_colors = clust_cols[cl_reord]
		scm2d@c_colors = clust_cols[cl_reord[scm2d@scl@clusts]]
	}
	return(scm2d)
}

# m - similarity mat
# clusts - cell assignment to clusters
# v1, v2 - list of cells
#
# returns for each cluster i, how many of the k nearst neighbors of it's cells come from each cluster.
.knn_on_clusts = function(m, clusts, v1, v2, K)
{
	a = apply(m[v1,v2], 1, function(x) names(tail(sort(x),n=K)))
	confu = table(clusts[as.character(unlist(a))],
			rep(clusts[colnames(a)],each=K))
	return(confu)
}

# construct a 2D cluster layout object using 2D graph projection of the confusion grap
#'
#' @importClassesFrom graph graphNEL
#' @importFrom graph plot addEdge addNode nodeRenderInfo
#' @importFrom Rgraphviz layoutGraph
#' @export
scp_compute_clust_knn_graph= function(scm2d, K=-1, T_edge=0.05,
				force_max_deg = 4, force_one_compo = F,
				T_edge_asym = F,
				k_expand_inout_factor=10,
				subselect_by_fpcor=T
				)
{
	restrict_in_degree = T

	scm_cl = scm2d@scl
	clust = scm_cl@clusts

	message("comp clust knn graph = will gen balanced confusion")
	scm2d = scp_compute_balanced_clust_confusion(scm2d, K, k_expand_inout_factor=k_expand_inout_factor)
	message("done gen balanced confusion")

	K = scm2d@k_confusion

	confu = scm2d@mat_bal_confusion

	csize = as.matrix(table(scm_cl@clusts))
	csize = pmax(csize, 20)
	csize2 = csize %*% t(csize)
	csize2 = csize2 / median(csize)**2
	confu = confu / csize2

	confu_p_from = confu/rowSums(confu)
	confu_p_to = t(confu)/colSums(confu)
	if(!is.na(force_max_deg)) {
		rank_fr = t(apply(confu_p_from, 1, rank))
		rank_to = t(apply(confu_p_to, 1, rank))
		rank2 = rank_fr * rank_to
		diag(rank2) = 1e+6
		amgraph = apply(rank2, 1, function(x) {  rank(-x) <= (1+force_max_deg) })
		mgraph = amgraph * ((confu_p_from + confu_p_to)>T_edge)
		if(restrict_in_degree) {
			amgraph2 = t(apply(rank2, 2, function(x) {  rank(-x) <= (1+force_max_deg) }))
			mgraph = mgraph * amgraph2
		}

		if(T_edge_asym) {
			mgraph = amgraph * (confu_p_from>T_edge)
			mgraph = amgraph * (t(confu_p_to)>T_edge)
		}
		mgraph = mgraph>0 | t(mgraph>0)
	} else {
		mgraph = (confu_p_from + confu_p_to) > T_edge
	}

	if(force_one_compo) {
		stop("connecting components is not suppurted anymore")
		d = as.matrix(scm_cl@feat_mat)
		sim_mat = as.matrix(cor(d))
		mgraph = .force_connected_components(mgraph, clust, K, sim_mat, add_multi_comp_edge=F)
	}
	if(subselect_by_fpcor & !is.na(force_max_deg)) {
		fp_cor = cor(log2(scm2d@scl@clust_fp))
		fp_rnk = t(apply(-fp_cor,1,function(x) rank(x)<3*force_max_deg))
		mgraph = mgraph * fp_rnk * t(fp_rnk)
	}
	scm2d@clust_graph = mgraph
	return(scm2d)
}
scp_project_clust_subgraph = function(sc2d, foc_clusts=NULL,
			T_edge = 0.05, T_edge_asym = F,
			force_max_deg = 4,
			add_neigh = T, blur=0.02, K_for_cells=NA)
{
	oldseed = .set_seed()
	if(length(sc2d@clust_graph) == 0) {
		sc2d = scp_compute_clust_knn_graph(sc2d,
				T_edge = T_edge, T_edge_asym = T_edge_asym,
				force_max_deg = force_max_deg)
	}
	if(is.na(K_for_cells)) {
		K_for_cells = round(sc2d@k_confusion)
	}
	mgraph = sc2d@clust_graph
	if(!is.null(foc_clusts)) {
		if(add_neigh) {
			nclusts = rowSums(mgraph[,foc_clusts])>0
			message("foc clust sz ", length(foc_clusts), " neighbors in total ", sum(nclusts), " clusts")
			mgraph = mgraph[nclusts, nclusts]
		} else {
			mgraph = mgraph[foc_clusts, foc_clusts]
		}
	}

	g2d = .scp_comp_clust_coord(mgraph)
	x_cl = g2d$x_cl
	y_cl = g2d$y_cl
	sc2d@x_cl = x_cl
	sc2d@y_cl = y_cl

	dbg_g2d <<- g2d

	blurx = blur*(max(x_cl) - min(x_cl))
	blury = blur*(max(y_cl) - min(y_cl))

	clust = as.character(sc2d@scl@clusts)

	cell_neigh_ids = t(apply(sc2d@mknnio_ij, 2, function(x) {
					names(x) = 1:length(x);
					neigh = x[x>0]
					neigh = as.numeric(names(head(sort(-neigh),n=2*K_for_cells)))
					if(length(neigh) < 2*K_for_cells) {
						neigh = c(neigh, rep(NA, 2*K_for_cells-length(neigh)))
					}
					return(neigh)
				 }))

	na_x = -1
	na_y = -1

	px0 = apply(cbind(1:nrow(cell_neigh_ids),cell_neigh_ids), 1,
			function(ids) {
				foc_clid = clust[ids[1]]
				cl_connect = mgraph[foc_clid,][clust[ids]]
				xs=x_cl[clust[ids]];
				xs=xs[!is.na(xs) & cl_connect];
				if(length(xs)>0) {
					return(mean(head(xs, K_for_cells)));
				} else {
					return(na_x);
				}
			})
	py0 = apply(cbind(1:nrow(cell_neigh_ids),cell_neigh_ids), 1,
			function(ids) {
				foc_clid = clust[ids[1]]
				cl_connect = mgraph[foc_clid,][clust[ids]]
				ys=y_cl[clust[ids]];
				ys=ys[!is.na(ys) & cl_connect];
				if(length(ys)>0) {
					return(mean(head(ys, K_for_cells)));
				} else {
					return(na_y);
				}
			})

	if(is.null(foc_clusts) & sum(px0 == na_x) > 0) {
		message("when projecting the entire graph, ", sum(px0 == na_x), " elemets are ill-positioned (will be shown on lower left)")
	}
	if(is.null(foc_clusts)) {
		base_x = min(px0[px0 != na_x],na.rm=T)
		base_y = min(py0[py0 != na_y],na.rm=T)
		max_x = max(px0[px0 != na_x],na.rm=T)
		base_x = base_x - (max_x-base_x)*0.1

		px0[px0 == na_x] = base_x
		py0[py0 == na_y] = base_y
	} else {
		px0[px0 == na_x] = NA
		py0[py0 == na_y] = NA
	}
	px = px0 + rnorm(mean=0, sd=blurx, n=length(px0))
	py = py0 + rnorm(mean=0, sd=blury, n=length(py0))

	sc2d@x = px
	sc2d@y = py

	proj_r2 = function(x) {
		ids = x[c(-1,-2)]
		cent_x = x[1]
		cent_y = x[2]
		offs_x = x_cl[clust[ids]]-cent_x
		offs_y = y_cl[clust[ids]]-cent_y
		x2 = offs_x * offs_x
		y2 = offs_y * offs_y
		return(sqrt(mean(x2+y2,na.rm=T)))
	}

	coord_v = apply(cbind(px0, py0, cell_neigh_ids), 1, proj_r2)
#	names(coord_v) = names(clust)
#	if (plot_graph) {
#	  plot(g, nodeAttrs = list(),"neato")
#	}

	.restore_seed(oldseed)

	sc2d@cell_2d_sd = coord_v
	return(sc2d)
}

.scp_comp_clust_coord = function(mgraph)
{
	N = nrow(mgraph)
	rEG <- new("graphNEL", nodes=as.character(1:N), edgemode="undirected")

	e = which(mgraph>0)
	n1 = ceiling((e)/N)
	n2 = 1+((e-1) %% N)

	rEG = addEdge(as.character(n1[n1!=n2]), as.character(n2[n1!=n2]), rEG, rep(1, length(n1[n1!=n2])))

	g = layoutGraph(rEG, layoutType="neato")
	x_cl = nodeRenderInfo(g)$nodeX
	y_cl = nodeRenderInfo(g)$nodeY
	names(x_cl) = rownames(mgraph)
	names(y_cl) = rownames(mgraph)
	return(list(g=g, x_cl=x_cl, y_cl=y_cl))
}

# add edges to connect all components of a graph
.force_connected_components <- function(mgraph, clust, K, sim_mat, add_multi_comp_edge) {
  compo = .graph_con_comp(mgraph)
  while(max(compo) > 1) {
    message("try to merge components of sizes ");
    print(table(compo))
    # print(sort(compo))
    compo_cells = compo[clust]
    names(compo_cells) = names(clust)
    #add bridges to each con comp
    compsize = table(compo)
    big_comp_i = which.max(compsize)
    for(comp_i in 1:max(compo)) {
      if(comp_i == big_comp_i) {
        next
      }
      comp_cells = which(compo_cells == comp_i)
      noncomp_cells = which(compo_cells != comp_i)
      tot_neigh = min(floor(length(comp_cells)*K/2), 2*K)
      #find k closest pairs
      c_neigh = list()
      c_neigh$nn.dists = t(apply(
        sim_mat[noncomp_cells, comp_cells],
        1,
        function(x) -tail(sort(x),n=round(K/2))))
      c_neigh$nn.idx = t(apply(
        sim_mat[noncomp_cells, comp_cells],
        1,
        function(x) {
          names(x) = 1:length(x);
          return(as.numeric(names(tail(sort(x),n=round(K/2)))))
        }
      ))
      #quick way to get the K'th best distance
      thresh = sort(c_neigh$nn.dists, partial=tot_neigh)[tot_neigh]
      #extracting the id's of the closest noncomp cells
      #we may have duplicates here, but it is OK
      neigh_cells = c_neigh$nn.idx[c_neigh$nn.dists <= thresh]
      neigh_cells_nm = names(noncomp_cells)[neigh_cells]
      #the clusts of the best comp neighbor cells
      src_cells = which(c_neigh$nn.dists <= thresh) %% length(comp_cells)
      src_cells_nm = names(comp_cells)[src_cells]

      neigh_clusts = clust[neigh_cells_nm]
      src_clusts = clust[src_cells_nm]

      potent_edge = sort(table(paste(neigh_clusts, src_clusts)))
      potent_edge = potent_edge/sum(potent_edge)
      edges = c(names(potent_edge)[which.max(potent_edge)])
      if(add_multi_comp_edge) {
        edges = c(edges, names(which(potent_edge > 0.1)))
      }
      for(e in edges) {
        n = unlist(strsplit(e, " "))
        message("connect comp, ", comp_i, " nodes ", n[1], " ", n[2])
        message("previous connected? ", mgraph[as.numeric(n[1]), as.numeric(n[2])])

        mgraph[as.numeric(n[1]), as.numeric(n[2])] = TRUE
        mgraph[as.numeric(n[2]), as.numeric(n[1])] = TRUE
      }
      #connect the top one, and all those
    }
    message("will recompute connect comp, total edges ", sum(mgraph))
    compo = .graph_con_comp(mgraph)
  }
  return(mgraph)
}

#' Plotting a the cluster confusion matrix
#'
#' @param scm2d an object of clust layout
#' @param outdir - file location
#' @param K (20) - the parameter
#'
#' @export
scp_plot_confusion_mat = function(scm2d, outdir=".",
					K=-1,
					width=1000, height=1000)
# defualt K has to be -1, o/w confusion matrix is re-calculated.
{
	scm2d = scp_compute_balanced_clust_confusion(scm2d, K, 
						k_expand_inout_factor=10)

	confu = scm2d@mat_bal_confusion
	K = scm2d@k_confusion
	scm2d = scp_compute_cell_order(scm2d, K)

	cl_ord = scm2d@cl_ord
  shades = colorRampPalette(c("white", "lightblue", "blue", "red", "yellow", "black"))(1000);
  png(sprintf("%s/confusion_k%d.png", outdir, K), w=width, h=height+50)
	layout(matrix(c(1,2), nrow=2), heights=c(height, 100))
	par(mar = c(0,8,5,8))
	image(log2(1+confu[cl_ord, cl_ord]), col=shades, xaxt='n', yaxt='n')
	mtext(rownames(confu)[cl_ord], at=seq(0,1,l=length(cl_ord)), side=2, las=2)
	mtext(rownames(confu)[cl_ord], at=seq(0,1,l=length(cl_ord)), side=3, las=2)
	mtext(rownames(confu)[cl_ord], at=seq(0,1,l=length(cl_ord)), side=4, las=2)

	par(mar = c(5,8,0,8))
	image(as.matrix(1:length(scm2d@cl_ord),nrow=1), col=scm2d@cl_colors[scm2d@cl_ord], yaxt='n', xaxt='n')
	mtext(rownames(confu)[cl_ord], at=seq(0,1,l=length(cl_ord)), side=1, las=2)
  	dev.off()

	c_confu = cor(log2(1/(confu+0.001)))
	diag(c_confu) = 0
  	png(sprintf("%s/confusion_cor_k%d.png", outdir, K), w=width,h=height+50)
	layout(matrix(c(1,2), nrow=2), heights=c(height, 100))
	par(mar = c(0,8,5,8))
  	image(c_confu[cl_ord, cl_ord], col=shades, xaxt='n', yaxt='n')
	mtext(rownames(confu)[cl_ord], at=seq(0,1,l=length(cl_ord)), side=2, las=2)
	mtext(rownames(confu)[cl_ord], at=seq(0,1,l=length(cl_ord)), side=3, las=2)
	mtext(rownames(confu)[cl_ord], at=seq(0,1,l=length(cl_ord)), side=4, las=2)

	par(mar = c(5,8,0,8))
	image(as.matrix(1:length(scm2d@cl_ord),nrow=1), col=scm2d@cl_colors[scm2d@cl_ord], yaxt='n', xaxt='n')
	mtext(rownames(confu)[cl_ord], at=seq(0,1,l=length(cl_ord)), side=1, las=2)
  	dev.off()
	return(scm2d)
}

#' Plotting a heatmap of clusters
#'
#'
#' @param scl2d
#' @param png_fn 
#' @param plot_cells 
#' @param genes_pool 
#' @param blacklist_genes 
#' @param width 
#' @param height 
#' @param text_cex 
#' @param alt_side_text 
#' @param per_clust_genes 
#' @param gene_min_fold 
#' @param gene_min_cov 
#' @param norm_type 
#' @param tot_umi_scale 
#' @param smooth_n 
#' @param top_marg 
#' @param lower_marg 
#' @param fp_shades 
#' @param cell_batch 
#' @param batch_color 
#' @param cell_ord 
#' @param genes 
#'
#' @export
#' 
#' 
#' 
scp_plot_mat = function(scl2d, png_fn=NULL,
		plot_cells = T,
		genes_pool = NULL,
		blacklist_genes = c(),
		remove_blacklist = T,
		width=2000, height=NA,
		text_cex = 1, alt_side_text=F,
		per_clust_genes=5,
		gene_min_fold = 2.5, gene_min_cov = 0.5,
		norm_type="linear", tot_umi_scale = 1000,
		smooth_n = 5,
		top_marg=c(0,13,5,20),
		lower_marg=c(5,13,0,20),
		fp_shades = colorRampPalette(c("white", "lightsteelblue", "blue", "darkblue", "purple", "brown", "black","yellow"))(1000),
		cell_batch = NULL,
		batch_color = NULL,
		cell_ord = NULL, 
		genes = NULL)
{
	scl = scl2d@scl

	if(is.null(cell_ord)) {
	  scl2d = scp_compute_cell_order(scl2d)
	  cell_ord = scl2d@c_ord
	  cl_ord = scl2d@cl_ord
	} else {
	  #need to test this!!!
	  cl_ord = order(tapply(1:length(cell_ord), scl@clusts[cell_ord], mean))
	}

	if (is.null(genes)) {
		good_marks = scc_get_ordered_marks(
			scl,
			cl_ord,
			genes_pool = genes_pool,
			blacklist_genes = blacklist_genes,
			remove_blist_genes = remove_blacklist,
			per_clust_genes = per_clust_genes,
			gene_min_fold = gene_min_fold,
			gene_min_cov = gene_min_cov
		)
	}
	else {
		good_marks = intersect(genes, rownames(scl2d@scl@clust_fp))
		if (!is.null(genes_pool)) {
			good_marks = intersect(good_marks, genes_pool)
		}
	}

	if(is.null(good_marks)) {
		if(!is.null(png_fn)) {
			png(png_fn)
			dev.off()
		}
		message("cannot find good marks to plot gene matrix")
		return(scl2d)
	}

	gene_folds = scl2d@scl@clust_fp
	clusts = scl@clusts

	if(!is.null(png_fn)) {
	  if(is.na(height)){
	    height = 16*length(good_marks) + 100
	  }
	  png(png_fn, w=width,h=height);
	}
	layout(matrix(c(1,2),nrow=2),heights=c(height, 100))
	par(mar=top_marg)

	if(plot_cells) {
		mat = scl2d@scl@scmat@mat[good_marks, names(clusts)]
		if(norm_type=="linear") {
			totu = colSums(scl2d@scl@scmat@mat)
			mat = t(t(mat)/totu)*tot_umi_scale
		} else {
			message("norm type ", norm_type, " is not supported")
		}

		lus_1 = log2(1+7*mat[good_marks, names(clusts)])
		lus = apply(lus_1 - apply(lus_1, 1, median),2, function(x) pmax(x,0))
		if (length(cell_ord) < width) {
			smooth_n = 1
		}
		lus_smoo = t(apply(lus[,cell_ord], 1, function(x) rollmean(x,smooth_n, fill=0)))
		image(t(lus_smoo), col=fp_shades, xaxt='n', yaxt='n')
		cell_x = rep(NA, length(cell_ord))
		cell_x[cell_ord] = 1:length(cell_ord)
		cl_x = tapply(cell_x, clusts, mean)/length(cell_ord)
		mtext(1:scl@nclust, side = 3, at=cl_x, las=2, line = 2, cex=text_cex)
		cl_x_b = tapply(cell_x, clusts, max)/length(cell_ord)
		abline(v=cl_x_b, lwd=2)
	} else {
		mat = log2(gene_folds[good_marks, cl_ord])
		mat = pmax(pmin(mat,3),-3)
		image(t(mat), col=fp_shades, xaxt='n', yaxt='n', zlim=c(-3,3))
		n = length(scl2d@cl_colors)
		mtext(scl2d@cl_ord, side = 3, at=seq(0,1,l=scl@nclust), las=2, line = 2, cex=text_cex)
	}

	g_n = length(good_marks)
	
	gene_cols = rep("black", g_n)
	gene_cols[is.element(good_marks, blacklist_genes)] = "red"
	
	if(alt_side_text) {
		odd = seq(1,g_n,2)
		even = seq(2,g_n,2)
		mtext(substr(good_marks[odd],1,8),
			at=seq(0,1,length.out=g_n)[odd],
			side=2, las=2, cex=text_cex, col=gene_cols[odd])
		mtext(substr(good_marks[even],1,8),
			at=seq(0,1,length.out=g_n)[even],
			side=4, las=2, cex=text_cex, col=gene_cols[even])
	} else {
		mtext(substr(good_marks,1,8),
			at=seq(0,1,length.out=g_n),
			side=2, las=2, cex=text_cex, col=gene_cols)
		mtext(substr(good_marks,1,8),
			at=seq(0,1,length.out=g_n),
			side=4, las=2, cex=text_cex, col=gene_cols)
	}

		par(mar=lower_marg)
	if(plot_cells) {
		image(as.matrix(1:length(cell_ord),nrow=1), col=scl2d@c_colors[cell_ord], xaxt='n', yaxt='n')
		mtext(1:scl@nclust, side = 1, at=cl_x, las=2, line = 2, cex=text_cex)
	} else {
		image(as.matrix(1:length(scl2d@cl_ord),nrow=1), col=scl2d@cl_colors[scl2d@cl_ord], yaxt='n', xaxt='n')
		mtext(scl2d@cl_ord, side = 1, at=seq(0,1,l=scl@nclust), las=2, line = 2, cex=text_cex)
	}
	if(!is.null(png_fn)) {
		write.table(x=as.matrix(good_marks), file=sprintf("%s.marks_labels", png_fn), quote=F)
	}

	if(!is.null(png_fn)) {
		dev.off()
	}

	if (!plot_cells) {
		marks_per_plot = 81
		for (i in seq(1, length(good_marks), by=marks_per_plot)) {
			marks = good_marks[seq(i, min(i + marks_per_plot - 1, length(good_marks)))]
			nx = floor(sqrt(marks_per_plot))
			ny = ceiling(marks_per_plot / nx)
			png(sprintf("%s_barplots_%d.png", sub(".png", "", png_fn), ceiling(i/marks_per_plot)), width=nx * 240, height=ny * 200, res=250, pointsize=3)
			layout(matrix(1:marks_per_plot, ny, nx, byrow=T))
			par(mar=c(3,3,4,1))
			lapply(marks, function(m) { barplot(log2(scl@clust_fp[m, scl2d@cl_ord]), border=NA, cex.names=1.5, cex.axis=1.5, las=2, col=scl2d@cl_colors[scl2d@cl_ord]); title(main=m, cex.main=3) })
			dev.off()
		}
	}
		
	
	return(scl2d)
}

#	if(!is.null(cell_batch)) {
#		png(file=sprintf("%s.batch_key.png", png_fn), w=width,h=200)
#		par(mar=c(0,0,0,0))
#		if(is.null(batch_color)) {
#		  batch_color = colorRampPalette(c("darkgray", "burlywood1","chocolate4","orange", "red",
#		                                  "purple", "blue", " cyan"))(max(cell_batch))
#			# batch_color = colorRampPalette(c("gray","red", "blue", "yellow"))(max(cell_batch))
#		}
#		nms_in_ord = colnames(mat)[cell_ord]
#		image(as.matrix(cell_batch[nms_in_ord]), col=batch_color)
#		dev.off()
#	}
#' @export
#'
scp_plot_metadata_factor = function(scl2d,
			clust_breakdown=F,
			meta_field_nm,
			height=500, width=1000,
			cell_ord=NULL,
			meta_cols = NULL,
			fname=NA,
			heatmap = F,
			heatmap.cex = 1.8)
{
	if(is.null(cell_ord)) {
		scl2d = scp_compute_cell_order(scl2d)
		cell_ord = scl2d@c_ord
	}
	#we can go back to use pheat map or anything, or just add a legend
	fact = factor(scl2d@scl@scmat@cell_metadata[cell_ord, meta_field_nm])
	meta = as.integer(fact)
	meta_names = levels(fact)
	if(is.null(meta_cols)) {
		meta_cols = colorRampPalette(c("gray","red", "blue", "yellow"))(max(meta))
	}
	png(fname, w=width, h= height)
	if(clust_breakdown) {
		bdown = table(meta, scl2d@scl@clusts[cell_ord])
		bdown = t(t(bdown)/colSums(bdown))
		if (heatmap) {
			par(mar = c(0,0,0,20))
			image(t(bdown[nrow(bdown):1,scl2d@cl_ord]),
						col = colorRampPalette(c("white", "navy", "black"))(100), breaks = seq(0, 1, l=101))
			par(las = 1)
			mtext(rev(meta_names),
						side = 4, at= (0:(nrow(bdown)-1))/(nrow(bdown)-1), cex = heatmap.cex)
		}else {
			barplot(bdown[,scl2d@cl_ord], col=meta_cols)
		}
	} else {
		image(as.matrix(meta), col=meta_cols)
	}
	# write.table()
	dev.off()
	return(scl2d)
}


#' @export
#'
#'
scp_plot_clust_2d = function(scl2d, fname = NULL, foc_clusts = NULL,
				T_edge = 0.05, T_edge_asym=F,
				add_neigh = T, blur=0.02, K_for_cells=NA,
				force_max_deg = force_max_deg,
				outcol="black", height=2000, width=2000, cex=1.5, expand_by = 50,
				plot_cells = T, plot_clusts = T, plot_edges = F, plot_markers = F, cell_cols = NULL)
{
	if(length(scl2d@cl_colors) == 0) {
		scl2d = scp_set_clust_cols(scl2d)
	}
	clust_cols = scl2d@cl_colors
	if (is.null(cell_cols)) {
		cell_cols = scl2d@c_colors
	}

	if(length(scl2d@clust_graph) == 0) {
		scl2d = scp_compute_clust_knn_graph(scl2d,
				T_edge = T_edge, T_edge_asym = T_edge_asym,
				force_max_deg = force_max_deg)
	}
	message("will project graph to 2d due to popular demand, remember its the wrong thing to do")
	scl2d = scp_project_clust_subgraph(scl2d, foc_clusts = foc_clusts,
					add_neigh = T, blur=blur, K_for_cells=NA)

	if(!is.null(fname)) {
		png(fname, h=height, w=width);
	}
	tp = 'p'
	if(plot_cells == F) {
		tp = 'n'
	}
	xmin = min(scl2d@x, scl2d@x_cl) - expand_by
	xmax = max(scl2d@x, scl2d@x_cl) + expand_by

	ymin = min(scl2d@y, scl2d@y_cl) - expand_by
	ymax = max(scl2d@y, scl2d@y_cl) + expand_by * ifelse(plot_markers, 3, 1)

	plot(scl2d@x, scl2d@y, pch=21, bg=cell_cols,
		col=ifelse(is.na(outcol), cell_cols, outcol), lwd=0.5, 
		cex=cex, type=tp,
		xlim=c(xmin,xmax),
		ylim=c(ymin,ymax));

	if(plot_edges) {
		e = which(scl2d@clust_graph>0)
		N = nrow(scl2d@clust_graph)
		n1 = ceiling((e)/N)
		n2 = 1+((e-1) %% N)
		xcl = scl2d@x_cl
		ycl = scl2d@y_cl
		segments(xcl[n1], ycl[n1], xcl[n2], ycl[n2], lwd=1)
	}
	# points(scl2d@x_cl, scl2d@y_cl, pch=21, bg =clust_cols, cex=cex*3)
	if (plot_clusts) {
		points(scl2d@x_cl, scl2d@y_cl, pch=21, bg =clust_cols[as.numeric(names(scl2d@x_cl))], cex=cex*3)
		text(scl2d@x_cl, scl2d@y_cl, names(scl2d@y_cl), cex=1.5)
	}

	if (plot_markers) {
		gmark = tapply(names(scl2d@mark_colors), scl2d@mark_groups, paste, collapse=", ")
		gcol = unique(data.frame(col=scl2d@mark_colors, group=scl2d@mark_groups))
		rownames(gcol) = gcol$group
		gmark = gmark[order(names(gmark))]
		legend("topleft", legend=gsub("_", " ", paste0(names(gmark), ": ", gmark)), pch=19, cex=2, col=as.character(gcol[names(gmark), 'col']), bty='n')
	}
	
#	grid(lwd=2,col="black",lty=1)
	if(!is.null(fname)) {
		dev.off()
	}
	return(scl2d)
}

#' @export
#'
#'
scp_plot_clust_2d_by_meta_field = function(scl2d, meta_field, odir, foc_clusts = NULL,
														 T_edge = 0.05, T_edge_asym=F,
														 add_neigh = T, blur=0.02, K_for_cells=NA,
														 force_max_deg = force_max_deg,
														 bg_col="darkgrey", fg_col="red", single_plot = T, panel_size = 200, cex=1.5, expand_by = 50)
{
	
	if(length(scl2d@clust_graph) == 0) {
		scl2d = scp_compute_clust_knn_graph(scl2d,
																				T_edge = T_edge, T_edge_asym = T_edge_asym,
																				force_max_deg = force_max_deg)
	}
	scl2d = scp_project_clust_subgraph(scl2d, foc_clusts = foc_clusts,
																		 add_neigh = T, blur=blur, K_for_cells=NA)
	
	xmin = min(scl2d@x, scl2d@x_cl) - expand_by
	xmax = max(scl2d@x, scl2d@x_cl) + expand_by
	
	ymin = min(scl2d@y, scl2d@y_cl) - expand_by
	ymax = max(scl2d@y, scl2d@y_cl) + expand_by
	
	c_by_f = split(scl2d@scl@scmat@cells, scl2d@scl@scmat@cell_metadata[, meta_field])
	
	dir.create(odir, showWarnings = F, recursive = T)
	
	if (single_plot) {
		ny = floor(sqrt(length(c_by_f)))
		nx = ceiling((length(c_by_f)/ny))
		
		.plot_start(sprintf("%s/all.png", odir), w=nx * panel_size, h=ny * panel_size)
		
		layout(matrix(1:(nx*ny), ny, nx, byrow=T))
		par(mar=c(0.5,0.5,2,0.5))
	}
	
	for (meta_field_v in names(c_by_f)) { 
		ind = is.element(scl2d@scl@scmat@cells, c_by_f[[meta_field_v]])
	
		if (!single_plot) {
			.plot_start(sprintf("%s/%s.png", odir, meta_field_v), w=panel_size, h=panel_size)
			par(mar=c(0.5, 0.5, 2, 0.5))
		}	
		
		plot(scl2d@x, scl2d@y, pch=19, col=bg_col, 
				 cex=cex, type='p',
				 xlim=c(xmin,xmax),
				 ylim=c(ymin,ymax), xlab="", ylab="", xaxt='n', yaxt='n')
		title(main=meta_field_v, cex=cex*3)
		
		fg_cols = if (is.null(fg_col)) { scl2d@c_colors[ind] } else { fg_col }
		points(scl2d@x[ind], scl2d@y[ind], pch=19, col=fg_cols, cex=cex)
		
		if (!single_plot) {
			dev.off()
		}
	}
	
	if (single_plot) {
		dev.off()
	}

	return(scl2d)
}

#' @importFrom KernSmooth bkde2D
#'
#' @export
#'
scp_plot_gene_2d = function(scl2d,
	gene_nm = NA, w, h, base_dir = NULL,
	vals_to_plot = NULL,
	vals_label = NULL,
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
	max_nchars = 15,
	rna_shades = colorRampPalette(c("white", "white", "lightgray", "darkorange1", "darkgoldenrod1", "darkgoldenrod4")),
	pt_shades = colorRampPalette(c("white", "white", "lightgray", "darkgray", "orange", "burlywood1", "chocolate4")))
{
	graph_x = scl2d@x
	graph_y = scl2d@y
	
	if(is.na(gene_nm)) {
		if(is.null(vals_to_plot)) {
			stop("plot gene 2d need either a gene name or values to project")
		}
		rna = vals_to_plot
		gene_nm = vals_label
	} else {
		rna = scl2d@scl@scmat@mat[gene_nm,names(graph_x)]
		
		rna_tot = colSums(scl2d@scl@scmat@mat[,names(graph_x)])
		med_n = median(rna_tot)
		rna_base = floor(med_n*rna/rna_tot)
		rna = rna_base + sapply(med_n*rna/rna_tot-rna_base, function(p) rbinom(1,1,p))
		
	}
	

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
		if(nchar(fnm) > max_nchars) {
			fnm = substr(gene_nm, 1, max_nchars)
		}
		fnm = sub("/", "_", gene_nm)
		fn = sprintf("%s/%s.png", base_dir, fnm)
		.plot_start(fn, h=h, w=w)
		par(mar=c(0,0,3,0))
	}
	rna_x = unlist(apply(cbind(graph_x,rna),1,function(x) rep(x[1],each=x[2]*reg_factor)))
	rna_y = unlist(apply(cbind(graph_y,rna),1,function(x) rep(x[1],each=x[2]*reg_factor)))

	fore <- bkde2D(rbind(pt_reg, cbind(rna_x, rna_y)),
		b=c(rngx/bw_bins, rngy/bw_bins),
		gridsize=c(500, 500))
	back <- bkde2D(rbind(pt_reg, cbind(graph_x, graph_y)),
		b=c(rngx/bw_bins, rngy/bw_bins),
		gridsize=c(500, 500))
	fore$fhat = fore$fhat * sum(back$fhat)/sum(fore$fhat)
#		smoothScatter(rna_x, rna_y, colramp=rna_col, xlim=xl, ylim=yl)

	lrs = log2(fore$fhat/back$fhat)
	if(median(lrs,na.rm=T)>-1) { #background regul is too dominant
		lrs = lrs - median(lrs,na.rm=T) - 1
	}
	# message("plot ", gene_nm, " tot ", sum(rna))
	lrs = pmax(pmin(lrs,4),-1)
	image(x=back$x1, y=back$x2, z=lrs,zlim=c(-1,4), col=rna_shades(1000),xaxt='n', yaxt='n', main=gene_nm)
	if(cont_levels > 0) {
		contour(x=back$x1, y=back$x2, z=lrs, nlevels=cont_levels, lwd=3,add=T, drawlabels=F)
	}

	low_umi_g = ifelse(low_umi == -1, floor(median(rna)), low_umi)
	high_rna_t = as.numeric(quantile(rna, top_color_q))
  #pt_cols = pt_shades[ifelse(rna <= low_umi_g, 1, ifelse(rna <= mid_umi, 2, ifelse(rna <= high_rna_t, 3,4)))]
	#pt_cols = pt_shades(1000)[round(1+999/length(rna) * rank(pmax(rna,min_rna_to_color), ties.method="min"))];
	pt_cols = pt_shades(1000)[.vals_to_n(rna, 1000, pmax(range(rna), min_rna_to_color))]

#			cex=ifelse(pt_cols=="white", 0.5, 1),
	points(graph_x,
		graph_y, pch=21,
		bg=pt_cols,
		cex=ifelse(rna>0, positive_psize, negative_psize),
		col=ifelse(pt_cols=="white", "black", "black"),
		lwd=0.5)
	#grid(lwd=3,col="black",lty=1)

	if(!is.null(base_dir)) {
		dev.off()
	}

	invisible(scl2d) # for competability
}

#
# Internal algs and utility functions
#

.graph_con_comp = function(amat)
{
	if(nrow(amat) > 500) {
		message("graph_con_comp work on small matrices only")
		retrun(NA)
	}
	diag(amat) = TRUE
	ps = amat/rowSums(amat)
	for(i in 1:20) {
		ps = ps %*% ps
	}
	hc = hclust(dist(ps>0))
	compo = cutree(hc, h=0.01)
	return(compo)
}

#
# Serialization - to be updated!!!
#

#' @export
setMethod(
	"tgscm_export",
	signature = "tgScMatClustLayout2D",
	definition =
	 function(.Object, file, supress_clust=F, ...) {
		if(!supress_clust) {
			tgscm_export(.Object@scl, file, ...)
		}
		if(!is.null(.Object@c_ord)) {
			c_ord_fn = sprintf("%s.c_ord", file);
			write.table(data.frame(x=.Object@c_ord), c_ord_fn, sep="\t", quote=F)
		}
		if(!is.null(.Object@x_cl)) {
			c_coord_fn = sprintf("%s.c_coord", file);
			cl_coord_fn = sprintf("%s.cl_coord", file);
			write.table(data.frame(x=.Object@x, y=.Object@y, sd=.Object@cell_2d_sd), c_coord_fn, sep="\t", quote=F)
			write.table(data.frame(x_cl=.Object@x_cl, y_cl=.Object@y_cl), cl_coord_fn, quote=F, sep="\t")
		}
	 }
)

#' Read a clustering layout from file.
#'
#' @param file Name of the file for inputting.
#'
#' @param scl an instance of the relevant \code{\linkS4class{tgScMatClust}}, optional.
#'
#' @export
#'
import_layout = function(file, scl=NULL, ...) {
  .Object = tgScMatClustLayout2D()
  if(is.null(scl)) {
    .Object@scl = import_clust(file)
  } else {
    .Object@scl = scl
  }

  c_ord_fn = sprintf("%s.c_ord", file);

  if(file.exists(g_ord_fn) & file.exists(c_ord_fn)) {
    c_ord = read.table(c_ord_fn, h=T, sep="\t", stringsAsFactors=F)
    .Object@c_ord = c_ord$x
  }
  c_coord_fn = sprintf("%s.c_coord", file);
  if(file.exists(c_coord_fn)) {
    c_coord = read.table(c_coord_fn, h=T, sep="\t", stringsAsFactors=F)
    .Object@x = c_coord$x
    names(.Object@x) = rownames(c_coord)
    .Object@y = c_coord$y
    names(.Object@y) = rownames(c_coord)
    .Object@cell_2d_sd = c_coord$sd
    names(.Object@cell_2d_sd) = rownames(c_coord)
  }
  cl_coord_fn = sprintf("%s.cl_coord", file);
  if(file.exists(cl_coord_fn)) {
    cl_coord = read.table(cl_coord_fn, h=T, sep="\t", stringsAsFactors=F)
    .Object@x_cl = cl_coord$x_cl
    names(.Object@x_cl) = rownames(cl_coord)
    .Object@y_cl = cl_coord$y_cl
    names(.Object@y_cl) = rownames(cl_coord)
  }
  return(.Object)
}

#=========.
#this should be moved to cluster.r
#=========.
scc_get_ordered_marks = function(scl, cl_ord,
		genes_pool=NULL, blacklist_genes=c(), remove_blist_genes = T,
		per_clust_genes=5,
		gene_min_fold = 2.5, gene_min_cov = 0.5)
{
	gene_folds = scl@clust_fp

	if(is.null(genes_pool)) {
		genes_pool = rownames(gene_folds)
	} else {
		genes_pool = intersect(genes_pool, rownames(gene_folds))
	}
	if(is.null(genes_pool) | length(genes_pool) == 0) {
		message("bad gene pool, does not intersect cluster footprint gene names")
		return(NULL)
	}
	if (remove_blist_genes) {
		genes_pool = setdiff(genes_pool, blacklist_genes)
	}
	mask_folds = gene_folds[genes_pool,]*(scl@clust_gcov[genes_pool,]>gene_min_cov)

	good_marks = unique(as.vector(unlist(
			apply(mask_folds,
				2,
				function(x)  {
				   names(head(sort(-x[x>gene_min_fold]),n=per_clust_genes)) })
		     )))

	if(is.null(good_marks) | length(good_marks) < 2) {
		message("no marks found")
		return(NULL)
	}

	.ord_genes = function(gf, cl_ord) {
		n_g = nrow(gf)
		
		genes_hc = hclust(dist(cor(t(gf))), "ward.D2")
		n_gmod = min(ceiling(n_g / min(8, n_g / 3)), n_g)
		hc_gmods = cutree(genes_hc, n_gmod)
		mean_e = apply(gf, 2, function(x)
			tapply(x, hc_gmods, mean))
		genes_ord = order(apply(mean_e[hc_gmods, cl_ord], 1, which.max))
		
		return(rownames(gf)[genes_ord])
	}
	
	if (remove_blist_genes) {
		marks = .ord_genes(gene_folds[good_marks, ], cl_ord)
	}
	else {
		marks = .ord_genes(gene_folds[setdiff(good_marks, blacklist_genes), ], cl_ord)
		bl_marks = intersect(good_marks, blacklist_genes)
		if (length(bl_marks) > 1) {
			bl_marks =  .ord_genes(gene_folds[bl_marks,], cl_ord)
		}
		marks = c(bl_marks, marks)
	}
	
	return(marks)
}

#' Plot heatmap of number of cells per cell module, breakdown by a metadata field, with additional metadata attributes
#'
#' @param sc_2d 
#' @param outdir 
#' @param by 
#' @param metadata_fn 
#' @param key_fields 
#' @param add_text_to 
#' @param key_dict 
#' @param na_col 
#' @param cols_list 
#' @param show_legends 
#' @param hm_pal 
#' @param n_shades 
#' @param cex_text 
#'
#' @return
#' @export
#'
#' @examples
plot_cells_breakdown_to_mods_by_feature = function(sc_2d, 
																									 outdir = get_param("outdir"),
																									 by_field = get_param("scp_cell_mod_tab_by"), 
																									 metadata_fn = get_param("scm_batch_table_fn"), 
																									 key_fields = get_param("scp_cell_mod_tab_key_fields"), 
																									 add_text_to = get_param("scp_cell_mod_tab_add_txt_to"),  
																									 key_dict = get_param("scp_cell_mod_tab_key_dict"),
																									 na_col = get_param("scp_cell_mod_tab_na_col"), 
																									 cols_list = get_param("scp_cell_mod_tab_keys_col_dict"), 
																									 show_legends = get_param("scp_cell_mod_tab_show_legend_for"), 
																									 hm_pal=RColorBrewer::brewer.pal(n=9, 'YlOrRd'), 
																									 n_shades=200, cex_text=1.5) 
{
	
	ab = read.table(metadata_fn, header=T, sep="\t", stringsAsFactors = F)
	ab = unique(ab[, unique(c(by_field, key_fields))])
	ab[ ab == ""] = NA
	
	scl = sc_2d@scl
	cl_tab = table(scl@scmat@cell_metadata[, by_field], scl@clusts)
	cl_tab_n = cl_tab / rowSums(cl_tab)
	
	ab = ab[ ab[, by_field] %in% rownames(cl_tab), ]
	
	rownames(ab) = ab[, by_field]
	
	ord = hclust(dist(cl_tab_n), method='ward.D2')$order
	
	cl_tab = cl_tab[ord, ]
	cl_tab_n = cl_tab_n[ord, ]
	
	ab = ab[rownames(cl_tab), ]
	
	.plot_start(sprintf("%s/cells_breakdown_to_cmods_by_%s.png", outdir, by_field), w=1600, h=1200)
	
	n_leg = length(show_legends)+1
	n_keys = length(key_fields)
	n_plots = n_keys + n_leg + 2
	
	m = rbind(matrix(1:n_keys, n_leg, n_keys, byrow=T), n_plots + 1)
	m = cbind(m, c(rep(n_keys+1, n_leg), n_keys+2))
	m = cbind(m, c(seq(n_keys+3, n_keys+2+n_leg), n_plots+1))
	
	layout(m, widths = c(rep(1, n_keys), 30, 6), heights = c(rep(20/n_leg, n_leg), 1))
	par(mar=c(1,0.5,12,0.5))
	
	# left keys
	for (key in rev(key_fields)) {
		v = ab[, key]
		if (key %in% names(key_dict)) {
			dict = key_dict[[key]]
			
			if (length(dict) == 1 && dict == "auto") {
				dict_vals = unique(v)
				dict = seq_along(dict_vals)
				names(dict) = dict_vals
				cols_list[[key]] = c(RColorBrewer::brewer.pal(n=9, 'Set1'), RColorBrewer::brewer.pal(n=12, 'Set3'))[1:length(dict)]
			}
			v2 = rep(NA, length(v))
			v2[ v %in% names(dict)] = unlist(dict[ v[v %in% names(dict)] ])
			v2[ is.na(v2) ] = max(v2, na.rm=T) + 1
			v = v2
		}
		else {
			v[ is.na(v) ] = max(v, na.rm=T) + 1
			v = v + 1
		}
		cols = c(cols_list[[key]], na_col)
		
		image(1, seq_along(v), t(seq_along(v)), col = cols[v], xaxt='n', yaxt='n')
		abline(h=seq(1, length(v))-0.5, col='grey90')
		mtext(key, side=3, line=1, las=2, at=1, cex=cex_text)
		if (key %in% add_text_to) {
			text(labels=names(dict[v]), x=1, y=1:length(v), col='grey30', cex=cex_text)
		}
	}
	
	# main matrix
	par(mar=c(1,6,12,40))
	image(1:ncol(cl_tab_n), 1:nrow(cl_tab_n), t(cl_tab_n[, sc_2d@cl_ord]), zlim=c(0, 1), col=colorRampPalette(hm_pal)(n_shades), xaxt='n', yaxt='n', xlab="", ylab="")
	abline(h=1:ncol(cl_tab_n) - 0.5, col='grey90')
	abline(v=1:nrow(cl_tab_n) - 0.5, col='grey90')
	for (x in 1:ncol(cl_tab)) {
		for (y in 1:nrow(cl_tab)) {
			text(x, y, cl_tab[, sc_2d@cl_ord][y,x], col='grey30', cex=cex_text)
		}
	}
	mtext(text=rowSums(cl_tab), side=2, at=1:nrow(cl_tab), las=2, cex=cex_text, line=1)
	mtext(text=colSums(cl_tab), side=3, at=1:ncol(cl_tab), las=2, cex=cex_text, line=1)
	mtext(text="#cells", side=3, at=ncol(cl_tab)/2, las=1, cex=cex_text, line=7)
	mtext(text=rownames(cl_tab), side=4, at=1:nrow(cl_tab), las=2, cex=cex_text, line=1)
	
	# bottom cl colors
	par(mar=c(3,6,1,40))
	image(1:ncol(cl_tab), 1, matrix(1:ncol(cl_tab), ncol=1), col=sc_2d@cl_colors[sc_2d@cl_ord], xaxt='n', yaxt='n', xlab="", ylab="")
	abline(v=1:ncol(cl_tab) - 0.5, col='grey90')
	mtext(text=sc_2d@cl_ord, side=1, at=1:ncol(cl_tab), las=2, cex=cex_text, line=1)
	
	# legends
	par(mar=c(4,4,4,16))
	
	image(t(as.matrix(1:n_shades)), col=colorRampPalette(hm_pal)(n_shades), xaxt='n', yaxt='n')
	mtext(text="%", side=3, cex=cex_text, line=1)
	mtext(text=c(0, 100), side=4, at=c(0, 1), cex=cex_text, las=2, line=1)
	
	for (key in show_legends) {
		labs = names(key_dict[[key]])
		cols = cols_list[[key]]
		image(1, 1:length(labs), t(1:length(labs)), col=cols, xaxt='n', yaxt='n', xlab="", ylab="")
		abline(h=1:length(labs) - 0.5, col='grey90')
		mtext(text=key, side=3, cex=cex_text, line=1)
		mtext(text=labs, side=4, at=1:length(labs), cex=cex_text, las=2, line=1)
	}
	
	dev.off()
}