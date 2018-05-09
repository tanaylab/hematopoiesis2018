#' generate a k_nn graph from a similarity matrix
#
#' @param sim_mat
#' @param m_knn
#'
#'
tg_simmat_to_knn_ordered = function(mat, k_knn)
{
	N = nrow(mat)
	m_knn_ordered = matrix(NA, nrow = k_knn, ncol = N)
	colnames(m_knn_ordered) = rownames(mat)
	for(i in 1:N) {
		knn = head(order(mat[i,], decreasing = T),n=k_knn+1)
		m_knn_ordered[,i] = knn[-1] # remove self
	}
	return(m_knn_ordered)
}

tg_knn_ordered_to_knn_adjs = function(knn_ord, k_knn)
{
	nodes = colnames(knn_ord)
	N = length(nodes)

	m_knn = matrix(0,nrow = N,ncol = N, dimnames = list(nodes, nodes))
	ranks = (1 - (0:k_knn+1)/(k_knn + 2))[-1]
	for(i in 1:N) {
		m_knn[knn_ord[1:k_knn, i], i] = ranks
	}
	return(Matrix(m_knn))
}

tg_knn_ordered_to_subset_knn_adjs = function(knn_ord, k_knn, nodes)
{
	N = length(nodes)

	m_knn = matrix(0,nrow = N,ncol = N, dimnames = list(nodes, nodes))

	subs_knn = knn_ord[,nodes]

	ranks = (1 - (0:k_knn+1)/(k_knn + 2))[-1]
	for(i in 1:N) {
		subset_neig = head(intersect(nodes, subs_knn[,i]),k_knn)
		m_knn[subset_neig, i] = ranks[1:length(subset_neig)]
	}
	return(Matrix(m_knn))
}


tg_knn_ordered_to_knn_balanced_adjs = function(knn_ord, k_knn, k_expand=4, nodes=NULL)
{
	subset = T
	if(is.null(nodes)) {
		nodes = colnames(knn_ord)
		subset = F
	}
	N = length(nodes)

	k_knn_potential = min(nrow(knn_ord), N)

	m_knn = matrix(k_knn_potential, 
			nrow = N,ncol = N, 
			dimnames = list(nodes, nodes))
	if(subset) {
		id_map = rep(-1, ncol(knn_ord))
		id_map[nodes] = 1:length(nodes)
		
		subs_knn = knn_ord[,nodes]
		for(i in 1:N) {
			subset_neig = id_map[subs_knn[,i]]
			subset_neig = subset_neig[subset_neig != -1]
			m_knn[subset_neig,i] = 1:length(subset_neig)
		}
	} else {
		message("fill m_knn with knn_ord")
		for(i in 1:N) {
			m_knn[knn_ord[1:k_knn_potential,i], i] = 1:k_knn_potential
		}
		message("done fill m_knn with knn_ord")
	}

	M = k_knn*k_knn*k_expand

	m_knn_io = pmax(-m_knn * t(m_knn) + M,0)

	A = nrow(m_knn_io)-k_knn*3	#no more than k_nn*3 contribs per row, or incoming neighbots
	m_knn_io_rows = t(apply(m_knn_io, 1, function(x) pmax(rank(x)-A,0)))
	A = nrow(m_knn_io)-k_knn	#no more than k_nn contribs per column, or outgoing neighbors
	m_knn_io = apply(m_knn_io_rows, 2, function(x) pmax(rank(x)-A,0))
	return(Matrix(m_knn_io))
}


tg_simmat_to_knn_balanced_adjs = function(mat, k_knn, k_expand=4)
{
	nodes = rownames(mat)
	N = length(nodes)
	k_knn10 = min(k_knn*10, N -1)
	return(tg_knn_ordered_to_knn_balanced_adjs(
		 tg_simmat_to_knn_ordered(mat, k_knn10), 
		 k_knn,
		 k_expand))
}


tg_piecewise_knn_graph_cover = function(mat, k_knn, 
					k_expand=4,
					is_raw_mat = F,
					min_clust_size, verbose=F, consolidate=T)
{
	if(is_raw_mat) {
		message("computing sim matrix using plain correlation")
		sim_mat = cor(t(mat), pairwise.complete.obs=T)
		colnames(sim_mat) = rownames(mat)
		rownames(sim_mat) = rownames(mat)
	} else {
		sim_mat = mat
	}
	if(nrow(sim_mat) != ncol(sim_mat)) {
		stop("running knn grahp cover with a non-square similarity matrix as input")
	}
	message("computing knn graph from similarity matrix")
	m_e = tg_simmat_to_knn_balanced_adjs(sim_mat, k_knn, k_expand)
	return(tg_piecewise_graph_cover(m_e, min_clust_size, verbose, consolidate))
}

tg_piecewise_knn_graph_cover_bootstrap = function(mat = NULL, 
			k_knn, 
			boot_ratio, 
			N_boot,
			full_knn_ord = NULL,
			compute_coclust = T,
			k_expand=4, verbose=T, min_clust_size)
{
	if(is.null(full_knn_ord)) {
		if(is.null(mat)) {
			stop("Either matrix or knn_ord should be provided to knn graph cover bootstrap")
		}
		message("computing sim mat using cor")
		sim_mat = cor(t(mat), pairwise.complete.obs=T)
		colnames(sim_mat) = rownames(mat)
		rownames(sim_mat) = rownames(mat)

		nodes = rownames(mat)
		N = length(nodes)
		k_knn10 = min(round(k_knn*10/boot_ratio), N -1)
		message("building knn order")
		full_knn_ord = tg_simmat_to_knn_ordered(sim_mat, k_knn10)
	} else {
		N = ncol(full_knn_ord)
	}

	boot_size = round(N * boot_ratio)
	if(compute_coclust) {
		tot_coclust = matrix(0, nrow=N, ncol=N)
		num_trials = matrix(0, nrow=N, ncol=N)
	} else {
		tot_coclust = NULL
		num_trials = NULL
	}
		
	all_clusts = list()
	all_nodes = list()
	for(i in 1:N_boot) {
		message("bootstrap", i)
		boot_nodes = sample(1:N, boot_size)
		boot_nodes = sort(boot_nodes)
		m_e = tg_knn_ordered_to_knn_balanced_adjs(full_knn_ord, k_knn, k_expand, boot_nodes)
		if (verbose) { message("done computing adj mat") }
		clusts = tg_piecewise_graph_cover(m_e, min_clust_size, verbose=verbose, consolidate=T)
		if(!compute_coclust) {
			all_clusts[[length(all_clusts)+1]] = clusts
			all_nodes[[length(all_nodes)+1]] = boot_nodes
		} else {
		   isclust_ci = diag(max(clusts))[,clusts]
		   coclust_ij = t(isclust_ci) %*% isclust_ci
		   if (verbose) { message("will add tot coclus ", paste(dim(coclust_ij),collapse=" ")) }
		   tot_coclust[boot_nodes, boot_nodes] = 
				tot_coclust[boot_nodes, boot_nodes] + 
				coclust_ij
		   if (verbose) { message("will add trials") }
		   num_trials[boot_nodes,boot_nodes ] = 
				num_trials[boot_nodes,boot_nodes] + 1
		}
		#collect co-clust data
	}
	return(list(clusts = all_clusts, nodes = all_nodes, coclust = tot_coclust, num_trials=num_trials))
}
#' Generic graph cluster-cover
#
#'
#' @param m_e the adjacency matrix of the graph to be clustered
#' @param min_clust_size  used to filter small clusters at the end of the process
#'
#' @importClassesFrom Matrix Matrix
#'

tg_piecewise_graph_cover = function(m_e, min_clust_size, verbose=F, consolidate=T)
{
	nodes = rownames(m_e)
	n_nodes = ncol(m_e)

	m_e = Matrix(m_e, sparse=T)

	m_e_01 = m_e>0
	m_e_01 = Matrix(m_e_01, sparse=T)

#now each entry represent the total weight of going
	m_e2 = t(m_e) %*% m_e

	m_e2 = Matrix(m_e2, sparse=T)

#outgoing mat will store the links from clustering seeds to their neighbors
	outgoing_mat = Matrix(0,
			nrow = nrow(m_e), ncol=ncol(m_e), sparse=T)


	rownames(outgoing_mat) = nodes
	colnames(outgoing_mat) = nodes
	outgoing_weight_d1 = rep(0,n_nodes)
	outgoing_weight_d2 = rep(0,n_nodes)
	outgoing_weight_d3 = rep(0,n_nodes)

	incoming_weight_d1 = rep(0,n_nodes)

	if(verbose) {
	      message("start seeding iterations")
	}
#we sample the next seed from those that are
#1. not yet reachable (d1==0)
#2. preferebly not directly adjacent to reachable (d2==0)
#3. preferebly are at distance 3 (d3>0)
# this is intended to enhance packing of the graph by the clusters

#outging_mat - G_i - a graph containing only outgoing edges from seeds
#N(seeds) - all nodes that are adjacent to seeds. (v s.t. deg_{G_i}(v)>0\}
#outgoing weight_d1 - indegree(G_i)
#uncovered neigh[i] - |N(i) \ N(seeds)|
#outgoing_weight_d2[j] - number of paths of the form seed->v<-j in the entire graph
#outgoing_weight_d3[j] - number of paths of the form seed->v<-u<-j in the entire graph
#

	seeds = c()
	uncovered_neigh = colSums(m_e_01)
	cands_gap = which(outgoing_weight_d2 == 0)
	while(length(cands_gap) > 0) {
		seed = cands_gap[sample(length(cands_gap), 1,
				prob=(1e-10+outgoing_weight_d3[cands_gap]))]
		if(verbose) {
		      message("sample gapped seed, ", seed)
		}

		seeds = c(seeds, seed)
		outgoing_mat[seed,] = m_e[,seed]

		new_covs = (outgoing_weight_d1 == 0) & (outgoing_mat[seed,] > 0)
		outgoing_weight_d1 = outgoing_weight_d1 + outgoing_mat[seed,]  #colSums(outgoing_mat)
		if(verbose) { message("new covs ", sum(new_covs)) }
		uncovered_neigh = uncovered_neigh - new_covs %*% m_e_01

		outgoing_weight_d2 = colSums(outgoing_mat %*% m_e)
		outgoing_weight_d3 = colSums(outgoing_mat %*% m_e2)
		if (verbose) message("sampled seed ", seed)

		cands_gap = which(outgoing_weight_d2 == 0)
	}


#everyhing at distance 2 is covered - so we prefer sampling objects with
#more free second neighbors
	while(sum(outgoing_weight_d1==0)>0
	& max(uncovered_neigh[outgoing_weight_d1==0]) >= min_clust_size/2) {
		cands = which(outgoing_weight_d1 == 0)
		#print(cands)
		#print(outgoing_weight_d2[cands])
		if(length(cands) == 1) {
			seed = cands;
		} else {
#		    if(length(cands) > 500) {
#			cands = sample(cands, 500, replace=F)
#		    }
#		    uncovered_neigh = apply(m_e[,cands], 2,
#				 function(x) {
#				    sum(x>0 & outgoing_weight_d1==0) } )
#		    max_uncovered_neigh = max(uncovered_neigh)
#
#		    names(uncovered_neigh) = cands
		    seed = cands[sample(length(cands), 1, prob=uncovered_neigh[cands]**3)]
#		    message("sample d1 with ", sum(outgoing_weight_d1 == 0), " uncovered, #seeds ", length(seeds), " seed ", seed, " with ", uncovered_neigh[as.character(seed)])
		    if(verbose) {
		      message(" seed ", length(seeds), " uncovered ", sum(outgoing_weight_d1 == 0))
		    }
		}
		seeds = c(seeds, seed)
		outgoing_mat[seed,] = m_e[,seed]
		new_covs = (outgoing_weight_d1 == 0) & (outgoing_mat[seed,] > 0)
		if(verbose) { message("new covs ", sum(new_covs)) }
		outgoing_weight_d1 = outgoing_weight_d1 + outgoing_mat[seed,]  #colSums(outgoing_mat)
		uncovered_neigh = uncovered_neigh - new_covs %*% m_e_01
		if(verbose) { message("after update, max uncov ", max(uncovered_neigh)) }
		if (verbose) message("sampled seed ", seed)
#		print(outgoing_weight_d1)
#		print(outgoing_weight_d2)
	}

	outgoing_weight_d2 = colSums(outgoing_mat %*% m_e)
	# message("done sampling seeds")

#	d2 = t(m_e) %*% t(outgoing_mat)

	if (verbose) message("done seeding")
#compute the cluster of each node by looking for the maximum weight seed neighbor
	clusts = apply(outgoing_mat, 2, function(x) {
			 ifelse(max(x)>0, which.max(x), -1)})
	names(clusts) = colnames(outgoing_mat)

	bad_clusts = .find_small_clusts(clusts, outgoing_mat, min_clust_size)
	# message("before killing small clusters ", sum(clusts==-1), " elements are unassigned")
	# message(length(bad_clusts), " small clusters to kill")

	if (length(bad_clusts) > 0) {
	  clusts[clusts %in% bad_clusts] = -1
	  outgoing_mat[bad_clusts,] = 0
	  # message("after killing small clusters ", sum(clusts==-1), " elements are unassigned")
	}

	orphan_members = names(which(clusts == -1))
	if(length(orphan_members) != 0) {
		outgoing_mat2 = outgoing_mat %*% m_e
		clusts2 = apply(outgoing_mat2, 2, which.max)
		names(clusts2) = colnames(outgoing_mat)
		clusts[orphan_members] = clusts2[orphan_members]
	}
	# message("after orphan reassignment by 2 chains ", sum(clusts==-1), " elements are unassigned")
	clusts = as.integer(as.factor(clusts))
	names(clusts) = colnames(outgoing_mat)

	if(consolidate == T) {
		clusts = .graph_clust_consolidate(clusts, m_e, min_clust_size)
	}

	return(clusts)
}

#cancel small clusts and reassign members
.find_small_clusts = function (clusts, outgoing_mat, min_clust_size) {
  cl_sz = table(clusts)
  bad_clusts = as.numeric(names(which(cl_sz<min_clust_size)))
  bad_clusts = bad_clusts[bad_clusts>0]
  return(bad_clusts)
}

# Not tested yet!
.graph_clust_consolidate = function(clusts, m_e, min_clust_size,
			min_diff = NA, max_iter = 40, start_cooling = 10)
{
	C = length(clusts)
	if(is.na(min_diff)) {
		min_diff = 10
	}
	diff = min_diff+1
	new_clusts = clusts
	ncls = max(clusts)
	iter = 0
	reg = 1
	cooling_rate = 1.05

	start_cooling = max(start_cooling, 5)

	diff_traj = rep(0,max_iter)

	while(diff > min_diff & iter < max_iter)	{
		clusts = new_clusts
		clusts1 = diag(ncls)[,clusts]    # 1 encoding of the clusters, rows= clusters, cols =nodes
		csizes = rowSums(clusts1); csizes[csizes==0] = 1 # add 1 to avoid division by 0
		clusts1 = clusts1/csizes

		#vector of forward clust assoc for each member (which clusts do I like?)
		forward_votes = t(clusts1 %*% m_e) *1000
		#vector of reverse clust assoc for each member (how many in the clust like me?)
		backward_votes = tcrossprod(m_e, clusts1) * 1000

		# #vector of forward clust assoc for each member (which clusts do I like?)
		# forward_votes_old = t(apply(m_e, 2, function(x) tapply(x, clusts, sum)/csize))
		# #vector of reverse clust assoc for each member (how many in the clust like me?)
		# backward_votes_old = t(apply(m_e, 1, function(x) tapply(x, clusts, sum)/csize))

		votes = backward_votes * forward_votes + 1e-10 * forward_votes  # + 1e-10 * forward_votes is for nodes that no one wants

		# add regularization to make nodes stay in their cluster
		if (iter > start_cooling) {
			if((diff_traj[iter-1]/diff_traj[iter-3]) > 0.75) {
				reg = reg * cooling_rate
				# message("cooling to ", reg)
			}
			if(reg > 1) {
				idx = nrow(votes)*(clusts-1) + 1:nrow(votes)
				votes[idx] = reg * votes[idx]
			}
		}

		new_clusts = apply(votes,1,which.max)

		diff = sum(new_clusts != clusts)
		diff_traj[iter] = diff
		message(" i ", iter, " d=",diff)
		iter=iter+1
	}
# TODO: add code to eliminate small clusters, reassign orphand, re-number clusters in case some are missing
#	bad_clusts = .find_small_clusts(clusts, outgoing_mat, min_clust_size)
#	if (length(bad_clusts) > 0) {
#	  new_clusts[clusts %in% bad_clusts] = -1
#	}
	message("done consolidating")

	return(new_clusts)
}

