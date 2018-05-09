# this is important for running this package from scripts
#' @import methods
NULL

# #' @example R/workflow_example.r

#' Single cell RNA-seq matrix
#'
#' a single cell RNA matrix interface. Does not do much beyond adding
#' specifically gene names, cell names and labels describing the type of data.
#' Methods for basic stat, gene selection, clustering are implemented separately
#'
#' @slot mat A sparse matrix containing the expression data
#' @slot Genes list of all genes
#' @slot Cells list of all cells
#' @slot stat_type Type of statistic. ("umi", "rpkm", "tpm" etc.)
#' @slot feature_type .
#' @slot cell_metadata dataframe with metadata on the cells. Rownames correspond to the cell names.
#'
#' @examples
#' # we first initialize a tgScMat object
#' scmat = tgScMat(mat=scrdb::toydata, stat_type="umi")
#'
#' # we next compute gene markers, while also dumping stats into a figure gset.png
#' # this stage likely requires intervention by tuning parameters
#' marks = tgscm_gene_select(scmat)
#'
#' # we next cluster using the knn graph approach and default parameters
#' scl = tgScMatClust(scmat=scmat, feat_mat=scmat@mat[marks,], alg_type="knn_graph")
#' scl_kmeans = tgScMatClust(scmat=scmat, feat_mat=scmat@mat[marks,], alg_type="kmeans")
#'
#' # we then compute the graph/cells 2d layout
#' scl_2d = tgScMatClustLayout2D(scl=scl, alg_type="graph", T_edge=0.08)
#'
#' # optionally you can save the data, clusteirng and layout for later loading (note that the sc mat
#' # will be duplicated - altohugh options to prevent this are available
#' tgscm_export(scl_2d, file="fig/toydata")
#'
#' # plotting the clusters and cells around them - with some default coloring scheme that will
#' # usually won't cut it (unless the model is really small)
#'
#' tgscm_plot_clust_2d(scl_2d)
#'
#' # we can now plot the density map for each marker genes
#' for(nm in marks) {
#'   tgscm_plot_gene_2d(scl_2d,gene_nm=nm, reg_factor=5, w=2000,h=2000,base_dir="fig/")
#' }
#'
#' @importClassesFrom Matrix Matrix dgCMatrix dgeMatrix
#' @import Matrix
#'
#### @export tgScMat
#### @exportClass tgScMat
tgScMat <- setClass(
        "tgScMat",
	slots = c(
	  mat = "dgCMatrix",
	  genes = "vector",
	  cells = "vector",
	  ncells = "numeric",
	  ngenes = "numeric",
	  stat_type = "character",
		feature_type = "character",
	  cell_metadata = "data.frame")
)

setMethod(
  "initialize",
  signature = "tgScMat",
  definition =
    function(.Object, mat=NULL, stat_type = "umi", feature_type = "all", cell_metadata = NULL, ...) {
      .Object@stat_type=stat_type
      .Object@feature_type=feature_type
      # if(!is.null(batch_map)) {
      #    .Object@batch_map = batch_map
      # }
      # mat
      if(!is.null(mat)) {
        if(!is.null(rownames(mat))) {
          .Object@genes=rownames(mat)
        } else {
          .Object@genes=seq(1,nrow(mat),1)
        }
        if(!is.null(colnames(mat))) {
          .Object@cells=colnames(mat)
        } else {
          .Object@cells=seq(1,ncol(mat),1)
        }
      	.Object@ngenes = dim(mat)[1]
      	.Object@ncells = dim(mat)[2]

      	.Object@mat = Matrix(as.matrix(mat))

      	# cell metadata
      	if(!is.null(cell_metadata)) {
      	  .Object@cell_metadata = .process_md(.Object, cell_metadata)
      	}
      }


      if(length(list(...)) > 0) {
        # ...
      }
      return(.Object)
    }
)

# Checks metadata
.process_md = function(.Object, cell_metadata) {
  if(is.null(rownames(cell_metadata))){
    if (dim(cell_metadata)[1]==length(.Object@cells)){
      rownames(cell_metadata) = .Object@cells
    } else {
      stop("Error in tgScMat construction. No cell names supplied for cell_metadata, and the number of cells differ")
    }
  } else {
    # select relevant cells
    if( length(intersect(rownames(cell_metadata), .Object@cells)) ==0) {
      stop("Error in tgScMat construction. Cell names in cell_metadata do not match expression matrix cell names")
    } else {
      if (length(setdiff(.Object@cells, rownames(cell_metadata)))>0) {
        warning("Some cells are missing metadata")
      }
      cell_metadata = cell_metadata[.Object@cells,]
      rownames(cell_metadata) = .Object@cells # in case some are missing
    }

  }
  return (cell_metadata)
}



setValidity("tgScMat", function(object) {
  # TBA
  return(TRUE)
})



# PUBLIC CONSTRUCTORS

#' Constract a tgScMat
#'
#' Constract a tgScMat gene expression matrix object
#'
#'
#' @param mat A matrix containing the expression data (will be stored as a sparse matrix)
#' @param cell_metadata dataframe with metadata on the cells. Rownames correspond to the cell names.
#' @param stat_type Type of statistic. ("umi", "rpkm", "tpm" etc.). Default: "umi"
#' @param feature_type Currently supporting UMI tables (feature_type = "all")

#'
#' @return A tgScMat gene expression matrix object
#'
#' @export
scm_new_matrix = function(mat,
													cell_metadata,
													stat_type = get_param("scm_stat_type"),
													feature_type = get_param("scm_feature_type")) {
	return(
		tgScMat(
			as.matrix(mat),
			stat_type = stat_type,
			feature_type = feature_type,
			cell_metadata = cell_metadata
		)
	)
}


#' Read a matrix from the output of a MARS-Seq run
#'
#' @param base_dir base directory of MARS pipeline
#' @param mars_batches all mars technical (amplification) batch codes. If null, will read all batches listed in batch_metadata.
#' @param batch_meta_file path to a tab delimited file describing the mars batches. First column must be batch name. Exactly one of \code{batch_meta_file}, \code{batch_meta} should be used.
#' @param batch_meta a data frame describing the mars batches.
#' Exactly one of \code{batch_meta_file}, \code{batch_meta} should be used.
#' @param skip_missing Skip batches that are present in the metadata file but a missing from the UMI directory? (default: T)
#' @param min_umis_n minimum number of umi's for retaining the cell in the matrix (default: 200)
#' @export
#'
#' @importFrom data.table fread
#'
scm_read_scmat_mars = function(base_dir = get_param("scm_mars_base_dir"),
				 mars_batches = get_param("scm_mars_batches"),
				 batch_meta_file = NULL,
				 batch_meta = NULL,
				 skip_missing = get_param("scm_skip_missing_batches"),
				 genes_by_first = T,
				 min_umis_n = get_param("scm_min_cell_umis"))
{
	batch_meta = .get_md(batch_meta, batch_meta_file)

	if (is.null(mars_batches)) {
		mars_batches = batch_meta[,1]
	} else {
		# check
		dif = setdiff(mars_batches, batch_meta[,1])
		if (length(dif) > 0) {
			stop ("Metadata missing for batches ", paste(dif, collapse=" "))
		}
	}
	# umis$Row.names = rownames(umis) done in fread_rownames
	cells_in_batch = c()
	good_batches = c()
	batches_stats = data.frame()
	
	for(i in 1:length(mars_batches)) {
		batch = as.character(mars_batches)[i]
		fn = sprintf("%s/%s.txt", base_dir, batch)
		if(file.exists(fn) | !skip_missing) {
			good_batches = c(good_batches, batch)
			bumis = fread_rownames(sprintf("%s/%s.txt", base_dir, batch), sep="\t")
			batches_stats = rbind(batches_stats, .calc_batch_stats(bumis, min_umis_n))
			
			if(length(cells_in_batch)==0) {
				umis = bumis
				cells_in_batch = c(ncol(umis) - 1)
			} else {
				if(genes_by_first) {
					if(nrow(umis) != nrow(bumis)) {
						message("row mismatch in loading mars batch ", batch, " reverting to join")
						umis = umis %>% left_join(bumis)
					} else {
						umis = cbind(umis,bumis[,-1])
					}
				} else {
					umis = umis %>% left_join(bumis)
				}
				cells_in_batch = c(cells_in_batch, ncol(bumis) - 1)
			}
			cat ("read", batch, ",dim ", paste(dim(umis), collapse=" "), "\n")
			rownames(umis) = umis$rowname
		} else {
			cat ("skip missing ", batch)
		}
	}
	no_ercc = grep("ERCC", rownames(umis), invert=T)
	umis = umis[no_ercc, -1] # remove ERCC & first col
	umis = umis[rowSums(umis)>1,]

	rownames(batch_meta) = batch_meta[,1]
	cell_metadata = batch_meta[rep(as.character(good_batches), times = cells_in_batch),]

	# filter small cells
	f = colSums(umis)>min_umis_n
	umis = umis[,f]
	cell_metadata = cell_metadata[f,]
	cat("Filtered", sum(!f), "small cells\n")

	rownames(cell_metadata) = colnames(umis)


	# report batch stats
	rownames(batches_stats) = good_batches
	write.table(batches_stats, sprintf("%s/batches_stats.txt", get_param("outdir")), quote=F, sep="\t")
	.plot_batch_stats(batches_stats, batch_meta , sprintf("%s/batches_stats.png", get_param("outdir")))

	return(tgScMat(as.matrix(umis), stat_type="umi", cell_metadata = cell_metadata))
}


# get metadata
.get_md = function(batch_meta, batch_meta_file) {
	if (is.null(batch_meta) & is.null(batch_meta_file)) {
		stop("Batch metadata is not specified")
	}
	if (!is.null(batch_meta) & !is.null(batch_meta_file)) {
		stop("Exactly one of batch_meta_file, batch_meta should be used.")
	}
	if (!is.null(batch_meta_file)) {
		batch_meta = read.delim(batch_meta_file)
	}

	batch_meta = as.data.frame(unclass(batch_meta)) # convert char to levels
	return (batch_meta)
}

#' calc batch stats
#'
#' @param bumis loaded umis table
#' @param min_count threshold for valid cells
#'
#' @return batch stats table
.calc_batch_stats = function(bumis, min_count=200) {
	
	if (colnames(bumis)[1] == 'rowname') {
		rownames(bumis) = bumis$rowname
		bumis = bumis[, -1]
	}
	
	ercc = grep("ERCC", rownames(bumis))
	
	if (length(ercc) == 0) {
		total_ercc = 0
	}
	else {
		total_ercc = sum(bumis[ercc, ])	
		bumis = bumis[-ercc, ]
	}
	
	tot_c = colSums(bumis)
	tot_valid_c = tot_c[tot_c >= min_count]
	
	data.frame(
		total = sum(bumis),
		total_ercc = total_ercc,
		n_lt50 = sum(tot_c < 50),
		n_50_100 = sum(tot_c >=  50 & tot_c < 100),
		n_100_200 = sum(tot_c >=  100 & tot_c < 200),
		n_200_400 = sum(tot_c >=  200 & tot_c < 400),
		n_400_800 = sum(tot_c >=  400 & tot_c < 800),
		n_800_1600 = sum(tot_c >=  800 & tot_c < 1600),
		n_gt1600 = sum(tot_c >= 1600), 
		n_valid = length(tot_valid_c),
		valid_q25 = quantile(tot_valid_c, 0.25),
		valid_q50 = quantile(tot_valid_c, 0.5),
		valid_q75 = quantile(tot_valid_c, 0.75)
		)
}

#' plot batches stats
#'
#' @param batches_stats collected by .calc_batch_stats 
#' @param batches_meta metadata
#' @param ofn output file name 
#' @param cells_breakdown_field color #cells per batch by this field 
#' @param cell_col_dict color dictionary for cells_breakdown_field values
.plot_batch_stats = function(batches_stats, batches_meta, ofn, cells_breakdown_field=get_param("scm_cells_breakdown_field"), cell_col_dict=get_param("scm_cell_col_dict")) 
{
	
	batches_stats = batches_stats[order(gsub("AB", "", rownames(batches_stats)), decreasing = T), ]
	batches_meta = batches_meta[rownames(batches_stats), ]
	
	b_tot = batches_stats$total
	names(b_tot) = rownames(batches_stats)
	
	main_cex = 1
	bp_cex = 0.5
	
	png(ofn, w=960, h=1200, pointsize=7, res=300)
	layout(matrix(1:4, 1, 4), w=c(2,1,5,1))
	
	type_cols = 'blue'
	if (!is.null(cell_col_dict) & !is.null(cells_breakdown_field)) {
		type_cols = unlist(cell_col_dict[batches_meta[, cells_breakdown_field]])
	}

	par(mar=c(4,5,2,0))
	barplot(b_tot, horiz=T, las=2, cex.axis=bp_cex, cex.names=bp_cex, border=NA, col=type_cols)
	title(main='total UMIs', cex.main=main_cex)
	if (!is.null(cell_col_dict) & !is.null(cells_breakdown_field)) {
		legend("topleft", legend=names(type_cols), bty='n', border=NA, fill=unlist(cell_col_dict), cex=bp_cex)
	}
	box(lwd=0.5)
	
	par(mar=c(4,1,2,0))
	barplot(batches_stats$total_ercc/batches_stats$total, cex.axis=bp_cex, horiz=T, border=NA, col='black')
	title(main='%ERCC', cex.main=main_cex)
	box(lwd=0.5)
	
	bd_cols = RColorBrewer::brewer.pal(n=7, 'Blues')
	barplot(t(batches_stats[,3:9]), horiz=T, border=NA, cex.axis=bp_cex, col=bd_cols, yaxt='n')
	title(main='#cells by #UMIs', cex.main=main_cex)
	legend("topleft", legend=gsub("_", "-", gsub("n_", "", colnames(batches_stats)[3:9])), bty='n', fill=bd_cols, cex=1.2 * bp_cex, pt.cex=1.2 * bp_cex, ncol=length(bd_cols), border=NA, x.intersp=0.8)
	box(lwd=0.5)
	
	par(mar=c(4,1,2,1))
	barplot(batches_stats$n_valid, horiz=T, border=NA, cex.axis=bp_cex, col='black')
	title(main='#valid', cex.main=main_cex)
	box(lwd=0.5)
	
	dev.off()
}

#' Read a matrix from the output of a 10x run. Extracting batch info from the cell barcodes (e.g. AACAGTGCGTGA-1 is in batch 1).
#'
#' @param sparse_matrix_ifn sparse matrix input file name. Expecting a tab delimited file with 3 columns and a header: 'row', 'column', 'value'. row and column fields are the row and column number of the value. Expecting 3 other files (if matrix file name is xxx.tab): xxx.rownames.tab, xxx.colnames.tab (each header by 'rowname' or 'colname'), with the row/column names. xxx.dims.tab (headed by 'dims') with 2 rows - the nmber of rows and the number of columns.
#' @param use_batches 
#' @param min_umis_n minimum number of umi's for retaining the cell in the matrix (default: 200)
#'
#' @export
#'
#' @importFrom data.table fread
#'
scm_read_scmat_10x = function(sparse_matrix_ifn = get_param("scm_10x_mat_ifn"),
															genes_ifn = get_param("scm_10x_mat_rownames_ifn"),
															cells_ifn = get_param("scm_10x_mat_colnames_ifn"),
															paralogs_policy = get_param("scm_10x_paralogs_policy"),
															use_batches = get_param("scm_10x_batches"),
															min_umis_n = get_param("scm_min_cell_umis"))
{
	# read sparse matrix with all genes and batches
	umis = fread_mm(fname = sparse_matrix_ifn, row.names = genes_ifn, col.names = cells_ifn)
	
	genes = read.table(genes_ifn, header=F, stringsAsFactors = F)
	colnames(genes) = c('id', 'name')
	rownames(genes) = genes$id
	
	unique_genes = names(which(table(genes$name) == 1))
	umis1 = umis[genes$name %in% unique_genes, ]
	rownames(umis1) = genes[rownames(umis1), 'name']
	if (paralogs_policy == 'remove') {
		umis = umis1
	}
	else if (paralogs_policy == 'sum') {
		non_unique_genes = setdiff(genes$name, unique_genes)
		umis2 = as.matrix(umis[genes$name %in% non_unique_genes, ])
		
		message(sprintf("summing up total of %d paralog genes into %d unique genes", nrow(umis2), length(non_unique_genes)))
		
		umis2s = apply(umis2, 2, function(x) { tapply(x, INDEX=genes[rownames(umis2), 'name'], FUN=sum) } )
		
		umis = rbind(umis1, umis2s)
	}
	else {
		stop(sprintf("Loading 10x data, unknown paralogs policy (%s), supprting remove/sum", paralogs_policy))
	}
	
	md = data.frame(row.names=colnames(umis), type='10x', batch=gsub(".*-", "", colnames(umis)))
	
	# filter batches if needed
	if (!is.null(use_batches)) {
		ind = md$batch %in% use_batches
		
		umis = umis[, ind]
		md = md[ind, ]
	}
	
	# calc batch stats
	batches_stats = data.frame()
	for(b in sort(unique(md$batch))) {
		batches_stats = rbind(batches_stats, .calc_batch_stats(umis[, md$batch == b], min_umis_n))
	}

	no_ercc = grep("ERCC", rownames(umis), invert=T)
	umis = umis[no_ercc, ] # remove ERCC & first col
	umis = umis[rowSums(umis)>1,]
	
	# filter small cells
	f = colSums(umis) > min_umis_n
	umis = umis[,f]
	md = md[f,]
	cat("Filtered", sum(!f), "small cells\n")
	
	# report batch stats
	rownames(batches_stats) = sort(unique(md$batch))
	write.table(batches_stats, sprintf("%s/batches_stats.txt", get_param("outdir")), quote=F, sep="\t")
	
	batch_meta = unique(md)
	rownames(batch_meta) = batch_meta$batch
	.plot_batch_stats(batches_stats, batch_meta , sprintf("%s/batches_stats.png", get_param("outdir")))
	
	return(tgScMat(as.matrix(umis), stat_type="umi", cell_metadata = md))			
}
	
#' Extract sub-matrix
#'
#' @param scmat A tgScMat object.
#' @param genes Genes range.
#' @param cells Cells range
#'
#' @export
scm_sub_mat = function(scmat, genes=NULL, cells=NULL) {
	if (class(scmat)[1] != "tgScMat") {
		stop("invalid scmat in scm_sub_mat")
	}
	if (is.null(genes)) {
		genes = scmat@genes
	}
	if (is.null(cells)) {
		cells = scmat@cells
	}

	if(length(cells)<2) {
		stop("At least 2 cells must be selected")
	}
	if(length(genes)<2) {
		stop("At least 2 genes must be selected")
	}
	meta = scmat@cell_metadata[cells,]
	# remove unused levels
	meta[] <- lapply(meta, function(x) if(is.factor(x)) factor(x) else x)
	return(tgScMat(mat = scmat@mat[genes, cells,drop=F], stat_type= scmat@stat_type, feature_type = scmat@feature_type,
								 cell_metadata= meta))
	#batch_map = scmat@batch_map[cells],

}


#' Read a tgScMat object that saved to disk.
#'
#' @param file Prefix of file names input (expecting %s.scmat and %s.scmat_cell_metadata and $s.scmat_attr to exist)
#'
#' @export
#'
scm_import_mat = function (file) {
	mat = fread_rownames(file=sprintf("%s.scmat", file), sep="\t", set_rownames = T)
	attr = as.list(read.table(file=sprintf("%s.scmat_attr", file), sep="\t", stringsAsFactors = F))
	cell_metadata = fread_rownames(file=sprintf("%s.scmat_cell_metadata", file), sep="\t", set_rownames = T)

	if(is.null(attr$stat_type)) {
		attr$stat_type = "undef"
	}
	if(is.null(attr$feature_type)) {
		attr$feature_type = "undef"
	}
	return(tgScMat(as.matrix(mat), stat_type=attr$stat_type, feature_type=attr$feature_type, cell_metadata=cell_metadata))
}


#'
#' Export a matrix to file.
#'
#' @param file Prefix of file names for outputting.
#'
#' @export
#'
#'
#'
setGeneric("tgscm_export",
           function(.Object, file,...) stnadrdGeneric("tgscm_export"))
setMethod(
  "tgscm_export",
  signature = "tgScMat",
  definition =
    function(.Object, file, ...) {
      write.table(as.matrix(.Object@mat), file=sprintf("%s.scmat", file), quote=F, sep="\t")
      scmat_attr = list(stat_type = .Object@stat_type, feature_type = .Object@feature_type)
      write.table(scmat_attr, file=sprintf("%s.scmat_attr", file), quote=F, sep="\t")
      # metadata
      write.table(as.matrix(.Object@cell_metadata), file=sprintf("%s.scmat_cell_metadata", file), quote=F, sep="\t")

    }
)
#'
#' Clean ambient noise from matrix
#' 
#' @param scmat A tgScMat object.
#' @param epsilon 
#' @param metadata_batch 
#' @param min_umi_n 
#'
#' @export
scm_remove_ambient_by_epsilon = function(scmat, 
																				 epsilon = get_param("scm_amb_noise_epsilon"), 
																				 metadata_batch = get_param("scm_batch_meta_attr"), 
																				 min_umi_n = get_param("scm_min_cell_umi"))
{
	#iter on batches
	#epsilon*U/N = something. the old heuristc is to remove below T as long as sum(ui<T) < epsilon_upper*U.The problem is that it creates.
	us = scmat@mat
	batch_factor = scmat@cell_metadata[,metadata_batch]
	for(b in unique(batch_factor)) {
		message("cleaning ", b)
		cells = which(batch_factor == b)
		ncell = length(cells)
		tot_g = rowSums(us[,cells])
		thresh_g = tot_g*epsilon/ncell
		thresh_u = us[,cells] < thresh_g & us[,cells] > 0
		us[,cells][thresh_u] = 0
	}
	scmat@mat = us
	if(min_umi_n > 0) {
		scmat = scm_sub_mat(scmat, scmat@genes, names(which(colSums(us)>min_umi_n)))
	}
	return(scmat)
}

####

setMethod(
	"show",
	signature = "tgScMat",
	definition =
		function(object) {
			cat("An object of class ", class(object), ",", sep = "")
			cat(" stat type ", object@stat_type, ", feat select ", object@feature_type, ".\n", sep = "")
			cat(length(object@cells), " cells by ", nrow(object@mat), " genes. median cell content ", median(colSums(object@mat)), ".\n", sep = "")
			invisible(NULL)
		}
)


#' (optional) Load FACS indices info and add it to cell_metadata table
#'
#' @param scmat 
#' @param base_dir 
#' @param facs_idx_batches 
#' @param wells2cells_fn 
#' @param wells2cells_batch_field 
#' @param wells2cells_cell_field 
#' @param wells2cells_well_field 
#'
#' @return
#' @export
#'
#' @examples
scm_add_facs_indices_from_tables = function(scmat,
																						base_dir = get_param("scm_facs_idx_base_dir"),
																						facs_idx_batches = get_param("scm_facs_idx_batches"), 
																						wells2cells_fn = get_param("scm_facs_idx_wells2cells_fn"),
																						
																						
																						
																						
																						wells2cells_batch_field = get_param("scm_facs_idx_batch_field"),
																						wells2cells_cell_field = get_param("scm_facs_idx_cell_field"),
																						wells2cells_well_field = get_param("scm_facs_idx_well_field"), 
																						skip_missing_batch_facs_indices = get_param("scm_facs_idx_skip_missing"))

{
	
	all_batches = unique(scmat@cell_metadata[, get_param("scm_batch_meta_attr")])
	
	if (is.null(facs_idx_batches)) {
		facs_idx_batches = all_batches
	} else {
		dif = setdiff(facs_idx_batches, all_batches)
		if (length(dif > 0)) {
			stop ("FACS indices supplied for batches not loaded ", paste(dif, collapse=" "))
		}
	}
	
	wells_dict = read.table(sprintf("%s/%s", base_dir, wells2cells_fn), header=T, sep="\t", stringsAsFactors = F)
	
	rownames(wells_dict) = wells_dict[, wells2cells_cell_field]
	stopifnot(all(is.element(rownames(scmat@cell_metadata), rownames(wells_dict))))
	wells_dict = wells_dict[ rownames(scmat@cell_metadata), c(wells2cells_batch_field, wells2cells_cell_field, wells2cells_well_field)]
	
	idx_tab = NULL
	
	
	for (batch in facs_idx_batches) {
		ifn = sprintf("%s/%s.txt", base_dir, batch)
		if (file.exists(ifn)) {
			facs_idx = read.table(ifn,
														header = T,
														sep = "\t",
														stringsAsFactors = F)
			facs_idx[, wells2cells_batch_field] = batch
			
			if (is.null(idx_tab)) {
				idx_tab = facs_idx
			}
			else {
				idx_tab[, setdiff(colnames(facs_idx), colnames(idx_tab))] = NA
				facs_idx[, setdiff(colnames(idx_tab), colnames(facs_idx))] = NA
				idx_tab = rbind(idx_tab, facs_idx[, colnames(idx_tab)])
			}
		}
		else if (!skip_missing_batch_facs_indices) {
			stop(sprintf("missing FACS indices file %s", ifn))
		}
	}
	message("loaded facs idxs from batches, merging with wells_cells...")
	m = merge(wells_dict, idx_tab, by=c(wells2cells_batch_field, wells2cells_well_field), all.x=T)
	rownames(m) = m[, wells2cells_cell_field]
	scmat@cell_metadata = cbind(scmat@cell_metadata, m[ rownames(scmat@cell_metadata), ])
	
	scmat
}
	
