#' Loading and saving scrdb objects
#'
#' @slot basedir character. Directory containing the objects for save and load.
#' @slot fn_prefix character. Filename prefix for the objects managed by the instance. 
#' @slot which_not_lazy vector. Object types not to lazy load (regenerate object even if it can be loaded from the db)
#'
#' @exportClass tgScDbMngr
#' 
tgScDbMngr <- setClass(
	"tgScDbMngr",
	slots = c(
		basedir = "character",
		fn_prefix = "character",
		which_not_lazy = "vector")
)

#' Print scdb object
#'
#' @param tgScGeneMods 
#'
#' @return
#' @export
#'
#' @examples
setMethod(
	"show",
	signature = "tgScDbMngr",
	definition =
		function(object) {
			cat(sprintf("tgScDbMngr\nbasedir: %s\nfile name prefix: %s\nobject not to lazy-load: %s\n", object@basedir, object@fn_prefix, paste(object@which_not_lazy, sep="", collapse=", ")))
			invisible(NULL)
		}
)

#' Initialize tgScDbMngr object
#'
#' @param basedir
#' @param fn_prefix
#' @param which_not_lazy
#'
#' @return invisible tgScDbMngr object
#' @export
#'
#' @examples
scdb_init = function(basedir = get_param("scdb_basedir"),
										 fn_prefix = get_param("scdb_fn_prefix"),
										 which_not_lazy = get_param("scdb_which_not_lazy"))
{
	scdb = new("tgScDbMngr")
	
	scdb@basedir = basedir
	scdb@fn_prefix = fn_prefix
	scdb@which_not_lazy = which_not_lazy
	
	return(scdb)
}

#' Saves an object
#'
#' @param scdb 
#' @param object 
#' @param obj_type 
#' @param overwrite_if_exists 
#'
#' @export
#'
#' @examples
scdb_save = function(scdb, object, obj_type, overwrite_if_exists = T) 
{
	fn = sprintf("%s/%s_%s.RData", scdb@basedir, scdb@fn_prefix, obj_type)
	
	if (file.exists(fn) && !overwrite_if_exists) {
		stop(sprintf("Cannot save %s object, file already exists (%s).", obj_type, fn))
	}
	
	save(object, file=fn)
}

#' Loads an object
#'
#' @param scdb 
#' @param obj_type 
#'
#' @return loaded object or NULL if it doesn't exist or marked not to lazy load
#' @export
#'
#' @examples
scdb_load = function(scdb, obj_type) 
{
	fn = sprintf("%s/%s_%s.RData", scdb@basedir, scdb@fn_prefix, obj_type)
	
	if (!is.element(obj_type, scdb@which_not_lazy) && file.exists(fn)) {
		return(local({ 
			load(fn)
			object
		}))
	}
	else {
		return(NULL)
	}
	
}

#' Loads umi matrix and cells metadata from tab-delimited files
#'
#' @param scdb 
#' @param mat_fn filename of umi table
#' @param metadata_fn filename of metadata information (row per cell)
#' @param obj_type type of scdb object - either raw_mat or clean_mat
#' @param overwrite_if_exists overwrite scdb object if it exists
#'
#' @return tgScMat
#' @export
#'
#' @examples
scdb_load_mat_from_tab = function(scdb, mat_fn, metadata_fn, obj_type, overwrite_if_exists = T) 
{
	message(sprintf("loading sc_mat from %s and %s", mat_fn, metadata_fn))
	mat = fread_rownames(file=mat_fn, sep="\t", set_rownames = T)
	cell_metadata = fread_rownames(file=metadata_fn, sep="\t", set_rownames = T)
	
	sc_mat = tgScMat(mat=as.matrix(mat), stat_type="umi", cell_metadata=cell_metadata)
	
	stopifnot(obj_type == "raw_mat" || obj_type == "clean_mat")
	
	scdb_save(scdb, sc_mat, obj_type, overwrite_if_exists)
	
	sc_mat
}

#' Sets lazy load value for an object type
#'
#' @param scdb 
#' @param obj_type 
#' @param lazy T/F
#'
#' @export
#'
#' @examples
scdb_set_lazy = function(scdb, obj_type, lazy) 
{
	if (lazy) {
		scdb@which_not_lazy = setdiff(scdb@which_not_lazy, obj_type)
	}
	else {
		scdb@which_not_lazy = union(scdb@which_not_lazy, obj_type)
	}
	
	invisible(scdb)
}

#' Get names of objects saved under the scdb instance
#'
#' @param scdb 
#' @param pattern to filter output object names
#'
#' @return vector of object names
#' @export
#'
#' @examples
scdb_get_existing_obj_names = function(scdb, pattern=NULL) 
{
	fns = list.files(path=scdb@basedir, pattern=paste0("^", scdb@fn_prefix), full.names=F)
	obj_nms = gsub(paste0("^", scdb@fn_prefix, "_"), "", gsub(".RData$", "", fns))
	if (!is.null(pattern)) {
		obj_nms = grep(pattern, obj_nms, value=T)
	}
	return(obj_nms)
}