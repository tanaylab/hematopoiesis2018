cells = unique(c(names(sc_pipe_plots(scdb_init(basedir="saved_work/tier1_clusts"))@scl@clusts),
	names(sc_pipe_plots(scdb_init(basedir="saved_work/tier2_clusts"))@scl@clusts),
	names(sc_pipe_plots(scdb_init(basedir="saved_work/tier3_clusts"))@scl@clusts),
	names(sc_pipe_plots(scdb_init(basedir="saved_work/core_clusts"))@scl@clusts),
	names(sc_pipe_plots(scdb_init(basedir="saved_work/myeloid_clusts"))@scl@clusts),
	names(sc_pipe_plots(scdb_init(basedir="saved_work/ugi_clusts"))@scl@clusts),
	names(sc_pipe_plots(scdb_init(basedir="saved_work/pu1_clusts"))@scl@clusts),
	names(sc_pipe_plots(scdb_init(basedir="saved_work/iv_clusts"))@scl@clusts),
	colnames(sc_pipe_mat_to_clean_mat(scdb_init(basedir = "saved_work/tier3_cytokines"))@mat),
	colnames(sc_pipe_mat_to_clean_mat(scdb_init(basedir = "saved_work/tier7_cytokines"))@mat)))

wells_cells = "~/amir/scdb/annotations/wells_cells.txt"
wells = read.delim(wells_cells, header = T, row.names = 1, stringsAsFactors = F)
B = unique(wells[cells, "Amp_batch_ID"])

index_fn = "annotations/batch_stats.txt"
index = read.delim(index_fn, stringsAsFactors = F)

dir.create("saved_work/all_cells/")
sc_raw_mat = sc_pipe_build_mat(scdb_init(basedir="saved_work/all_cells"),
        index_fn = index_fn,
        batch_meta_attr = "amplification.batch",
        base_dir= "output/umi.tab/",
        mars_batches = B,
        outdir = ".",
	min_umi_n = 0)


umis_all = as.matrix(sc_raw_mat@mat)
B = unique(wells[cells, "Amp_batch_ID"])
#all_cells = rownames(wells)[ wells$Amp_batch_ID %in% B & wells$Number_of_cells == 1]
#umis_all = umis_all[,all_cells]
all_cells = colnames(umis_all)

amp_batches = "~/amir/scdb/annotations/amp_batches.txt"
batches = read.delim(amp_batches, stringsAsFactors = F)

read_stats = read.delim("annotations/read_stats.txt", sep = "", stringsAsFactor = F)
read_stats_old = read.delim("annotations/read_stats_old.txt", stringsAsFactor = F)[,-1]
old_wells = setdiff(read_stats_old$well_id, read_stats$well_id)
read_stats = rbind(read_stats, read_stats_old[ read_stats_old$well_id %in% old_wells,])
read_stats = read_stats[ read_stats$well_id %in% all_cells, c("well_id", "total")]
read_stats = unique(read_stats)
X = read_stats[ read_stats$well_id %in% names(which(table(read_stats$well_id) > 1)),]
X = X[ order(X$well_id),]
X = X[seq(1,nrow(X),2),]
read_stats = read_stats[ read_stats$well_id %in% names(which(table(read_stats$well_id) == 1)),]
read_stats = rbind(read_stats,X)
rownames(read_stats) = read_stats$well_id

project_cells = wells[ wells$Amp_batch_ID %in% B, c("well_coordinates", "plate_ID", "Subject_ID", "Cell_barcode", "Number_of_cells", "Amp_batch_ID")]
project_cells$well = rownames(project_cells)
project_cells = merge(project_cells, batches[,c("Amp_batch_ID", "Seq_batch_ID", "Pool_barcode", "Experiment_ID")], by.x = "Amp_batch_ID", by.y = "Amp_batch_ID")
rownames(project_cells) = project_cells$well
project_cells$filtered_in = 1 * (project_cells$well %in% cells)
project_cells$umicount = 0
project_cells[all_cells, "umicount"] = colSums(umis_all[,all_cells])

project_cells$readcount = 0
project_cells[all_cells, "readcount"] = read_stats[all_cells, "total"]

batches_stats = read.delim("annotations//amp_batches_stats.txt")
batches_stats_new = read.delim("annotations/amp_batches_stats_new.txt")
batches_stats_3 = read.delim("annotations/amp_batches_stats_3.txt")
batches_stats = rbind(batches_stats_3, batches_stats_new, batches_stats)
project_batches_stats = batches_stats[unlist(lapply(sapply(paste0("^", B, "$"), grep, batches_stats$amp_batch), "[[", 1)),]
project_batches_stats = merge(project_batches_stats, unique(batches[,c("Amp_batch_ID", "Experiment_ID", "Seq_batch_ID")]), by.x = "amp_batch", by.y = "Amp_batch_ID")
project_batches_stats = merge(project_batches_stats, unique(wells[,c("Amp_batch_ID", "Subject_ID")]), by.x = "amp_batch", by.y = "Amp_batch_ID")
project_batches_stats = merge(project_batches_stats, tapply(project_cells$filtered_in, factor(project_cells$Amp_batch_ID), sum), by.x = "amp_batch", by.y = 0)
colnames(project_batches_stats)[ncol(project_batches_stats)] = "ncells"
project_batches_stats = project_batches_stats[,c(1,17,18,19,20,2,6,8,12)]

write.table(project_batches_stats, sep = "\t", col.names = NA, file = "results/TableS1.txt")
write.table(project_cells, sep = "\t", col.names = NA, file = "results/project_cells.txt")
