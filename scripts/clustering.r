
# cluster tiers 1-3
tier1_B = rownames(batch_stats)[ batch_stats$tier == 1]
tier2_B = rownames(batch_stats)[ batch_stats$tier == 2 & batch_stats$noise.estimation <= 0.05]
tier3_B = rownames(batch_stats)[ batch_stats$tier == 3 & batch_stats$treatment == "-"]

scdb = sc_pipeline("saved_work/tier1_clusts", "results/tier1_clusts", tier1_B, index_fn,
	"output/umi.tab/", batch_meta_attr = "amplification.batch", mark_blacklist_terms = bad_genes)
scdb = sc_pipeline("saved_work/tier2_clusts", "results/tier2_clusts", tier2_B, index_fn,
	"output/umi.tab/", batch_meta_attr = "amplification.batch", mark_blacklist_terms = bad_genes)
scdb = sc_pipeline("saved_work/tier3_clusts", "results/tier3_clusts", tier3_B, index_fn,
	"output/umi.tab/", batch_meta_attr = "amplification.batch", mark_blacklist_terms = bad_genes)

perform_bootstrap(scdb, 50, "results/tier3_clusts", mark_blacklist_terms = bad_genes, analysis_only = T)

# Cytokine stimulation
tier3_epo_B = rownames(batch_stats)[ batch_stats$tier == 3 & batch_stats$treatment == "Epo" & batch_stats$noise.estimation <= 0.05]
tier3_gcsf_B = rownames(batch_stats)[ batch_stats$tier == 3 & batch_stats$treatment == "G-CSF" & batch_stats$noise.estimation <= 0.05]
tier3_ectl_B = rownames(batch_stats)[ batch_stats$tier == 3 & batch_stats$treatment == "Epo (control)" & 
	     batch_stats$noise.estimation <= 0.05]
tier3_gctl_B = rownames(batch_stats)[ batch_stats$tier == 3 & batch_stats$treatment == "G-CSF (control)" & 
	     batch_stats$noise.estimation <= 0.05]
tier7_B = setdiff(rownames(batch_stats)[ batch_stats$tier == 7 & batch_stats$treatment %in% c("-", "Epo", "G-CSF")], c("AB3975", "AB3976"))

scdb = sc_pipeline("saved_work/tier3_cytokines", "results/tier3_cytokines", c(tier3_epo_B, tier3_gcsf_B, tier3_ectl_B, tier3_gctl_B), index_fn,
        "output/umi.tab/", batch_meta_attr = "amplification.batch",
        mark_blacklist_terms = bad_genes,
	read_and_clean_only = T)

scdb = sc_pipeline("saved_work/tier7_cytokines", "results/tier7_cytokines", tier7_B, index_fn,
        "output/umi.tab/", batch_meta_attr = "amplification.batch",
        mark_blacklist_terms = bad_genes, read_and_clean_only = T)

sc_2d = sc_pipe_plots(scdb_init(basedir = "saved_work/tier3_clusts"))
sc_cl = sc_2d@scl

# The core and myeloid models (Figures 2 and 4)
tier3_annotations = read.table("annotations/tier3_annotations.txt", sep = "\t", stringsAsFactors = F, header = T)
tier3_annotations = tier3_annotations[!is.infinite(tier3_annotations$threshold),]
markers = tier3_annotations[,1]; tier3_cols = c(tier3_annotations[,2], "darkgray"); thresh = tier3_annotations[,3]
above_exp = sc_cl@clust_fp[ markers, ] > thresh
clust_ass = max.col(t(above_exp), ties.method = "first"); clust_ass[colSums(above_exp) == 0] = length(markers) + 1
names(clust_ass) = colnames(sc_cl@clust_fp)

cluster_genes = scr_create_gene_modules_table(as.matrix(sc_cl@scmat@mat), clust_ass[ sc_cl@clusts],
              clusters = seq_along(markers), K = 50, Z = 2)
cluster_genes[ rowSums(cluster_genes) > 1] = F
colnames(cluster_genes) = markers
gene_map = max.col(cbind(cluster_genes, T), ties.method = "first")
names(gene_map) = rownames(cluster_genes)

core_B = setdiff(rownames(batch_stats)[ batch_stats$tier %in% c(3,5:7) & batch_stats$treatment == "-" & !is.na(batch_stats$noise.estimation) &
       batch_stats$noise.estimation <= 0.05], c("AB3975", "AB3976"))

bad_genes = c("Atpase6", cc_genes, ribo_genes, rownames(modules[ modules$module == 1,]), "Hba-a2", "Hbb-b1", "Beta-s", "Malat1")
scdb = sc_pipeline("saved_work/core_clusts", "results/core_clusts", core_B, index_fn,
        "output/umi.tab/", batch_meta_attr = "amplification.batch", 
	mark_blacklist_terms = bad_genes,
	read_and_clean_only = T)
sc_clean_mat = sc_pipe_mat_to_clean_mat(scdb)
zscore = score_on_gene_prograns(as.matrix(sc_clean_mat@mat), gene_map)
colnames(zscore) = markers

thresh = c(10,13,13,30,25, 15,12,100,20,22, 20,20,10,20,10)
above_thr = t(t(zscore) > thresh)
core_cells = names(which(rowSums(above_thr) == 0))

scdb = sc_pipeline("saved_work/core_clusts", "results/core_clusts", core_B, index_fn,
	"output/umi.tab/", batch_meta_attr = "amplification.batch",
        mark_blacklist_terms = bad_genes, cells = core_cells)
perform_bootstrap(scdb, 50, "results/core_clusts", 
	mark_blacklist_terms = bad_genes)

mye_B = rownames(batch_stats)[ batch_stats$tier %in% 3:5 & batch_stats$treatment %in% c("-")]
bad_genes = c("Atpase6", cc_genes, ribo_genes, "Camp", "Ltf", "Lcn2", "S100a8", "S100a9")
scdb = sc_pipeline("saved_work/myeloid_clusts", "results/myeloid_clusts", mye_B, index_fn,
        "output/umi.tab/", batch_meta_attr = "amplification.batch",
        mark_blacklist_terms = bad_genes,
        read_and_clean_only = T)
sc_clean_mat = sc_pipe_mat_to_clean_mat(scdb)
zscore = score_on_gene_prograns(as.matrix(sc_clean_mat@mat), gene_map)
colnames(zscore) = markers
thresh = c(10,13,13,30,Inf,15,12,100,20,Inf,20,Inf,Inf,Inf,Inf)
above_thr = t(t(zscore) > thresh)
mye_cells = names(which(rowSums(above_thr) == 0))
scdb = sc_pipeline("saved_work/myeloid_clusts", "results/myeloid_clusts", mye_B, index_fn,
        "output/umi.tab/", batch_meta_attr = "amplification.batch",
        mark_blacklist_terms = bad_genes, cells = mye_cells)
perform_bootstrap(scdb, 50, "results/myeloid_clusts",
	mark_blacklist_terms = bad_genes)

sc_2d = sc_pipe_plots(scdb_init(basedir = "saved_work/tier3_clusts"))
#reposition_coords = read.delim("annotations/tier3_reposition_coords.txt", stringsAsFactor=F)
#for (i in seq_len(nrow(reposition_coords))) { sc_2d = reposition_cc(sc_2d, coords = as.numeric(reposition_coords[i,]))$sc_2d}
sc_cl = sc_pipe_plots(scdb_init(basedir = "saved_work/myeloid_clusts"))@scl
dir.create("temp")
map = proj_ds_on_graph(sc_2d, umis = as.matrix(sc_cl@scmat@mat),
            bg_cells = NULL, bw = 30, fn = "temp/proj_mye.png", 
            reg = 10, bg_reg = 3, cex = 3, lwd = 2, markers = intersect(rownames(sc_cl@scmat@mat), rownames(sc_2d@scl@feat_mat)))
map$x_cl = tapply(map$coords$x, sc_cl@clusts, mean)
map$y_cl = tapply(map$coords$y, sc_cl@clusts, mean)

contam = c("Prss34", "Prg2", "Mcpt8", "Vpreb1", "Vpreb3", "Car1", "Mt2", "Klf1", "Ccl5", "Pf4","Apoe","Cd79b")
contam_val = round(pmin(apply(sc_cl@clust_fp[contam,],2,max),4))
good_clusts = names(which(contam_val <= 2))

mye_2d = sc_pipe_plots(scdb_init(basedir = "saved_work/myeloid_clusts"))
mye_cl = sc_cl; 
mye_cells = names(mye_cl@clusts)[ mye_cl@clusts %in% good_clusts]
mye_cl@scmat@mat = mye_cl@scmat@mat[,mye_cells]
mye_cl@clusts = mye_cl@clusts[ mye_cells];
mye_cl@clust_fp = .calc_clusts_fp(mye_cl)
mye_2d@scl = mye_cl
mye_2d@x = map$coords$x[ mye_cells]; mye_2d@y = map$coords$y[mye_cells]
mye_2d@x_cl = map$x_cl[good_clusts]; mye_2d@y_cl = map$y_cl[good_clusts]
sc_2d = mye_2d
save(sc_2d, file = "saved_work/myeloid_clusts/scrdb_data_projected_2d.RDa")

# CRISPR clustering
ugi_B = rownames(batch_stats)[ batch_stats$sequencing.batch.unique.ID %in% c("SB134", "SB136", "SB166", "SB167", "SB170", "SB171") &
        batch_stats$virus.mix != "mix 4"]
scdb = sc_pipeline("saved_work/ugi_clusts", "results/ugi_clusts", ugi_B, index_fn,
        "output/umi.tab/", batch_meta_attr = "amplification.batch", mark_blacklist_terms = bad_genes)

pu1_B = rownames(batch_stats)[ batch_stats$treatment %in% c("CRISPR_PU1", "CRISPR_CTL") &
	batch_stats$tier %in% c(3, "Ly6g+")]
iv_B = rownames(batch_stats)[ batch_stats$treatment %in% c("CRISPR_PU1", "CRISPR_CTL") &
        !(batch_stats$tier %in% c(3, "Ly6g+", "d14"))]

scdb = sc_pipeline("saved_work/pu1_clusts", "results/pu1_clusts", pu1_B, index_fn,
        "output/umi.tab/", batch_meta_attr = "amplification.batch", mark_blacklist_terms = bad_genes)
scdb = sc_pipeline("saved_work/iv_clusts", "results/iv_clusts", iv_B, index_fn,
        "output/umi.tab/", batch_meta_attr = "amplification.batch", mark_blacklist_terms = bad_genes)
