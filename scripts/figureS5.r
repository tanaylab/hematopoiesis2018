dir.create("supp_figures/figureS5")
message("generation Supplementary Fig 5")

tier3_annotations = read.table("annotations/tier3_annotations.txt", sep = "\t", stringsAsFactors = F, header = T)
tier3_annotations = tier3_annotations[!is.infinite(tier3_annotations$threshold),]
markers = tier3_annotations[,1]; tier3_cols = c(tier3_annotations[,2], "darkgray"); thresh = tier3_annotations[,3]
sc_clean_mat = sc_pipe_mat_to_clean_mat(scdb_init(basedir = "saved_work/core_clusts"))
zscore = score_on_gene_prograns(as.matrix(sc_clean_mat@mat), gene_map)
colnames(zscore) = markers

thresh = c(10,13,13,30,25, 15,12,100,20,22, 20,20,10,20,10)

dir.create("supp_figures/figureS5/lin_depletion")
linsx = c(1,2,3,4,5,14,12,9); linsy = c(7,11,11,6,10,15,13,8)
for(i in 1:length(linsx)) {
  linx = linsx[i]; liny = linsy[i]
  png(paste0("supp_figures/figureS5/lin_depletion/", markers[linx], "_", markers[liny], ".png"), height = 1000, width = 1000)
  plot(zscore[,linx] - min(zscore[,linx]) + 1, zscore[,liny] - min(zscore[,liny]) + 1,
        pch = 20, col = "chocolate4", log = "xy", cex = 2, axes = F,
	xlab = "", ylab = "", xlim = c(1, 3e2), ylim = c(1,3e2))
  axis(side = 1, at = 10^(0:2)); axis(side = 4, , at = 10^(0:2))
  abline(v = thresh[linx] - min(zscore[,linx]) + 1, h = thresh[liny] - min(zscore[,liny]) + 1, col = "red", lwd = 3)
  dev.off()
}

#################

sc_2d = sc_pipe_plots(scdb_init(basedir = "saved_work/tier3_clusts"))
sc_cl = sc_2d@scl

tier3_annotations = read.table("annotations/tier3_annotations.txt", sep = "\t", stringsAsFactors = F, header = T)
tier3_annotations = tier3_annotations[!is.infinite(tier3_annotations$threshold),]
markers = tier3_annotations[,1]; tier3_cols = c(tier3_annotations[,2], "darkgray"); thresh = tier3_annotations[,3]
above_exp = sc_cl@clust_fp[ markers, ] > thresh
clust_ass = max.col(t(above_exp), ties.method = "first"); clust_ass[colSums(above_exp) == 0] = length(markers) + 1
names(clust_ass) = colnames(sc_cl@clust_fp)

meta_order = c(16,2,11,3,14,15,13,12,4,6,10,5,7,8,1,9)
clust_ord = order(factor(clust_ass, levels = meta_order))

png("supp_figures/figureS5/lindep_markers_tracks.png", height = 2000, width = 2000)
par(mfrow = c(3,3), mar = c(0,4,1,0))
invisible(sapply(c("Tfrc", "Cd24a", "Fcgr3", "Itgax", "Cd22", "H2-Aa", "Siglech", "Cd19", "Klrb1c;Klrb1b"), function(x)
        barplot(sc_cl@clust_fp[x, clust_ord], col = tier3_cols[clust_ass[ clust_ord]], names.arg = rep("", length(clust_ass)),
        border = NA, ylim = c(1, max(sc_cl@clust_fp[x, ])))))
dev.off()

#################

scdb = scdb_init(basedir="saved_work/core_clusts")
sc_2d = sc_pipe_plots(scdb)
sc_cl = sc_2d@scl
lfp = log2(sc_cl@clust_fp)
cell_stats = sc_cl@scmat@cell_metadata

bad_marks = c("Ccl5", "Vpreb3", "Vpreb1", "Camp")
outliers = names(which(apply(lfp[bad_marks,], 2, max) > 2))
good_clusts = setdiff(colnames(lfp), outliers)
good_cells = names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% good_clusts))
markers = c("Ccl5", "Pf4", "Siglech",  "Cd74", "Vpreb1", "Car1","Ly86","Lmo4","Apoe", "Dntt", "Mpo")
tier3_cols = c("darkorchid4", "darksalmon", "darkcyan", "cyan3", "navyblue", "indianred3", "limegreen",
        "#aa790b", "pink3", "dodgerblue3", "palegreen2", "darkgray")
thresh = rep(4,length(markers)); names(thresh) = markers
above_exp = sc_cl@clust_fp[ markers, ] > thresh
clust_ass = max.col(t(above_exp), ties.method = "first"); clust_ass[colSums(above_exp) == 0] = length(markers) + 1
names(clust_ass) = colnames(sc_cl@clust_fp)

chc = hclust(dist(cor(lfp[c("Hlf", "Flt3", "Irf8", "Gata2", "Apoe", "Ifitm1"),good_clusts])), "ward.D2")
chc = as.hclust(reorder(as.dendrogram(chc),
                             lfp["Apoe",] - lfp["Irf8",],
                             agglo.FUN=max))

meta_order = c(8,6,9,2,12,10,11,7,3)
clust_ord = good_clusts[order( factor(clust_ass[good_clusts], levels = meta_order), chc$order)]

nms = choose_genes_from_clust(sc_cl, clust_ord, 5, 2, bad_genes = bad_genes, max_num = 40, ord = "max.col")
png("supp_figures/figureS5/core_expression.png", height=1500, width = 1500)
par(mar = rep(0.5,4), fig = c(0,1,0.05,1))
cell_ord = plot_sc_heatmap(sc_cl, nms, clust_ord)
par(mar = rep(0.5,4), fig = c(0,1,0,0.05), new = T)
image(matrix(matrix(clust_ass[ sc_cl@clusts[cell_ord]])), axes = F, col = tier3_cols, zlim = c(1, max(clust_ass)))
dev.off()
write.table(rev(nms), row.names = F, quote = F, file = "supp_figures/figureS5/core_expression.txt")

load("saved_work/core_clusts/scrdb_data_cocluster.Rda")
dimnames(cc) = list(names(sc_cl@clusts), names(sc_cl@clusts))
png("supp_figures/figureS5/core_co_clustering.png", height = 1500, width = 1500)
par(mar = rep(0.5,4), fig = c(0.05,1,0.05,1))
image(log(1 + cc[cell_ord, cell_ord]),
        col = colorRampPalette(c("white", "green", "orange2", "red4", "black"))(1000), axes = F)
cls = cumsum(table(factor( sc_cl@clusts[good_cells], levels = clust_ord))) / length(good_cells)
abline(h = cls, v = cls)
par(mar = rep(0.5,4), fig = c(0,0.05,0.05,1), new = T)
image(t(matrix(clust_ass[ sc_cl@clusts[cell_ord]])), axes = F, col = tier3_cols, zlim = c(1, max(clust_ass)))
par(mar = rep(0.5,4), fig = c(0.05,1,0,0.05), new = T)
image(matrix(clust_ass[ sc_cl@clusts[cell_ord]]), axes = F, col = tier3_cols, zlim = c(1, max(clust_ass)))
dev.off()
zlim = quantile(log(1 + cc), c(0,1))
message("coclustering zlim: ", zlim[1], " - ", zlim[2])
