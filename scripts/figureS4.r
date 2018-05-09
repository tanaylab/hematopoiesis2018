dir.create("supp_figures/figureS4")
message("generation Supplementary Fig 4")


sc_2d = sc_pipe_plots(scdb_init(basedir = "saved_work/tier3_clusts"))
sc_cl = sc_2d@scl
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

sc_clean_mat = sc_pipe_mat_to_clean_mat(scdb_init(basedir = "saved_work/core_clusts"))
zscore = score_on_gene_prograns(as.matrix(sc_clean_mat@mat), gene_map)
colnames(zscore) = markers

thresh = c(10,13,13,30,25, 15,12,100,20,22, 20,20,10,20,10)

dir.create("supp_figures/figureS4/FigS4A")
linsx = c(1,2,3,4,5,14,12,9); linsy = c(7,11,11,6,10,15,13,8)
for(i in 1:length(linsx)) {
  linx = linsx[i]; liny = linsy[i]
  png(paste0("supp_figures/figureS4/FigS4A/", markers[linx], "_", markers[liny], ".png"), height = 1000, width = 1000)
  plot(zscore[,linx] - min(zscore[,linx]) + 1, zscore[,liny] - min(zscore[,liny]) + 1,
	pch = 20, col = "chocolate4", log = "xy", cex = 2, axes = F,
	xlab = "", ylab = "", xlim = c(1, 3e2), ylim = c(1,3e2))
  axis(side = 1, at = 10^(0:2)); axis(side = 4, , at = 10^(0:2))
  abline(v = thresh[linx] - min(zscore[,linx]) + 1, h = thresh[liny] - min(zscore[,liny]) + 1, col = "red", lwd = 3)
  dev.off()
}

#######################

scdb = scdb_init(basedir="saved_work/core_clusts")
sc_2d = sc_pipe_plots(scdb)
sc_cl = sc_2d@scl
lfp = log2(sc_cl@clust_fp)
cell_stats = sc_cl@scmat@cell_metadata
umis = as.matrix(sc_cl@scmat@mat)
umis_n = sweep(umis,2,colSums(umis),"/") * 1000
foc = log(1 + 7 * umis_n)

######################


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
png("supp_figures/figureS4/FigS4B.png", height=1500, width = 1500)
par(mar = rep(0.5,4), fig = c(0,1,0.05,1))
cell_ord = plot_sc_heatmap(sc_cl, nms, clust_ord)
par(mar = rep(0.5,4), fig = c(0,1,0,0.05), new = T)
image(matrix(matrix(clust_ass[ sc_cl@clusts[cell_ord]])), axes = F, col = tier3_cols, zlim = c(1, max(clust_ass)))
dev.off()
write.table(rev(nms), row.names = F, quote = F, file = "supp_figures/figureS4/FigS4B.txt")

load("saved_work/core_clusts/scrdb_data_cocluster.Rda")
dimnames(cc) = list(names(sc_cl@clusts), names(sc_cl@clusts))
png("supp_figures/figureS4/FigS4C.png", height = 1500, width = 1500)
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

#############################

sample_dist = table(cell_stats$tier == 7, sc_cl@clusts)
dist_n = sample_dist / rowSums(sample_dist)
rel_clusts = names(which(sample_dist[2,] >= 10))
rel_clusts = rel_clusts[ order(-dist_n[2, rel_clusts])]
rel_cells = rownames(cell_stats)[cell_stats$tier == 7 & sc_cl@clusts %in% rel_clusts]
s_clust = head(rel_clusts,1)

genes = c("Hlf", "Ifitm1", "Ifitm3", "Serpina3f;Serpina3g", "Apoe", "Gata2","Pf4")
m = t(apply(umis_n[genes, rel_cells], 1, tapply, sc_cl@clusts[rel_cells], sum))

sizes = table(sc_cl@clusts[rel_cells])
m = sweep(m,2,as.vector(sizes),"/") * min(sizes)
IM = log2(sc_cl@clust_fp[genes,rel_clusts])
zlim = max(abs(IM))

png("supp_figures/figureS4/FigS4D.png", height = 1000, width = 1200)
par(fig=c(0,1,0,0.7), mar = c(0.5,3,0.5,0.5))
image(t(IM), axes = F, col = colorRampPalette(c("blue", "white", "red"))(1000), zlim = c(-zlim, zlim))
par(fig = c(0,1,0.7,1), mar = c(0.5,3,0.5,0.5), new = T, lwd = 6)
barplot(dist_n[2, rel_clusts], names.arg = rep("", length(rel_clusts)), xaxs = "i", col = "gray40")
dev.off()
write.table(rev(gsub(";.*", "", genes)), row.names = F, col.names = F, quote = F, file = "supp_figures/figureS4/FigS4D.txt")
message("S clusters zlim: ", zlim)

png("supp_figures/figureS4/FigS4D_cb.png", height = 100, width = 1000)
par(mar = rep(0,4))
image(matrix(1:100), axes = F, col = colorRampPalette(c("blue", "white", "red"))(1000))
dev.off()

tier7_cells = rownames(cell_stats)[cell_stats$tier  == 7]
s_cells = rownames(cell_stats)[cell_stats$tier  == 7 & sc_cl@clusts == s_clust]

g1 = s_cells
g2 = setdiff(good_cells, g1)
bad_genes = c("Atpase6", cc_genes, ribo_genes, rownames(modules[ modules$module == 1,]), "Hba-a2", "Hbb-b1", "Beta-s", "Malat1")
good_genes = setdiff(rownames(umis), bad_genes)
x = rowSums(umis_n[good_genes, g2]) / length(g2) * min(length(g1), length(g2))
y = rowSums(umis_n[good_genes, g1]) / length(g1) * min(length(g1), length(g2))
z = (y + 10) / (x + 10)
disp_genes = names(which(abs(log2(z)) > 1))

png("supp_figures/figureS4/FigS4E.png", height = 1000, width = 1000)
xy_scatter_genes(x,y,bad_genes, disp_genes)
dev.off()

png("supp_figures/figureS4/FigS4E_text.png", height = 1000, width = 1000)
xy_scatter_genes(x,y,bad_genes, disp_genes, cols = c(NA, NA), text = T)
dev.off()

#################

s_genes = setdiff(read.table("results/s_genes.txt", stringsAsFactors=F, header = T)[[1]], ribo_genes)
map = ifelse(rownames(foc) %in% s_genes, "S", ifelse(rownames(foc) %in% cc_genes, "P", "none"))
val_mat = apply(foc, 2, tapply, map, sum)

png("supp_figures/figureS4/FigS4F.png", height = 1500, width = 1500)
shades = colorRampPalette(c("white","white", "white", "wheat1", "orange", "red3", "purple4", "blue"))
s_score = val_mat["S",]; p_score = val_mat["P",]
smoothScatter(s_score, p_score, log = "", colramp = shades, axes = F, xaxs = "i", xlab = "", ylab = "")
points(s_score, p_score, pch = 20)
abline(v = quantile(s_score, c(0.6,0.8,0.95,1)), lwd = 8, lty = 2)
axis(2); axis(1)
dev.off()

################

shades = colorRampPalette(c("blue", "White", "red"))(1000)
all_tfs = read.table("annotations/mouse.tf.list.symbols.txt", stringsAsFactors=F)[[1]]
all_tfs = intersect(rownames(lfp), all_tfs)
C = cor(t(lfp), t(lfp[c("Apoe", "Flt3"),]))
tfs = unique(c("Zfpm1", "Gfi1b", "Gata1", "Gata2", "Pbx1;Pbx3", "Meis1", "Klf1", "Hlf", "Lmo2", "Lmo4",
	"Sfpi1", "Cebpa", "Irf8", "Satb1", "Mta3",
	names(head(sort(C[all_tfs,1],T),25)), names(head(sort(C[all_tfs,2],T),25))))
effectors = unique(c("Car1", "Apoe", "Ifitm1", "Ifitm3", "Gpr56", "Serpina3f;Serpina3g", "Mpo", "Dntt",
	"Flt3", "Csf1r", "Cd34", "Cd52", "Eltd1", "Prtn3", "Fcgr3",
	names(head(sort(C[setdiff(rownames(C), all_tfs),1],T),25)),
	names(head(sort(C[setdiff(rownames(C), all_tfs),2],T),25))))


C1 = cor(t(lfp[tfs,])); diag(C1) = NA
hc = hclust(as.dist(1-C1), "ward.D2")
IM1 = C1[rev(hc$order), hc$order]
C2 = cor(t(lfp[tfs,]), t(lfp[effectors,]))
IM2 = C2[rev(hc$order),]
IM2 = IM2[, rev(order(max.col(t(IM2))))]
zlim = max(abs(cbind(C1,C2)), na.rm = T)
png("supp_figures/figureS4/FigS4J_left.png", height = 1500, width = 1500)
par(mar = rep(0,4))
image(t(IM1), col = shades, zlim = c(-zlim,zlim), axes = F)
dev.off()
png("supp_figures/figureS4/FigS4J_right.png", height = 1500, width = 1500)
par(mar = rep(0,4))
image(t(IM2), col = shades, zlim = c(-zlim,zlim),     axes = F)
dev.off()
write.table(rev(rownames(IM2)), row.names = F, quote = F, col.names = F, file = "supp_figures/figureS4/FigS4J_tfs.txt")
write.table(colnames(IM2), row.names = F, quote = F, col.names = F, file = "supp_figures/figureS4/FigS4J_effectors.txt")
message("TF correlation zlim: ", zlim)
