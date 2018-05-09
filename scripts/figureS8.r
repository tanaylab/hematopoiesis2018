dir.create("supp_figures/figureS8")
message("generation Supplementary Fig 8")

pu1_col = "#9d1a7d"

sc_2d = sc_pipe_plots(scdb_init(basedir="saved_work/pu1_clusts"))
sc_cl = sc_2d@scl
nms = choose_genes_from_clust(sc_cl, nms_per_clust=10, nms_thresh=3, max_num=50, bad_genes = bad_genes, ord = "max.col")
cell_stats = sc_cl@scmat@cell_metadata
comb = with(cell_stats, paste0(treatment, ".", tier)); names(comb) = rownames(cell_stats)

message(length(sc_cl@clusts), " cells")
png("supp_figures/figureS8/FigS8C.png", height = 1200, width = 1000)
par(fig = c(0,1,0.1,1), mar = c(1,0,0,0))
cell_ord = plot_sc_heatmap(sc_cl, nms)
par(fig = c(0,1,0,0.1), new = T, mar = rep(0,4))
sample_dist = table(comb, sc_cl@clusts)
X = sample_dist / rowSums(sample_dist); X = sweep(X,2,colSums(X),"/")
image(t(X[, as.character(sc_cl@clusts[ cell_ord])]), col = colorRampPalette(c("white", "gold3", "red4"))(1000), axes = F)
dev.off()
write.table(rev(nms), row.names = F, quote = F, col.names = F, file = "supp_figures/figureS8/FigS8C.txt")

####################

pu1_cl = sc_cl;
pu1_umis = as.matrix(pu1_cl@scmat@mat)
pu1_n = sweep(pu1_umis,2,colSums(pu1_umis),"/") * 1000

bad_marks = c("Vpreb1", "Car2", "Hba-a2", "Beta-s", "Prss34", "Gstm1", "Cd74")
bad_clusts = names(which(apply( sc_cl@clust_fp[ bad_marks,],2,max) > 4))
good_clusts = setdiff(colnames(sc_cl@clust_fp), bad_clusts)
good_cells = names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% good_clusts))

lfp = log2(sc_cl@clust_fp)
wt_cells = rownames(cell_stats)[ cell_stats$treatment == "CRISPR_CTL"]
ga = "Ltf"; gb = "Ccl6"; a = lfp[ga,]; b = lfp[gb,]
g1 = intersect(wt_cells, names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% names(which(b > 2)))))
g2 = setdiff(intersect(wt_cells, names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% good_clusts))), g1)
x = rowSums(pu1_n[,g2]) / length(g2) * min(length(g1), length(g2))
y = rowSums(pu1_n[,g1]) / length(g1) * min(length(g1), length(g2))
genes = setdiff(scr_chi_square_diff_genes(pu1_umis, g1 = g1,g2 = g2, pval = 1e-3, fdr = T), c(cc_genes, ribo_genes))
z = (y + 10) / (x + 10)

mature_genes = names(which(log2(z[genes]) > 1.5))

png("supp_figures/figureS8/FigS8D.png", height = 1000, width = 700)
barplot(rev(head(sort(log2(z[mature_genes]),T),20)), horiz = T, col = "black", las = 2, names.arg = rep("",20))
dev.off()
write.table(names(head(sort(log2(z[mature_genes]),T),20)), row.names = F, quote = F, file = "supp_figures/figureS8/FigS8D.txt")

####################

sc_2d = sc_pipe_plots(scdb_init(basedir="saved_work/iv_clusts"))
sc_cl = sc_2d@scl
nms = choose_genes_from_clust(sc_cl, nms_per_clust=10, nms_thresh=3, max_num=50, bad_genes = bad_genes, ord = "max.col")
cell_stats = sc_cl@scmat@cell_metadata
comb = with(cell_stats, paste0(treatment, ".", cytokines, ".", tier)); names(comb) = rownames(cell_stats)

message(length(sc_cl@clusts), " cells")
png("supp_figures/figureS8/FigS8E.png", height = 1200, width = 1000)
par(fig = c(0,1,0.1,1), mar = c(1,0,0,0))
cell_ord = plot_sc_heatmap(sc_cl, nms)
par(fig = c(0,1,0,0.1), new = T, mar = rep(0,4))
sample_dist = table(comb, sc_cl@clusts)
X = sample_dist / rowSums(sample_dist); X = sweep(X,2,colSums(X),"/")
image(t(X[, as.character(sc_cl@clusts[ cell_ord])]), col = colorRampPalette(c("white", "gold3", "red4"))(1000), axes = F)
dev.off()
write.table(rev(nms), row.names	= F, quote = F,	col.names = F, file = "supp_figures/figureS8/FigS8E.txt")

###############

lfp = log2(sc_cl@clust_fp)
mon_clusts = names(which(lfp["S100a4",] > 3))
mon_cells = names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% mon_clusts))
X = table(comb, names(comb) %in% mon_cells)
X = (X / rowSums(X))[c(2,4),2]
png("supp_figures/figureS8/FigS8G.png", height = 1000, width = 700)
par(lwd=8)
barplot(X, col = c("gray60", pu1_col))
dev.off()

umis = as.matrix(sc_cl@scmat@mat)
umis_n = sweep(umis,2,colSums(umis),"/") * 1000
mon_genes = intersect(rownames(umis), mon_genes)
g1 = intersect(mon_cells, names(which(comb == "CRISPR_PU1.GM-CSF.d9")))
g2 = intersect(mon_cells, names(which(comb == "CRISPR_CTL.GM-CSF.d9")))
x = rowSums(umis_n[,g2]) / length(g2) * min(length(g1), length(g2))
y = rowSums(umis_n[,g1]) / length(g1) * min(length(g1), length(g2))
lim = log2(c(10, max(c(x[mon_genes], y[mon_genes]) + 10)))
png("supp_figures/figureS8/FigS8H.png", height = 1000, width = 1000)
plot(log2(x[mon_genes] + 10), log2(y[mon_genes] + 10), pch = 20, col = "navyblue", axes = F, xlab = "", ylab = "", xlim = lim, ylim = lim, cex = 3)
abline(coef = c(0,1), col = "red")
axis(1); axis(2)
dev.off()

png("supp_figures/figureS8/FigS8H_text.png", height = 1000, width = 1000)
plot(log2(x[mon_genes] + 10), log2(y[mon_genes] + 10), pch = 20, col = "navyblue", axes = F, xlab = "", ylab = "")
abline(coef = c(0,1), col = "red")
text(log2(x[mon_genes] + 10), log2(y[mon_genes] + 10), mon_genes, col = "red")
axis(1); axis(2)
dev.off()
