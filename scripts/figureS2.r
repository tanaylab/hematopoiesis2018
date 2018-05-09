dir.create("supp_figures/figureS2")
message("generation Supplementary Fig 2")

sc_2d = sc_pipe_plots(scdb_init(basedir = "saved_work/tier1_clusts"))
sc_cl = sc_2d@scl

bad_clusts = c(); good_clusts = setdiff(colnames(sc_cl@clust_fp), bad_clusts)
good_cells = names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% good_clusts))
nms = choose_genes_from_clust(sc_cl, good_clusts, 7,3,bad_genes=bad_genes, max_num = 60, ord = "max.col")

cells = names(sc_cl@clusts)
gate_dist = table(sc_cl@clusts[cells], lin_gate$gate[cells]) + 1
gate_dist = log2((gate_dist[,2] / sum(gate_dist[,2])) / (rowSums(gate_dist) / sum(gate_dist)))
gate_dist[ is.infinite(gate_dist)] = NA

png("supp_figures/figureS2/FigS2B.png", height = 2000, width = 3300)
par(fig = c(0,1,0.05,1), mar = c(1,0,0,0))
cell_ord = plot_sc_heatmap(sc_cl, nms, good_clusts)
par(fig = c(0,1,0,0.05), new = T)
zlim = max(abs(gate_dist))
image(matrix(gate_dist[as.character(sc_cl@clusts[cell_ord])]), col = colorRampPalette(c("blue", "gray80", "red"))(1000), axes = F, zlim = c(-zlim,zlim))
dev.off()
write.table(rev(nms), row.names = F, col.names = F, quote = F, file = "supp_figures/figureS2/FigS2B.txt")
message("tier1 enrichment zlim: ", zlim)

sc_2d = sc_pipe_plots(scdb_init(basedir = "saved_work/tier2_clusts"))
sc_cl = sc_2d@scl
bad_clusts = c(); good_clusts = setdiff(colnames(sc_cl@clust_fp), bad_clusts)
good_cells = names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% good_clusts))
nms = choose_genes_from_clust(sc_cl, good_clusts, 7,3,bad_genes=bad_genes, ord = "max.col")

cells = names(sc_cl@clusts)
gate_dist = table(sc_cl@clusts[cells], ckit_gate$gate[cells]) + 1
gate_dist = log2((gate_dist[,2] / sum(gate_dist[,2])) / (rowSums(gate_dist) / sum(gate_dist)))
gate_dist[ is.infinite(gate_dist)] = NA

png("supp_figures/figureS2/FigS2C.png", height = 2000, width = 3300)
par(fig = c(0,1,0.05,1), mar = c(1,0,0,0))
cell_ord = plot_sc_heatmap(sc_cl, nms, good_clusts)
par(fig = c(0,1,0,0.05), new = T)
zlim = max(abs(gate_dist))
image(matrix(gate_dist[as.character(sc_cl@clusts[cell_ord])]), col = colorRampPalette(c("blue", "gray80", "red"))(1000), axes = F, zlim = c(-zlim, zlim))
dev.off()
write.table(rev(nms), row.names = F, col.names = F, quote = F, file = "supp_figures/figureS2/FigS2D.txt")
message("tier2 enrichment zlim: ", zlim)

sc_2d = sc_pipe_plots(scdb_init(basedir = "saved_work/tier3_clusts"))
sc_cl = sc_2d@scl

bad_clusts = c(); good_clusts = setdiff(colnames(sc_cl@clust_fp), bad_clusts)
good_cells = names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% good_clusts))
nms = choose_genes_from_clust(sc_cl, good_clusts, 3,3,bad_genes=bad_genes, must_haves = c("F13a1", "Csf1r", "Fcnb", "C1qb"), ord = "max.col")
cell_ord = plot_sc_heatmap(sc_cl, nms, good_clusts, fname = "supp_figures/figureS2/FigS2D.png")
write.table(rev(nms), row.names = F,    col.names = F, quote = F, file = "supp_figures/figureS2/tier3_heatmap.txt")
