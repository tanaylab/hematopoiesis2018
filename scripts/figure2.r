###############
# Figure 2
##############

message("generating figure 2")
dir.create("figures/figure2")

sc_2d = sc_pipe_plots(scdb_init(basedir = "saved_work/tier3_clusts"))
sc_cl = sc_2d@scl

tier3_annotations = read.table("annotations/tier3_annotations.txt", sep = "\t", stringsAsFactors = F, header = T)
tier3_annotations = tier3_annotations[!is.infinite(tier3_annotations$threshold),]
markers = tier3_annotations[,1]; tier3_cols = c(tier3_annotations[,2], "darkgray"); thresh = tier3_annotations[,3]
above_exp = sc_cl@clust_fp[ markers, ] > thresh
clust_ass = max.col(t(above_exp), ties.method = "first"); clust_ass[colSums(above_exp) == 0] = length(markers) + 1
names(clust_ass) = colnames(sc_cl@clust_fp)

secondary = rep(0, length(clust_ass)); names(secondary) = names(clust_ass)
for (i in 1:length(markers)) {
  clusts = which(clust_ass == i)
  if (length(clusts) < 2) {next}
  marker = markers[i]
  hc = hclust(dist(t(sc_cl@clust_fp[rownames(sc_cl@feat_mat), clusts])), method = "ward.D2")
  d = as.hclust(reorder(as.dendrogram(hc),
              sc_cl@clust_fp[marker, clusts],
              agglo.FUN=mean))
  secondary[clusts] = rank(factor(clusts, levels = clusts[ d$order]))
}

clusts = which(clust_ass == length(markers) + 1)
marker = "Ifitm1"
hc = hclust(dist(t(sc_cl@clust_fp[rownames(sc_cl@feat_mat), clusts])), method = "ward.D2")
d = as.hclust(reorder(as.dendrogram(hc),
                      sc_cl@clust_fp[marker, clusts],
                      agglo.FUN=mean))
secondary[clusts] = rank(factor(clusts, levels = clusts[ rev(d$order)]))

x_cl = sc_2d@x_cl
y_cl = sc_2d@y_cl
filt = names(x_cl)

filt_cells = names(sc_cl@clusts)
scr_umis = as.matrix(sc_cl@scmat@mat)
umis_norm = t(t(scr_umis) / colSums(scr_umis)) * 1000
foc = log2(1 + 7 * umis_norm[,filt_cells])
markers = c('Ifitm1', 'Ifitm3', 'Txnip', 'Pf4', 'Cd34', 'Flt3', 'Gata1', 'Car2', 'Mt2', 'Car1', 'Hba-a2', 'Hbb-b1', 'Fcgr3', 'Ctsg', 'Mpo', 'Ly6c2', 'Ly86', 'F13a1', 'Csf1r',
	'Lyz1;Lyz2', 'Elane', 'Fcnb', 'Gstm1', 'Camp', 'Ltf', 'Cd63', 'Prss34', 'Prg2', 'Cd74', 'H2-Aa', 'Ly6d', 'Siglech', 'Dntt',
	'Vpreb1', 'Vpreb3', 'Cd79b', 'Fcrla', 'Ccl5', 'Nkg7', 'C1qb')
meta_order = c(16,2,11,3,14,15,13,12,4,6,10,5,7,8,1,9)
clust_ord = order(factor(clust_ass, levels = meta_order), secondary)
cls = cumsum(table(factor(clust_ass[ sc_cl@clusts], levels = meta_order))) / length(filt_cells)

png("figures/figure2/Fig2A.png", height = 2000, width = 2000)
par(mar = rep(0,4))
cell_ord = plot_sc_heatmap(sc_cl, rev(markers), clust_ord, draw_cls = F)
abline(v = cls, lwd = 2, col = "gray40", lty = 2)
dev.off()
write.table(markers, row.names = F, col.names = F, quote = F, file = "figures/figure2/Fig2A.txt")

png("figures/figure2/Fig2A_cb.png", width = 1500, height = 500)
par(mar = rep(0,4))
image(matrix(1:130), col = colorRampPalette(c("white", "orange", "tomato","mediumorchid4", "midnightblue"))(130), axes = F)
dev.off()

png("figures/figure2/Fig2A_cluster_bar.png", width = 1500, height = 200)
par(mar = rep(0,4))
image(matrix(as.numeric(factor(clust_ass[sc_cl@clusts[cell_ord]], levels = meta_order))), col = tier3_cols[meta_order], axes = F)
abline(v = cls,lwd=2) # cluster seperators
dev.off()

###############

bad_cells= c()
P = .graph_con_comp(sc_2d@clust_graph); large_comp = which.max(table(P))
good_clusts = names(which(P == large_comp))
good_cells = setdiff( names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% good_clusts)), bad_cells)
outline = scr_find_outline(sc_2d, reg = 0.75, cells = good_cells)

###############

png("figures/figure2/Fig2B.png", width = 1200, height = 1500)
plot(sc_2d@x, sc_2d@y, axes = F, xlab = "", ylab = "",
     bg = tier3_cols[clust_ass[ as.numeric(sc_cl@clusts)]], pch = 21, cex = 2)
lines(outline[,1], outline[,2], lwd = 4)
dev.off()

###############

dir.create("figures/figure2/Fig2C")
con_2d = sc_2d
con_2d@scl@scmat@mat = con_2d@scl@scmat@mat[,good_cells]
con_2d@x = con_2d@x[good_cells]; con_2d@y = con_2d@y[good_cells]
invisible(sapply(c("Pf4", "Car1", "Hba-a2", "Lyz1;Lyz2", "Gstm1", "Prss34", "Siglech", "Cd74", "Ly86", "F13a1", "Csf1r",
      "Fcnb", "Cebpe", "Camp", "Ngp", "Ccl5"), function(x) plot_gene_2d(con_2d, x, 2000, 2000,
      "figures/figure2/Fig2C", reg_factor = 10, pt_shades = colorRampPalette(c("black")), 
      outline = outline, positive_psize = 1, negative_psize = 1)))

rna_shades = colorRampPalette(c("white", "white", "lightgray", "darkorange1", "darkgoldenrod1", "darkgoldenrod4"))(100)
png("figures/figure2/Fig2C/Fig2C_cb.png", height = 100, width =1000)
par(mar = rep(0,4))
image(matrix(1:100), axes = F, col = rna_shades)
dev.off()

###############

tier3_gated_cells = intersect(good_cells, rownames(wells))
dir.create("figures/figure2/Fig2D/")
for (gate in setdiff(unique(tier3_gates$gate), c("other", NA))) {
  png(paste0("figures/figure2/Fig2D/", gate, ".png"), width = 1500, height = 1500)
  plot(sc_2d@x[tier3_gated_cells], sc_2d@y[tier3_gated_cells], pch = 21, 
       bg =  ifelse(tier3_gates[ tier3_gated_cells, "gate"] == gate, "chocolate3", "white"), 
       cex = ifelse(tier3_gates[ tier3_gated_cells, "gate"] == gate, 3, 0),
      axes = F, xlab = "", ylab = "")
  points(outline[,1], outline[,2], type = "l", lwd = 4)
  dev.off()
}

###############

clusts = sc_cl@clusts
filt_cells = names(clusts)
umis = as.matrix(sc_cl@scmat@mat)
cc_umis = colSums(umis[cc_genes, filt_cells]) / colSums(umis[, filt_cells]) * 1000
cc_by_cluster = tapply(cc_umis, clusts, median)
markers = c("Ccl5", "Pf4", "Hba-a2","Prss34", "Siglech", "Prg2", "Vpreb1", "Fcrla", "C1qb", "Cd74", "Car1", "Ltf", "Fcnb", "Ly86", "Gstm1")
metac = c(markers, "core")[ clust_ass]
names(tier3_cols) = c(markers, "core")

dir.create("figures/figure2/Fig2E")
cols = tier3_cols#[-c(1,10)]
metac_u = unique(metac)
ylim = c(0, max(sapply(seq_along(metac_u), function(i) max(density(log10(cc_umis[ metac[clusts] == metac_u[i]]))$y))))
for(i in 1:length(metac_u)) {
  png(paste0("figures/figure2/Fig2E/cc_", metac_u[i], ".png"), width = 1000, height = 1000)
  par(lwd = 8)
  plot(density(log10(cc_umis[ metac[clusts] == metac_u[i]])), col = "black", main = "", axes = F, xlab = "", ylab = "", xlim = c(0,3), ylim = ylim)
  polygon(density(log10(cc_umis[ metac[clusts] == metac_u[i]])), col = tier3_cols[ metac_u[i]], border = "black")
  axis(1, at = c(1,2)); axis(2)
  dev.off()
}
