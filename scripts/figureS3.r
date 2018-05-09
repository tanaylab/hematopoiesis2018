dir.create("supp_figures/figureS3")
message("generation Supplementary Fig 3")

sc_2d = sc_pipe_plots(scdb_init(basedir = "saved_work/tier3_clusts"))
sc_cl = sc_2d@scl
umis = as.matrix(sc_cl@scmat@mat)

gstat = umi_gene_stat(sc_cl@scmat)
marks = rownames(sc_cl@feat_mat)
umis_n = t(t(umis)/colSums(umis))
score = gstat$niche_norm[marks] - gstat$sz_cor_norm[marks]
marks_to_clust = tail(marks[order(score)],n=800)

png("supp_figures/figureS3/FigS3A.png", height = 1000, width = 1000)
  plot(log2(gstat$tot), gstat$niche_stat, cex=1.5, pch=19,
       xlab ="log2 total expression", ylab="niche score")
  points(log2(gstat[marks,"tot"]), gstat[marks,"niche_stat"], cex=2, col="red", pch=19)
  grid()
dev.off()

png("supp_figures/figureS3/FigS3B.png", height = 1000, width = 1000)
  plot(log2(gstat$tot), gstat$sz_cor, cex=1.5, pch=19,
       xlab ="log2 total expression", ylab="size correlation")
  points(log2(gstat[marks,"tot"]), gstat[marks,"sz_cor"], cex=2, col="red", pch=19)
  grid()
dev.off()

png("supp_figures/figureS3/FigS3C.png", height = 1000, width = 1000)
  f=gstat$ds_log_varmean != 0
  plot(log2(gstat$ds_mean[f]), gstat$ds_log_varmean[f], cex=1.5, pch=19,
       xlab = "log2(mean down-sampled)", ylab="log2(var/mean down-sampled)")
  points(log2(gstat[marks,"ds_mean"]), gstat[marks,"ds_log_varmean"],
         cex=2, col="red", pch=19)
  grid()
dev.off()

tier3_annotations = read.table("annotations/tier3_annotations.txt", sep = "\t", stringsAsFactors = F, header = T)[-16,]
markers = tier3_annotations[,1]; tier3_cols = c(tier3_annotations[,2], "darkgray"); thresh = tier3_annotations[,3]
above_exp = sc_cl@clust_fp[ markers, ] > thresh
clust_ass = max.col(t(above_exp), ties.method = "first"); clust_ass[colSums(above_exp) == 0] = length(markers) + 1
names(clust_ass) = colnames(sc_cl@clust_fp)

meta_order = c(16,2,11,3,14,15,13,12,4,6,10,5,7,8,1,9)
clust_ord = order(factor(clust_ass, levels = meta_order))

coc_shades = colorRampPalette(c("white", "green", "orange2", "red4", "black"))(1000)
load("saved_work/tier3_clusts/scrdb_data_cocluster.Rda")
dimnames(cc) = list(names(sc_cl@clusts), names(sc_cl@clusts))
IM = log(1 + cc[order(factor(sc_cl@clusts, levels = clust_ord)), order(factor(sc_cl@clusts, levels = clust_ord))])
png("supp_figures/figureS3/FigS3D.png", height = 1500, width = 1500)
par(mar = rep(0.5,4), fig = c(0.05,1,0.05,1))
image(IM, col = coc_shades, axes = F)
cls = cumsum(table(factor( clust_ass[sc_cl@clusts], levels = meta_order))) / length(sc_cl@clusts)
abline(h = cls, v = cls, lty = 2, lwd = 2, col = "gray40")
par(mar = rep(0.5,4), fig = c(0,0.05,0.05,1), new = T)
image(t(matrix(clust_ass[ sc_cl@clusts[colnames(IM)]])), axes = F, col = tier3_cols, zlim = c(1, max(clust_ass)))
par(mar = rep(0.5,4), fig = c(0.05,1,0,0.05), new = T)
image(matrix(clust_ass[ sc_cl@clusts[colnames(IM)]]), axes = F, col = tier3_cols, zlim = c(1, max(clust_ass)))
dev.off()
zlim = quantile(IM, c(0,1))
message("coclustering zlim: ", zlim[1], " - ", zlim[2])

png("supp_figures/figureS3/FigS3D_cb.png", height = 100, width = 1000)
par(mar = rep(0,4))
image(matrix(1:100), axes = F, col = coc_shades)
quantile(IM, c(0,1))
dev.off()

png("supp_figures/figureS3/FigS3F.png", height = 1000, width = 2500)
par(mfrow = c(3,5), mar = c(0,6,1,0))
invisible(sapply(meta_order[-1], function(x)
        barplot(sc_cl@clust_fp[markers[x], clust_ord], col = tier3_cols[clust_ass[ clust_ord]], names.arg = rep("", length(clust_ass)), border = NA)))
dev.off()

marker.exp = as.data.frame(t(umis[ markers, ]))
expressed = marker.exp > 2
X = expressed[,meta_order[-1]]
png("supp_figures/figureS3/FigS3G.png", width = 2000, height = 2000)
par(mar = rep(0,4))
image(expressed[ order(max.col(cbind(0,X), ties.method = "first")), rev(meta_order[-1])], col = c("white", "black"), axes = F)
dev.off()

######################

batch2sort = read.delim("annotations/batch2sort.txt", stringsAsFactor = F)
cell_stats = sc_cl@scmat@cell_metadata
cell_stats$well= rownames(cell_stats)
cell_stats$map = sc_cl@clusts
#cell_stats = merge(cell_stats, batch2sort[,c(1,5,6)], by.x = "amplification.batch", by.y = "batch")

cell_stats = with(cell_stats, cell_stats[order(mouse.specimen.unique.ID, sort.batch, sort.plate),])
rownames(cell_stats) = cell_stats$well
tier3_batches = with(cell_stats, unique(as.vector(amplification.batch[ order(mouse.specimen.unique.ID)])))
batch_dist = with(cell_stats, table(factor(amplification.batch), map))
clust_factor = factor(clust_ass, levels = meta_order)
cols = tier3_cols[meta_order]
batch_dist = batch_dist[ tier3_batches, order(clust_factor)]

png("supp_figures/figureS3/FigS3E.png", width = 2000, height = 2000)
par(mar = rep(0,4))
barplot(t(batch_dist / rowSums(batch_dist)), las =2 , col = cols[ sort(clust_factor)], names.arg = rep("", nrow(batch_dist)), xaxs = "i")
dev.off()

png("supp_figures/figureS3/FigS3E_cb.png", width = 1000, height = 200)
par(mar = c(1,rep(0,3)), mfrow = c(2,1))
#par(mar = rep(0,4))
batch2subject = unique(cell_stats[,c("amplification.batch" ,"mouse.specimen.unique.ID", "sort.batch", "sort.plate")])
with(batch2subject, image(matrix(as.numeric(factor(mouse.specimen.unique.ID))), col = rainbow(length(unique(mouse.specimen.unique.ID)))))
with(batch2subject, image(matrix(as.numeric(factor(sort.batch))), col = rainbow(length(unique(sort.batch)))))
dev.off()

#######################

igh_nm = grep("Igh", rownames(umis), value = T)[1]
pairs_x = c("Flt3", "Vpreb1", "Vpreb1", "Flt3", "Flt3", "Flt3", "Vpreb1")
pairs_y = c("Dntt", "Vpreb3", "Igll1", "Vpreb1", "Igll1", igh_nm, "Dntt")
genes = unique(c(pairs_x, pairs_y))

fp = t(apply(umis[genes,], 1, tapply,  sc_cl@clusts, sum)) /
		as.vector(tapply(colSums(umis), sc_cl@clusts, sum)) * 1000

dir.create("supp_figures/figureS3/Fig3H")
good_clusts = colnames(sc_cl@clust_fp)
bad_clusts = c(); good_clusts = setdiff(colnames(sc_cl@clust_fp), bad_clusts)
fp = sc_cl@clust_fp[genes,]
for (i in seq_along(pairs_x)) {
	x = pairs_x[i]; y = pairs_y[i]
	png(paste0("supp_figures/figureS3/Fig3H/", x, "-", gsub(";.*", "", y), ".png"), height = 1000, width = 1000)
	plot(fp[x, good_clusts], fp[y, good_clusts], cex = 5, pch = 21, bg = tier3_cols[ clust_ass[ good_clusts]],
		lwd = 3, axes = F, xlab = "", ylab = "")
	axis(1); axis(2)
	dev.off()
}

####################

dir.create("supp_figures/figureS3/Fig3I")
plot_virtual_facs(wells, "Sca1", "cKit", "supp_figures/figureS3/Fig3I/all.png", gates = list(lsk$polygon, lk$polygon))
plot_virtual_facs(wells, "CD34", "Flt3", "supp_figures/figureS3/Fig3I/lsk.png", filter = lsk$gate, gates = list(lt$polygon, st$polygon, mpp$polygon))
plot_virtual_facs(wells, "CD34", "FcgR", "supp_figures/figureS3/Fig3I/lk.png", filter = lk$gate, gates = list(cmp$polygon, mep$polygon, gmp$polygon))
plot_virtual_facs(wells, "IL7Ra", "cKit", "supp_figures/figureS3/Fig3I/clp.png", gates = list(clp$polygon))

###################


P = .graph_con_comp(sc_2d@clust_graph); large_comp = which.max(table(P))
regs = rep(0.97, 5);
regs[c(1,2)] = 0.75;
outlines = sapply(seq_len(max(P)), function(x) scr_find_outline(sc_2d, reg = regs[x], cells = names(sc_cl@clusts)[ P[ sc_cl@clusts] == x]))
tier3_gated_cells = names(sc_cl@clusts)
dir.create("supp_figures/figureS3/Fig3J/")
for (gate in setdiff(unique(tier3_gates$gate), c("other", NA))) {
  png(paste0("supp_figures/figureS3/Fig3J/", gate, ".png"), width = 1500, height = 1500)
  plot(sc_2d@x[tier3_gated_cells], sc_2d@y[tier3_gated_cells], pch = 21,
       bg =  ifelse(tier3_gates[ tier3_gated_cells, "gate"] == gate, "chocolate3", "white"),
       cex = ifelse(tier3_gates[ tier3_gated_cells, "gate"] == gate, 3, 0),
      axes = F, xlab = "", ylab = "")
  sapply(seq_len(max(P)), function(x) lines(outlines[[ x * 2 - 1]], outlines[[x * 2]], lwd = 4))
  dev.off()
}
