dir.create("supp_figures/figureS6")
message("generation Supplementary Fig 6")

load("saved_work/myeloid_clusts/scrdb_data_projected_2d.RDa")
sc_cl = sc_2d@scl
mye_umis = as.matrix(sc_cl@scmat@mat)
mye_n = sweep(mye_umis,2,colSums(mye_umis), "/") * 1000
mye_foc = log(1 + 7 * mye_n)

sc_clean_mat = sc_pipe_mat_to_clean_mat(scdb_init(basedir = "saved_work/myeloid_clusts"))
all_umis = as.matrix(sc_clean_mat@mat)

lin_ord = c(2,9,3,12,13,11,10,4,7,6,5,8,1,14,15)
markers = c("Ccl5", "Pf4", "Hba-a2","Prss34", "Siglech",  "Cd74", "Prg2", "Vpreb1", "Car1", "Ltf", "Fcnb", "Ly86", "Gstm1", "Myl4", "Fcrla")[lin_ord]
mature_cols = c("darkorchid4", "darksalmon", "indianred4", "goldenrod1",
           "darkcyan", "cyan3", "goldenrod3", "navyblue", "indianred3", "darkolivegreen4", "darkgreen", "limegreen", "green4", "dodgerblue3", "darkslateblue")[lin_ord]

all.exp = as.data.frame(t(all_umis[ markers, ]))
wmax = max.col(all.exp, ties.method = "first")
expressed = rowSums(all.exp > 1)
identity = ifelse(expressed, markers[wmax], "none")
all_i = table(sc_clean_mat@cell_metadata$tier, factor(identity, levels = c("none", markers)))
mye_cells = colnames(mye_umis)
mye_i = table(sc_clean_mat@cell_metadata[mye_cells, "tier"], factor(identity[mye_cells], levels = c("none", markers)))

dir.create("supp_figures/figureS6/FigS6B")
for (tier in rownames(all_i)) {
        png(paste0("supp_figures/figureS6/FigS6B/", tier, ".png"), height=700, width = 1000)
        par(lwd=6)
	barplot(t(rbind(mye_i[tier,], all_i[tier,]) / sum(all_i[tier,])), horiz = T, col = c("gray80", mature_cols))
        dev.off()

}

###################

dir.create("supp_figures/figureS6/FigS6C")
plot_virtual_facs(wells, "Csf1r", "Flt3", "supp_figures/figureS6/FigS6C/all.png", gates = list(flt3_gate$polygon, csf1r_gate$polygon, dp_gate$polygon))
plot_virtual_facs(wells, "SiglecH", "Flt3", "supp_figures/figureS6/FigS6C/flt3.png", filter = flt3_gate$gate, gates = list(pdc_gate$polygon))
plot_virtual_facs(wells, "Ly6C", "cKit", "supp_figures/figureS6/FigS6C/csf1r.png",filter = csf1r_gate$gate,  gates = list(cmop_gate$polygon))
plot_virtual_facs(wells, "CD11c", "cKit", "supp_figures/figureS6/FigS6C/dp.png", filter = dp_gate$gate, gates = list(predc_gate$polygon, cdp_gate$polygon, mdp_gate$polygon))

###################

load("gates/mye_clust_gates.Rda")
X = data.frame(x = sc_2d@x_cl, y = sc_2d@y_cl)
#hscs = gate_facs(X, "x", "y", n = 4)
#mpos = gate_facs(X, "x", "y", n = 4, colgate = cbind(hscs$gate,1))
#dntts = gate_facs(X, "x", "y", n = 4, colgate = cbind(hscs$gate,mpos$gate))
#neuts = gate_facs(X, "x", "y", n = 4, colgate = cbind(hscs$gate,mpos$gate, dntts$gate))
#monos = gate_facs(X, "x", "y", n = 4, colgate = cbind(hscs$gate,mpos$gate, dntts$gate, neuts$gate))
#pdcs = gate_facs(X, "x", "y", n = 4, colgate = cbind(hscs$gate,mpos$gate, dntts$gate, neuts$gate, monos$gate))
#cdcs = gate_facs(X, "x", "y", n = 4, colgate = cbind(hscs$gate,mpos$gate, dntts$gate, neuts$gate, monos$gate, pdcs$gate))
#dcps = gate_facs(X, "x", "y", n = 4, colgate = cbind(hscs$gate,mpos$gate, dntts$gate, neuts$gate, monos$gate, pdcs$gate, cdcs$gate))
#ct = hscs$gate + 2*mpos$gate + 3*dntts$gate + 4*neuts$gate + 5*monos$gate + 6*dcps$gate + 7*pdcs$gate + 8*cdcs$gate
#names(ct) = names(sc_2d@x_cl)
#save(ct, file = "gates/mye_clust_gates.Rda")
#plot(sc_2d@x_cl, sc_2d@y_cl, pch = 21, cex = 2, bg = rainbow(8)[ct])

igh_nm = grep("Igh", rownames(mye_umis), value = T)[1]
left_genes = c("Ifitm1", "Ifitm1", "Ifitm1", "Elane", "Cd34", "Cd34", "Dntt", "Mpo")
right_genes = c("Cd34", "Mpo", igh_nm, "Ltf", "Lyz1;Lyz2", "Cd7", "Siglech", "Cd74")
more_genes = choose_genes_from_clust(sc_cl, nms_per_clust = 12, nms_thresh = 1.7, bad_genes = bad_genes)
#more_genes = setdiff(, c("Camp", "S100a8", "S100a9", "Ngp", "Ltf", "Lcn2", cc_genes, ribo_genes))
foc_more = mye_foc[more_genes,]
hc = hclust(dist(foc_more), "ward.D2")
more_genes = more_genes[ hc$order]
total_ord = rep(0, length(ct)); names(total_ord) = names(ct)
zlim = c(min(foc_more), max(foc_more))
dir.create("supp_figures/figureS6/FigS6E")
clust_titles = rep(0, length(ct)); names(clust_titles) = names(ct)
good_clusts = names(ct)
for (gate in 1:max(ct)) {
  gate_clusts = names(which(ct == gate))
  gate_cells = names(sc_cl@clusts)[ sc_cl@clusts %in% gate_clusts]
  if (length(gate_clusts) > 1) {
	  genes = names(head(sort(apply(sc_cl@clust_fp[more_genes, gate_clusts],1,max), decreasing = T),25))
	  genes = genes[ apply(sc_cl@clust_fp[genes, gate_clusts],1,max) > 2.5]
  } else {
	  genes = names(head(sort(sc_cl@clust_fp[more_genes, gate_clusts],T),25))
	  genes = genes[ sc_cl@clust_fp[genes, gate_clusts] > 2.5]
  }
  genes = genes[ order(factor(genes, levels = more_genes))]
  clust_ord = names(sort(sc_cl@clust_fp[right_genes[gate],gate_clusts] - sc_cl@clust_fp[left_genes[gate],gate_clusts]))
  clust_titles[ clust_ord] = seq_along(clust_ord)
  total_ord[clust_ord] = 1:length(good_clusts)
  png(paste0("supp_figures/figureS6/FigS6E/", gate, ".png"), width = 1.5 * length(gate_cells) , height = 50 * length(genes))
  par(mar = c(0,0,0,0))
  X = foc_more[genes, gate_cells[ order(factor(sc_cl@clusts[gate_cells], levels = clust_ord))]]
  hc = hclust(dist(X), "ward.D2")
  hc$order = 1:length(genes)
  image(t(X[hc$order,]), zlim = zlim,
	col = colorRampPalette(c("white", "orange", "tomato","mediumorchid4", "midnightblue"))(130), axes = F)
  mtext(gsub(";.*", "", genes[ hc$order]), side = 2, las = 2, at = seq(0,1,length.out = length(genes)))
  mtext(gsub(";.*", "", genes[ hc$order]), side = 4, las = 2, at = seq(0,1,length.out = length(genes)))
  cls = cumsum(table(factor(sc_cl@clusts[gate_cells], levels = clust_ord))) / length(gate_cells)
  abline(v = cls, lwd = 8)
  #abline(h = gene_cls, lwd = 2, lty = 2)#, h = length(unique_genes) / length(genes))
  mtext(clust_ord, at = rowMeans(cbind(c(0,cls[-length(cls)]), cls)), side = 1)
  #mtext(good_clusts, at = rowMeans(cbind(c(0,cls[-max(model$MAP)]), cls)), side = 3)
  dev.off()
  write.table(rev(gsub(";.*", "", genes[ hc$order])), row.names = F, quote = F,col.names = F, file =
                paste0("supp_figures/figureS6/FigS6E/",gate, ".txt"))
  write.table(clust_ord, sep = "\t", row.names = F, quote = F,col.names = F, file =
                paste0("supp_figures/figureS6/FigS6E/",gate, "_ord.txt"))
}

graph = melt(sc_2d@clust_graph[ names(sc_2d@x_cl), names(sc_2d@x_cl)])
graph = graph[ graph$value == 1 & graph$X2 > graph$X1,]
png("supp_figures/figureS6/FigS6E.png", width = 1700, height = 1700)
plot(sc_2d@x, sc_2d@y, pch = 20, col = "darkgray", ylim = c(0,620), axes = F, xlab = "", ylab = "", cex = 3)
segments(sc_2d@x_cl[ as.character(graph$X1)], sc_2d@y_cl[ as.character(graph$X1)], sc_2d@x_cl[ as.character(graph$X2)], sc_2d@y_cl[ as.character(graph$X2)], lwd = 2.5)
points(sc_2d@x_cl, sc_2d@y_cl, pch = 21, bg = "white", col = rainbow(8, v = 0.8)[ct], cex = 9, lwd = 4)
text(sc_2d@x_cl, sc_2d@y_cl, clust_titles[names(sc_2d@x_cl)], col = rainbow(8, v = 0.8)[ct], cex = 3)
dev.off()

################

tier3_2d = sc_pipe_plots(scdb_init(basedir = "saved_work/tier3_clusts"))
tier3_cl = tier3_2d@scl
bad_cells= c()
P = .graph_con_comp(tier3_2d@clust_graph); large_comp = which.max(table(P))
good_clusts = names(which(P == large_comp))
good_cells = setdiff( names(which(tier3_cl@clusts > 0 & tier3_cl@clusts %in% good_clusts)), bad_cells)
outline = scr_find_outline(tier3_2d, reg = 0.75, cells = good_cells)

mye_gated_cells = intersect(names(sc_cl@clusts), rownames(wells))
dir.create("supp_figures/figureS6/Fig6D/")
xlim = quantile(outline[,1], c(0,1))
ylim = quantile(outline[,2], c(0,1))
for (gate in setdiff(unique(mye_gates$gate), c("other", NA))) {
  png(paste0("supp_figures/figureS6/Fig6D/", gate, ".png"), width = 1500, height = 1500)
  plot(tier3_2d@x[good_cells], tier3_2d@y[good_cells], type = "n", axes = F, xlab = "", ylab = "",
        xlim = xlim, ylim = ylim)
  points(sc_2d@x[mye_gated_cells], sc_2d@y[mye_gated_cells], pch = 21,
       bg =  ifelse(mye_gates[ mye_gated_cells, "gate"] == gate, "chocolate3", "white"),
       cex = ifelse(mye_gates[ mye_gated_cells, "gate"] == gate, 3, 0))
  points(outline[,1], outline[,2], type = "l", lwd = 4)
  dev.off()
}

###################

lfp = log2(sc_cl@clust_fp)
all_tfs = read.table("annotations/mouse.tf.list.symbols.txt", stringsAsFactors=F)[[1]]
all_tfs = intersect(rownames(lfp), all_tfs)
mon_genes = read.table("results/mon_genes.txt", stringsAsFactors=F, header = T, sep = "\t")[[1]]
neut_genes = read.table("results/neut_genes.txt", stringsAsFactors=F, header = T)[[1]]
mon_clusts = log2(tapply(colSums(mye_umis[mon_genes,]), sc_cl@clusts, sum) / tapply(colSums(mye_umis), sc_cl@clusts, sum) * 1000)
neut_clusts = log2(tapply(colSums(mye_umis[neut_genes,]), sc_cl@clusts, sum) / tapply(colSums(mye_umis), sc_cl@clusts, sum) * 1000)
png("supp_figures/figureS6/FigS6F.png", height = 1000, width = 1000)
plot(mon_clusts, neut_clusts, pch = 21, cex = 3, lwd = 2, bg = "gray60", axes = F, xlab = "", ylab = "")
axis(1); axis(2)
dev.off()


C = cor(t(lfp), cbind(mon_clusts, neut_clusts))
tfs = intersect(all_tfs, names(which(apply(C,1,max) > 0.5)))
tfs = tfs[ order(C[tfs,1] - C[tfs,2])]
tfs = tfs[c(seq_len(20), length(tfs) - rev(seq_len(20)) + 1)]
png("supp_figures/figureS6/FigS6G.png", height = 1500, width = 1000)
barplot(t(C[tfs,]), beside = T, col = c("mediumspringgreen", "darkolivegreen4"), border = NA, las = 2, horiz = T, space = c(0.1,0.1), names.arg = rep("", length(tfs)))
dev.off()
write.table(rev(tfs), row.names = F, quote = F, file = "supp_figures/figureS6/FigS6G.txt")