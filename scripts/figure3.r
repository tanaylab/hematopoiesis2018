#######################
# Figure 3
######################

message("generating figure 3")
dir.create("figures/figure3")

scdb = scdb_init(basedir="saved_work/core_clusts")
sc_2d = sc_pipe_plots(scdb)
sc_cl = sc_2d@scl
lfp = log2(sc_cl@clust_fp)
umis = as.matrix(sc_cl@scmat@mat)
umis_n = sweep(umis,2,colSums(umis),"/") * 1000
foc = log(1 + 7 * umis_n)
cell_stats = sc_cl@scmat@cell_metadata

####################

sc_clean_mat = sc_pipe_mat_to_clean_mat(scdb)
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
core_cells = colnames(umis)
core_i = table(sc_clean_mat@cell_metadata[core_cells, "tier"], factor(identity[core_cells], levels = c("none", markers)))

dir.create("figures/figure3/Fig3A")
for (tier in rownames(all_i)) {
	png(paste0("figures/figure3/Fig3A/", tier, ".png"), height=700, width = 1000)
	par(lwd=6)
	barplot(t(rbind(core_i[tier,], all_i[tier,]) / sum(all_i[tier,])), horiz = T, col = c("gray80", mature_cols))
	dev.off()
}


###################

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

outline = scr_find_outline(sc_2d, reg = 0.6, cells = good_cells)

library("reshape")
graph = melt(sc_2d@clust_graph[ names(sc_2d@x_cl), names(sc_2d@x_cl)])
graph = graph[ graph$value == 1 & graph$X2 > graph$X1,]
graph = graph[ graph$X1 %in% good_clusts & graph$X2 %in% good_clusts,]
png("figures/figure3/Fig3B.png", width = 1500, height = 1500)
plot(sc_2d@x[good_cells], sc_2d@y[good_cells], pch = 20, col = tier3_cols[ clust_ass[ sc_cl@clusts[good_cells]]], axes = F, xlab = "", ylab = "", cex = 3)
segments(sc_2d@x_cl[ as.character(graph$X1)], sc_2d@y_cl[ as.character(graph$X1)], sc_2d@x_cl[ as.character(graph$X2)], sc_2d@y_cl[ as.character(graph$X2)], 
	col = "gray20", lwd = 5)
points(sc_2d@x_cl[good_clusts], sc_2d@y_cl[good_clusts], pch = 21, bg = tier3_cols[clust_ass[good_clusts]], cex = 8, lwd = 2.5, col = "black")
lines(outline[,1], outline[,2], lwd = 6)
dev.off()

####################

sample_dist = table(cell_stats$tier == 7, sc_cl@clusts)
dist_n = sample_dist / rowSums(sample_dist)
rel_clusts = names(which(sample_dist[2,] >= 10))
rel_clusts = rel_clusts[ order(-dist_n[2, rel_clusts])]
rel_cells = rownames(cell_stats)[cell_stats$tier == 7 & sc_cl@clusts %in% rel_clusts]
s_clust = head(rel_clusts,1)

tier7_cells = rownames(cell_stats)[cell_stats$tier  == 7]
s_cells = rownames(cell_stats)[cell_stats$tier	== 7 & sc_cl@clusts == s_clust]

png("figures/figure3/Fig3C.png", height = 1000, width = 1000)
cells = intersect(tier7_cells, good_cells)
plot(sc_2d@x[good_cells], sc_2d@y[good_cells], type = "n", axes=F, xlab = "", ylab = "")
lines(outline[,1], outline[,2], lwd = 4)
points(sc_2d@x[setdiff(cells, s_cells)], sc_2d@y[setdiff(cells, s_cells)], pch = 21, cex = 2.5, bg = "gray80", lwd = 2)
points(sc_2d@x[intersect(cells, s_cells)], sc_2d@y[intersect(cells, s_cells)], pch = 21, cex = 2.5, bg = "gray40", lwd = 2)
dev.off()

####################

m = umis[,good_cells]
a = .downsamp(m,800)
f = a["Hlf",]>0
f_cands = rowSums(a[,f])>5  #just to save cpu
targs = t(apply(a[f_cands,],1,function(x) { fe = fisher.test(x>0, a["Hlf",]>0); return(c(fe$estimate, fe$p.value))}))

umis_n = t(t(umis) / colSums(umis)) * 1000
g1 = names(which(f)); g2 = names(which(!f))
x = rowSums(umis_n[,g2]) / length(g2) * min(length(g1), length(g2))
y = rowSums(umis_n[,g1]) / length(g1) * min(length(g1), length(g2))
z = (y + 10) / (x + 10)

disp_genes = intersect(names(which(log2(z) > 0)), names(which(p.adjust(targs[,2], "fdr") < 1e-3)))
s_genes = setdiff(disp_genes, bad_genes)
write.table(s_genes, row.names = F, quote = F, file = "results/s_genes.txt")

####################

png("figures/figure3/Fig3D.png", height = 1000, width = 700)
par(mar = c(4,6,1,1), lwd = 2)
barplot(rev(head(sort(log2(z[setdiff(disp_genes, "Hlf")]),decreasing = T),20)), 
	horiz = T, col = "black", las = 2, names.arg = rep("",20))
write.table(names(head(sort(log2(z[setdiff(disp_genes, "Hlf")]),decreasing = T),20)), row.names = F, col.names = F, quote = F, file = "figures/figure3/s_genes_bp.txt")
dev.off()

####################

s_score = colSums(log(1 + 7 * umis_n[ s_genes, ]))

png("figures/figure3/Fig3E.png", height = 1200, width = 1500)
plot(density(s_score[good_cells]), col = "gray80", main = "", axes = F, xlab = "", ylab = "", lwd = 13, xaxs = "i", ylim = c(0,0.055))
lines(density(s_score[s_cells]), col = "gray40", lwd=13)
abline(v = quantile(s_score, c(0.6,0.8,0.95,1)), lwd = 3, lty = 2)
axis(2); axis(1)
dev.off()

###################

map = ifelse(rownames(foc) %in% s_genes, "S", ifelse(rownames(foc) %in% cc_genes, "P", "none"))
val_mat = apply(foc, 2, tapply, map, sum)
b = 20
palette = c("white", "white", "#ddd1ba", "#7cc5e9","#ccd600", "#f8a921", "#d73531", "#cb0b78")
good_cells = good_cells[sample(length(good_cells))]
mat = val_mat
dir.create("figures/figure3/Fig3F-G")
for (val in c("S", "P")	) {
	vals = mat[val,]
	norm_val = as.numeric(cut(vals, unique(quantile(vals, (0:b)/b)), include.lowest = T));
	names(norm_val) = names(vals)
	cols = colorRampPalette(palette)(max(norm_val))
	png(paste0("figures/figure3/Fig3F-G/", val, ".png"), height = 1200, width = 1200)
	plot(sc_2d@x[good_cells], sc_2d@y[good_cells], cex = 1 + 0.7 * round(norm_val[good_cells] / max(norm_val) * 5), 
		pch = 21, bg = cols[norm_val[good_cells]], axes = F, xlab = "", ylab = "")
	lines(outline[,1], outline[,2], lwd = 4)
	dev.off()
}
png("figures/figure3/Fig3F-G/Fig3F-G_cb.png", height = 200, width = 1000)
par(mar = rep(0,4))
image(matrix(1:100), col = colorRampPalette(palette)(100), axes = F)
dev.off()
