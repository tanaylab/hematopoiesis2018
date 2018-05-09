message("generating figure 6")

dir.create("figures/figure6")
load("saved_work/myeloid_clusts/scrdb_data_projected_2d.RDa")
sc_cl = sc_2d@scl
graph = melt(sc_2d@clust_graph[ names(sc_2d@x_cl), names(sc_2d@x_cl)])
colnames(graph) = c("Var1", "Var2", "value")
graph = graph[ graph$value == 1 & graph$Var2 > graph$Var1,]

mye_umis = as.matrix(sc_cl@scmat@mat)
mye_n = sweep(mye_umis,2,colSums(mye_umis), "/") * 1000
mye_foc = log(1 + 7 * mye_n)
s_genes = setdiff(read.table("results/s_genes.txt", stringsAsFactors=F, header = T)[[1]], ribo_genes)
s_score = colSums(mye_foc[s_genes,])
s_clust = tapply(s_score, sc_cl@clusts, median)
gradient = colorRampPalette(c("gold", "chocolate3", "black"))(101)
s_val = 101 - round((s_clust - min(s_clust)) / (max(s_clust) - min(s_clust)) * 100)
ylim = quantile(sc_2d@y, c(0,1))
png("figures/figure6/Fig6A.png", width = 1700, height = 1700)
plot(sc_2d@x, sc_2d@y, pch = 20, col = "darkgray", ylim = ylim, axes = F, xlab = "", ylab = "", cex = 3)
segments(sc_2d@x_cl[ as.character(graph$Var1)], sc_2d@y_cl[ as.character(graph$Var1)], 
	sc_2d@x_cl[ as.character(graph$Var2)], sc_2d@y_cl[ as.character(graph$Var2)], lwd = 2.5)
points(sc_2d@x_cl, sc_2d@y_cl, pch = 21, bg = gradient[s_val[names(sc_2d@x_cl)]], cex = 8, lwd = 2.5, col = "black")
zlim = quantile(s_clust, c(0,1))
message("s score zlim: ", zlim[1], " - ", zlim[2])
dev.off()

png("figures/figure6/Fig6A_cb.png", height=100, width = 1000)
par(mar = rep(0,4))
image(matrix(1:100), axes = F, col = gradient)
dev.off()

#########

tier3_2d = sc_pipe_plots(scdb_init(basedir = "saved_work/tier3_clusts"))
tier3_cl = tier3_2d@scl

#reposition_coords = read.delim("annotations/tier3_reposition_coords.txt", stringsAsFactor=F)
#for (i in seq_len(nrow(reposition_coords))) { tier3_2d = reposition_cc(tier3_2d, coords = as.numeric(reposition_coords[i,]))$sc_2d}

bad_cells= c()
P = .graph_con_comp(tier3_2d@clust_graph); large_comp = which.max(table(P))
good_clusts = names(which(P == large_comp))
good_cells = setdiff( names(which(tier3_cl@clusts > 0 & tier3_cl@clusts %in% good_clusts)), bad_cells)
outline = scr_find_outline(tier3_2d, reg = 0.75, cells = good_cells)

png("figures/figure6/Fig6A_small.png", width = 1200, height = 1500)
plot(tier3_2d@x, tier3_2d@y, pch = 20, col = "darkgray", axes = F, xlab = "", ylab = "", cex = 3)
dev.off()

nms = c("Gpr56", "Dntt", "Flt3", "Ly6d", "Siglech", "Cd7", "Itgb7", "NAAA;Naaa", "Csf1r", "F13a1", "Elane",
	"Fcnb", "Irf8", "Cebpa", "Cebpe", "Sfpi1", "Satb1")
dir.create("figures/figure6/Fig6A-B")
trash = sapply(nms, function(x) 
	plot_gene_2d(sc_2d, x, 2000, 2000, "figures/figure6/Fig6A-B/", reg_factor = 10, positive_psize=1, negative_psize=1,
	outline = outline))

#############
message("Monocyte gradient")
lfp = log2(sc_cl@clust_fp)
ga = "Csf1r"; gb = "Ly86"; a = lfp[ga,]; b = lfp[gb,];
g1 = names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% names(which(a > 1 & b > 1))))
g2 = setdiff(names(sc_cl@clusts), g1) 
genes = setdiff(scr_chi_square_diff_genes(mye_umis, g1 = g1, g2 = g2, pval = 1e-9), c(cc_genes, ribo_genes))
x = rowSums(mye_n[,g2]) / length(g2) * min(length(g2), length(g1))
y = rowSums(mye_n[,g1]) / length(g1) * min(length(g2), length(g1))
z = (y + 10) / (x + 10)
mon_genes = setdiff(names(which(log2(z[genes]) > 1)), c("NAAA;Naaa", "Prss34"))
mon_cells = names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% names(which( a > 0))))
write.table(mon_genes, row.names = F, quote = F, file = "results/mon_genes.txt")

png("figures/figure6/Fig6D.png", height = 1000, width = 1000)
barplot(rev(head(sort(log2(z[mon_genes]),T),20)), horiz = T, col = "black", las = 2, names.arg = rep("",20))
dev.off()
write.table(names(head(sort(log2(z[mon_genes]),T),20)), row.names = F, col.names = F, quote = F, file = "figures/figure6/Fig6D.txt")

mon_score = colSums(mye_foc[mon_genes,])
mon_cells = names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% names(which( a > 0))))
#mon_cells = names(which(mon_score[mon_cells] > 0))
mon_clusts = names(table(sc_cl@clusts[mon_cells]))
mon_grad = c("cornsilk3", "#91aa83", "#49c97d", "mediumspringgreen")

clust_score = tapply(colSums(mye_umis[mon_genes,]), sc_cl@clusts, sum) / tapply(colSums(mye_umis), sc_cl@clusts, sum) * 1000
mon_breaks = c(0,10,18,23,1000)
clust_cut = cut(clust_score, mon_breaks); names(clust_cut) = names(clust_score)
score_cut = clust_cut[ as.character(sc_cl@clusts)]; names(score_cut) = names(sc_cl@clusts)

genes_shades = colorRampPalette(c("white", "orange", "tomato","mediumorchid4", "midnightblue"))(1000)
genes = setdiff(union(c("Serpina3f;Serpina3g", "Gpr56", "Ifitm1", "Eltd1", "Cd34", "Elane",
      "Ly86", "Csf1r", "Irf8"), names(head(sort(log2(z[mon_genes]),T),20))), 
      c("Ly6c1", "Ccdc109b", "Fn1", "Lgals1", "Glepp1;Ptpro", "L1cam"))

m = sc_to_bulk(sc_cl, sc_cl@clusts, bad_genes, mon_cells, min_comb= 0, choose_genes = F)
IM = log2(10 + m[genes,])
IM2 = t(scale(t(IM)))[, order(clust_score[colnames(IM)])]
genes = genes[ order(max.col(IM2))]

IM = mye_foc[genes, mon_cells[order(score_cut[mon_cells], mon_score[mon_cells])]]
cls = cumsum(table(score_cut[mon_cells])) / length(mon_cells)
png("figures/figure6/Fig6I.png", height = 2000, width = 1000)
par(mar = c(0,0,0,0))
image(t(IM), col = genes_shades, axes = F)
abline(v = cls, lwd=3)
dev.off()
write.table(rev(gsub(";.*", "", rownames(IM))), row.names = F, quote = F, col.names = F, file = "figures/figure6/Fig6I.txt")
zlim = quantile(IM, c(0,1))
message("monocyte zlim: ", zlim[1], " - ", zlim[2])

png("figures/figure6/Fig6I_cb.png", height = 100, width = 1000)
par(mar = rep(0,4))
image(matrix(as.numeric(score_cut[colnames(IM)])), axes = F, col = mon_grad)
abline(v = cls, lwd=3)
dev.off()

png("figures/figure6/Fig6F.png", height = 1000, width = 1000)
plot(sc_2d@x, sc_2d@y, ylim = c(0,650), pch = 20, col = "darkgray", cex = 2, axes = F, xlab = F, ylab = "")
points(sc_2d@x_cl[mon_clusts], sc_2d@y_cl[mon_clusts], pch = 21, bg = mon_grad[clust_cut[mon_clusts]], cex = 8, lwd = 3)
dev.off()

genes = c("Fcer1g",  "F13a1", "Klf4", "Irf5")
dir.create("figures/figure6/Fig6J")

#plot(tapply(colSums(mye_umis[,mon_cells]), sc_cl@clusts[mon_cells], sum) / as.vector(table(sc_cl@clusts[mon_cells])), fp)
for(gene in genes) {
        png(paste0("figures/figure6/Fig6J/", gene, ".png"), height = 700, width = 1000)
	fp = tapply(mye_umis[gene,], sc_cl@clusts, sum) /
                tapply(colSums(mye_umis), sc_cl@clusts, sum) * 1000
	plot(clust_score[mon_clusts], fp[mon_clusts], cex = 8, pch = 21, bg = mon_grad[ as.numeric(clust_cut[mon_clusts])], 
		log = "x", 
		lwd = 3, axes = F, xlab = "", ylab = "")
	axis(1); axis(2)
        dev.off()
}


#########

message("Neutrophil gradient")
ga = "Camp"; gb = "Elane"; a = lfp[ga,]; b = lfp[gb,]
#plot(a,b, log = "xy")
g1 = names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% names(which(a > 4))))
g2 = setdiff(names(sc_cl@clusts), g1) #names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% bottom_clusts))
genes = setdiff(scr_chi_square_diff_genes(mye_umis, g1 = g1, g2 = g2, pval = 1e-9), c(cc_genes, ribo_genes))
x = rowSums(mye_n[,g2]) / length(g2) * min(length(g2), length(g1))
y = rowSums(mye_n[,g1]) / length(g1) * min(length(g2), length(g1))
z = (y + 10) / (x + 10)
neut_genes = names(which(log2(z[genes]) > 2))
neut_cells = names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% names(which( a > 4 | b > 0))))
neut_grad = c("cornsilk3", colorRampPalette(c("green3", "darkgreen", "darkolivegreen4"))(3))
write.table(neut_genes, row.names = F, quote = F, file = "results/neut_genes.txt")

png("figures/figure6/Fig6C.png", height = 1000, width = 1000)
barplot(rev(head(sort(log2(z[neut_genes]),T),20)), horiz = T, col = "black", las = 2, names.arg = rep("",20))
dev.off()
write.table(names(head(sort(log2(z[neut_genes]),T),20)), col.names = F, row.names = F, quote = F, file = "figures/figure6/Fig6C.txt")

neut_score = colSums(mye_foc[neut_genes,])
#neut_cells = names(which(score[neut_cells] > 0))
neut_clusts = names(table(sc_cl@clusts[neut_cells]))
clust_score = tapply(colSums(mye_umis[neut_genes,]), sc_cl@clusts, sum) / tapply(colSums(mye_umis), sc_cl@clusts, sum) * 1000

neut_breaks = c(0,5,10,50,1000)
clust_cut = cut(clust_score, neut_breaks); names(clust_cut) = names(clust_score)
score_cut = clust_cut[ as.character(sc_cl@clusts)]; names(score_cut) = names(sc_cl@clusts)

genes_shades = colorRampPalette(c("white", "orange", "tomato","mediumorchid4", "midnightblue"))(1000)
genes = setdiff(union(c("Serpina3f;Serpina3g", "Gpr56", "Ifitm1", "Eltd1", "Cd34", "Elane",
      "Fcnb", "Gstm1", "Mpo"), names(head(sort(log2(z[neut_genes]),T),20))),
      c("6430548M08Rik", "AK028782", "1100001G20Rik", "Lyz1;Lyz2", "Mapk13", "Mmp9"))

m = sc_to_bulk(sc_cl, sc_cl@clusts, bad_genes, neut_cells, min_comb= 0, choose_genes = F)
IM = log2(10 + m[genes,])
IM2 = t(scale(t(IM)))[, order(clust_score[colnames(IM)])]
genes = genes[ order(max.col(IM2))]

IM = mye_foc[genes, neut_cells[order(score_cut[neut_cells], neut_score[neut_cells])]]
cls = cumsum(table(score_cut[neut_cells])) / length(neut_cells)
png("figures/figure6/Fig6G.png", height = 2000, width = 1000)
par(mar = c(0,0,0,0))
image(t(IM), col = genes_shades, axes = F)
abline(v = cls, lwd=3)
dev.off()
write.table(rev(gsub(";.*", "", rownames(IM))), row.names = F, quote = F, col.names = F, file = "figures/figure6/Fig6G.txt")
zlim = quantile(IM, c(0,1))
message("neutrophil zlim: ", zlim[1], " - ", zlim[2])

png("figures/figure6/Fig6E.png", height = 1000, width = 1000)
plot(sc_2d@x, sc_2d@y, ylim = c(0,650), pch = 20, col = "darkgray", cex = 2, axes = F, xlab = F, ylab = "")
points(sc_2d@x_cl[neut_clusts], sc_2d@y_cl[neut_clusts], pch = 21, bg = neut_grad[clust_cut[neut_clusts]], cex = 8, lwd = 3)
dev.off()

genes = c("Elane", "Ngp", "Cebpe", "Gfi1")
dir.create("figures/figure6/Fig6H")

#plot(tapply(colSums(mye_umis[,neut_cells]), sc_cl@clusts[neut_cells], sum) / as.vector(table(sc_cl@clusts[neut_cells])), fp)
for(gene in genes) {
        png(paste0("figures/figure6/Fig6H/", gene, ".png"), height = 700, width = 1000)
        fp = tapply(mye_umis[gene,], sc_cl@clusts, sum) /
                tapply(colSums(mye_umis), sc_cl@clusts, sum) * 1000
        plot(clust_score[neut_clusts], fp[neut_clusts], cex = 8, pch = 21, bg = neut_grad[ as.numeric(clust_cut[neut_clusts])],
                log = "x",
                lwd = 3, axes = F, xlab = "", ylab = "")
        axis(1); axis(2)
        dev.off()
}

