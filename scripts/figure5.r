##################
# Figure 5 and Supplementary Fig 5
##################

message("generating figure 5")
sc_2d = sc_pipe_plots(scdb_init(basedir = "saved_work/tier3_clusts"))
sc_cl = sc_2d@scl


dir.create("figures/figure5")
dir.create("supp_figures/figureS5")

#reposition_coords = read.delim("annotations/tier3_reposition_coords.txt", stringsAsFactor=F)
#for (i in seq_len(nrow(reposition_coords))) { sc_2d = reposition_cc(sc_2d, coords = as.numeric(reposition_coords[i,]))$sc_2d}

tier3_annotations = read.table("annotations/tier3_annotations.txt", sep = "\t", stringsAsFactors = F, header = T)
markers = tier3_annotations[,1]; tier3_cols = c(tier3_annotations[,2], "darkgray"); thresh = tier3_annotations[,3]
above_exp = sc_cl@clust_fp[ markers, ] > thresh
clust_ass = max.col(t(above_exp), ties.method = "first"); clust_ass[colSums(above_exp) == 0] = length(markers) + 1
names(clust_ass) = colnames(sc_cl@clust_fp)

bad_cells= c()

P = .graph_con_comp(sc_2d@clust_graph); large_comp = which.max(table(P))
good_clusts = names(which(P == large_comp))
good_cells = setdiff( names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% good_clusts)), bad_cells)
outline = scr_find_outline(sc_2d, reg = 0.75, cells = good_cells)

con_2d = sc_2d
con_2d@scl@scmat@mat = con_2d@scl@scmat@mat[,good_cells]
con_2d@scl@clusts = con_2d@scl@clusts[good_cells]
con_2d@x = con_2d@x[good_cells]; con_2d@y = con_2d@y[good_cells]

invisible(sapply(c("Epor", "Csf3r"), function(x) plot_gene_2d(con_2d, x, 2000, 2000,
      "figures/figure5", reg_factor = 15, pt_shades = colorRampPalette(c("black")),
      outline = outline, positive_psize = 0.5, negative_psize = 0.5)))

######################

sc_clean_mat = sc_pipe_mat_to_clean_mat(scdb_init(basedir="saved_work/tier3_cytokines"))
pert_umis = as.matrix(sc_clean_mat@mat)
pert_cells = colnames(pert_umis)
cell_stats = sc_clean_mat@cell_metadata

#outline = scr_find_outline(prog_2d, reg = 0.5)
fg_cells = rownames(cell_stats)[ cell_stats$treatment == "Epo"]
bg_cells = rownames(cell_stats)[ cell_stats$treatment == "Epo (control)"]

shades = colorRampPalette(c("gray40", "gray40", "orange", "red"))(101)
epo_map = proj_ds_on_graph(con_2d, 50, umis = pert_umis[, union(fg_cells, bg_cells)], outline = outline,
            bg_cells = bg_cells, bw = 30, fn = "figures/figure5/Fig5C.png", reg = 10,
            bg_reg = 5, cex = 3.5, lwd = 2, clust_shades = shades)

fg_cells = rownames(cell_stats)[ cell_stats$treatment == "G-CSF"]
bg_cells = rownames(cell_stats)[ cell_stats$treatment == "G-CSF (control)"]

gcsf_map = proj_ds_on_graph(con_2d, 50, umis = pert_umis[, union(fg_cells, bg_cells)], outline = outline,
            bg_cells = bg_cells, bw = 30, fn = "figures/figure5/Fig5E.png", reg = 10,
            bg_reg = 5, cex = 3.5, lwd = 2, clust_shades = shades)


png("figures/figure5/colorbar.png", height = 100, width = 1000)
par(mar = rep(0,4))
image(matrix(1:100), axes = F, col = shades)
dev.off()

cells = rownames(cell_stats)
map = proj_ds_on_graph(sc_2d, 50, umis = pert_umis[, cells], outline = NULL,
            bg_cells = NULL, bw = 30, fn = "temp/all.png", reg = 10,
            bg_reg = 5, cex = 3.5, lwd = 2)

tier3_dist = table(as.vector(cell_stats$treatment), factor(c(markers, "core")[ clust_ass[map$clust ]], levels = c(markers, "core")))

fates = c(11,3,14,15,13,7)
cols = tier3_cols[fates]
dist_n = tier3_dist / rowSums(tier3_dist)
fg_samples = c(1,3); bg_samples = c(2,4)
enr = dist_n[ fg_samples,] / dist_n[ bg_samples,]
se = sqrt(1/tier3_dist[fg_samples,] - 1/rowSums(tier3_dist[fg_samples,]) +
            1/tier3_dist[bg_samples,] -1/ rowSums(tier3_dist[bg_samples,]))
enr = enr[,fates]; se = se[,fates]
alpha = 0.05
z = qnorm(1-alpha/2)
ci_l = enr*exp(-z*se)
ci_h = enr*exp(+z*se)

library("gplots")
dir.create("figures/figure5/Fig5DF")
for (samp in rownames(enr)) {
  png(paste0("figures/figure5/Fig5DF/", samp, ".png"), height = 1000, width = 1000)
  par(lwd=6)
  barplot2(log2(enr[samp,]), horiz = F, col = cols,names.arg = "",
           las = 2, plot.ci = T, ci.l = log2(ci_l[samp,]), ci.u = log2(ci_h[samp,]),
           ci.lwd = 3, lwd = 3, ylim = c(-2.5,3))
  dev.off()
}

message("Epo p-vals:")
cells = rownames(cell_stats)[ cell_stats$treatment %in% c("Epo", "Epo (control)")]
X = p.adjust(sapply(markers[fates], function(x)
fisher.test(c(markers, "core")[ clust_ass[map$clust[cells]]] == x, as.vector(cell_stats[cells, "treatment"]))$p.value), "fdr")
write.table(cbind(X, X < 0.05, X < 1e-3, X < 1e-5), sep = "\t", quote = F)

message("G-CSF p-vals:")
cells = rownames(cell_stats)[ cell_stats$treatment %in% c("G-CSF", "G-CSF (control)")]
X = p.adjust(sapply(markers[fates], function(x)
fisher.test(c(markers, "core")[ clust_ass[map$clust[cells]]] == x, as.vector(cell_stats[cells, "treatment"]))$p.value), "fdr")
write.table(cbind(X, X < 0.05, X < 1e-3, X < 1e-5), sep = "\t", quote = F)

############

fates = c("Car1", "Hba-a2", "Ly86", "Gran", "Vpreb1")
names(tier3_cols) = c(markers, "none")
tier3_cols["Gran"] = "green4"

comb = as.vector(cell_stats$treatment); names(comb) = rownames(cell_stats)
comb[ comb %in% c("Epo (control)", "G-CSF (control)")] = "CTL"

cells = names(comb) #rownames(cell_stats)[ cell_stats$treatment %in% c("Epo", "Epo (control)")]
identity = c(markers, "core"); names(identity) = identity
identity[c("Gstm1", "Fcnb")] = "Gran"
#mouse_dist = table(as.vector(cell_stats[cells, "mouse.specimen.unique.ID"]), factor(c(markers, "core")[ clust_ass[map$clust[cells] ]], levels = fates))
mouse_dist = table(as.vector(cell_stats[cells, "mouse.specimen.unique.ID"]), factor(identity[ clust_ass[map$clust[cells] ]], levels = fates))
dist_n = mouse_dist / rowSums(mouse_dist)
dist_n = dist_n[,fates]
mouse2treatment = with(cell_stats[cells,], table(as.vector(mouse.specimen.unique.ID), comb))
mouse_ass = colnames(mouse2treatment)[ max.col(mouse2treatment)]; names(mouse_ass) = rownames(mouse2treatment)

dir.create("supp_figures/figureS5/FigS5D")

for (fate in fates) {
	png(paste0("supp_figures/figureS5/FigS5D/", fate, ".png"), height = 1000, width = 700)
	plot(runif(length(mouse_ass), -0.2, 0.2) + as.numeric(factor(mouse_ass)), dist_n[,fate], pch = 21, bg = tier3_cols[fate],
		xlim = c(0.5,3.5), cex = 8, axes = F, xlab = "", ylab = "")
	segments((1:3) - 0.3, tapply(dist_n[,fate], factor(mouse_ass), mean), (1:3) + 0.3, col = "black", lwd = 10)
	axis(2)
	dev.off()
}

############

fates = c("Car1", "Hba-a2", "Ly86", "Gstm1", "Fcnb", "Vpreb1")
cols = tier3_cols; names(cols) = markers
cells = names(map$clust)[ c(markers, "core")[clust_ass[map$clust ]] %in% fates]
clusts = paste0(c(markers, "core")[  clust_ass[map$clust[cells] ]], ".", 
	ifelse(cell_stats[cells, "treatment"] %in% c("G-CSF (control)", "Epo (control)"), "Control", as.vector(cell_stats[cells, "treatment"])))
names(clusts) = cells

X = expand.grid(fates, c("Control", "Epo", "G-CSF"))
X$comb = with(X, paste0(Var1, ".", Var2))

comb_mat = sc_clean_mat
comb_mat@mat = comb_mat@mat[,cells]
comb_cl = sc_cl; comb_cl@scmat = comb_mat
comb_cl@clusts = clusts
m = sc_to_bulk(comb_cl, comb_cl@clusts, bad_genes, min_comb = 0, choose_genes = F) 

X = merge(X, t(m[c("Ifitm1","Hlf"),]),by.x = "comb", by.y = 0, all.x = T)
X = X[order(factor(X$Var1, levels = fates)),]
gene = "Ifitm1"
z1 = (X[X$Var2 == "Epo", gene] + 10) / (X[X$Var2 == "Control", gene] + 10)
z2 = (X[X$Var2 == "G-CSF", gene] + 10) / (X[X$Var2 == "Control", gene] + 10)
ylim = quantile(log2(c(z1,z2)),c(0,1))
png("figures/figure5/Fig5G.png", height = 1000, width = 1000)
par(lwd = 6)
barplot(log2(z1), col = cols[fates], ylim = ylim)
dev.off()
png("figures/figure5/Fig5H.png", height = 1000, width = 1000)
par(lwd = 6)
barplot(log2(z2), col = cols[fates], ylim = ylim)
dev.off()

fates = c("Car1", "Hba-a2", "Ly86", "Gstm1", "Fcnb", "Vpreb1")
X = cbind(pert_umis["Ifitm1",], colSums(pert_umis) - pert_umis["Ifitm1",])
pvals = rep(0, length(fates)); names(pvals) = fates
for (fate in fates) {
	cells = rownames(cell_stats)[ cell_stats$treatment %in% c("Epo", "Epo (control)") & c(markers, "none")[ clust_ass[map$clust]] == fate]
	pvals[fate] = fisher.test(apply(X[cells,], 2, tapply, as.vector(cell_stats[cells, "treatment"]), sum))$p.value
}
message("Ifitm1 Epo pvals:")
padj = p.adjust(pvals, "fdr")
write.table(cbind(padj, padj < 0.05, padj < 1e-3, padj < 1e-5), sep = "\t", quote = F)

pvals = rep(0, length(fates)); names(pvals) = fates
for (fate in fates) {
	cells = rownames(cell_stats)[ cell_stats$treatment %in% c("G-CSF", "G-CSF (control)") & c(markers, "none")[ clust_ass[map$clust]] == fate]
	pvals[fate] = fisher.test(apply(X[cells,], 2, tapply, as.vector(cell_stats[cells, "treatment"]), sum))$p.value
}
message("Ifitm1	G-CSF pvals:")
padj = p.adjust(pvals, "fdr")
write.table(cbind(padj, padj < 0.05, padj < 1e-3, padj < 1e-5), sep = "\t", quote = F)

#########################

sc_clean_mat = sc_pipe_mat_to_clean_mat(scdb_init(basedir = "saved_work/tier7_cytokines"))
tier7_umis = as.matrix(sc_clean_mat@mat)
tier7_n = sweep(tier7_umis, 2, colSums(tier7_umis), "/") * 1000
tier7_foc = log(1 + 7 * tier7_n)
cell_stats = sc_clean_mat@cell_metadata

####################
scdb = scdb_init(basedir="saved_work/core_clusts")
sc_2d = sc_pipe_plots(scdb)
sc_cl = sc_2d@scl

map = proj_ds_on_graph(sc_2d, umis = tier7_umis, outline = NULL,
            bg_cells = NULL, bw = 30, fn = "temp/tier7_all.png", 
            reg = 10, bg_reg = 5, cex = 2, lwd = 2, markers = intersect(rownames(sc_cl@feat_mat), rownames(tier7_umis)))


s_genes = read.table("results/s_genes.txt", stringsAsFactors=F, header = T)[[1]]
s_score = colSums(log(1 + 7 * tier7_n[s_genes, ]))

png("supp_figures/figureS5/FigS5E.png", height = 1000, width = 700)
par(lwd=6)
boxplot(s_score ~ cell_stats$treatment, axes = F, col = "gray40")
axis(2)
dev.off()
message("tier7 Epo vs. control p = ", ks.test(s_score[cell_stats$treatment == "Epo"], s_score[cell_stats$treatment == "-"])$p.value)
message("tier7 G-CSF vs. control p = ", ks.test(s_score[cell_stats$treatment == "G-CSF"], s_score[cell_stats$treatment == "-"])$p.value)

cells = names(which(map$clust == s_clust))
good_genes = setdiff(rownames(tier7_n), c("Hba-a2", "Beta-s", "Hbb-b1"))
m = t(apply(tier7_n[good_genes,cells], 1, tapply, cell_stats[cells, "treatment"], sum))
sizes = as.vector(table(cell_stats[cells, "treatment"]))
m = sweep(m,2,sizes,"/") * min(sizes)

x = log2((m[,3] + 10) / (m[,1] + 10)); y = log2((m[,2] + 10) / (m[,1] + 10))
disp_genes = names(which(pmax(abs(x), abs(y)) > 1))
message(length(disp_genes), " changed genes")
write.table(table(x = x[disp_genes] > 0, y = y[disp_genes] > 0) / length(disp_genes), sep = "\t", quote = F, col.names=NA)
png("figures/figure5/Fig5K.png", height = 1000, width = 1000)
plot(x,y, pch = 20, col = ifelse(names(x) %in% disp_genes, "red", "navyblue"), cex = 1 + (names(x) %in% disp_genes), axes = F, xlab = "", ylab = "")
grid(lty = 2, col = "gray40", lwd = 2)
abline(h = 0, v = 0, col = "black", lwd = 4)
axis(1); axis(2)
dev.off()

png("figures/figure5/Fig5K_text.png", height = 2000, width = 2000)
plot(x,y, pch = 20, type = "n")
abline(h = 0, v = 0, col = "red", lty = 2)
text(x[disp_genes], y[disp_genes], gsub(";.*", "", disp_genes), col = "red")
dev.off()


##############

lfp = log2(sc_cl@clust_fp)
bad_marks = c("Ccl5", "Vpreb3", "Vpreb1", "Camp")
outliers = names(which(apply(lfp[bad_marks,], 2, max) > 2))
good_clusts = setdiff(colnames(lfp), outliers)
good_cells = names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% good_clusts))
outline = scr_find_outline(sc_2d, reg = 0.6, cells = good_cells)


fg_cells = rownames(cell_stats)[ cell_stats$treatment == "Epo"]
bg_cells = rownames(cell_stats)[ cell_stats$treatment == "-"]

epo_map = proj_ds_on_graph(sc_2d, umis = tier7_umis[, union(fg_cells, bg_cells)], outline = outline,
            bg_cells = bg_cells, bw = 30, fn = "figures/figure5/Fig5I.png", 
            reg = 10, bg_reg = 5, cex = 5.5, lwd = 2, 
            clust_shades = colorRampPalette(c("gray40", "gray40", "orange", "red"))(101), 
	    markers = intersect(rownames(sc_cl@feat_mat), rownames(tier7_umis)))

fg_cells = rownames(cell_stats)[ cell_stats$treatment == "G-CSF"]
bg_cells = rownames(cell_stats)[ cell_stats$treatment == "-"]

gcsf_map = proj_ds_on_graph(sc_2d, umis = tier7_umis[, union(fg_cells, bg_cells)], outline = outline,
            bg_cells = bg_cells, bw = 30, fn = "figures/figure5/Fig5J.png",
            reg = 10, bg_reg = 5, cex = 5.5, lwd = 2, 
            clust_shades = colorRampPalette(c("gray40", "gray40", "orange", "red"))(101),
	    markers = intersect(rownames(sc_cl@feat_mat), rownames(tier7_umis)))

