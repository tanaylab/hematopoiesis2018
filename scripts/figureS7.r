dir.create("supp_figures/figureS7")
message("generation Supplementary Fig 7")

pu1_col = "#9d1a7d"

tier3_2d = sc_pipe_plots(scdb_init(basedir = "saved_work/tier3_clusts"))
tier3_cl = tier3_2d@scl

#reposition_coords = read.delim("annotations/tier3_reposition_coords.txt", stringsAsFactor=F)
#for (i in seq_len(nrow(reposition_coords))) { tier3_2d = reposition_cc(tier3_2d, coords = as.numeric(reposition_coords[i,]))$sc_2d}

scdb = scdb_init(basedir="saved_work/ugi_clusts")
sc_2d = sc_pipe_plots(scdb)
sc_cl = sc_2d@scl

umis = as.matrix(sc_cl@scmat@mat)
ugi_list = read.delim("annotations/crispseq_count.txt", stringsAsFactor = F)

good_grna = c("CTL_LacZ", "CTRL_LacOp", "CTRL_Oris", "Cebpa_1", "Csf1r_1", "Csf1r_2", "Fcgr3_1", "Fcgr3_2", "Flt3_2_a", "Flt3_2_b",
          "Irf8-C", "Klf4_1", "Klf4_2", "Spi1_0")

cell_stats = sc_cl@scmat@cell_metadata
clusts = sc_cl@clusts
cells = names(clusts)
well2ugi = as.matrix(table(factor(ugi_list$well, levels = cells), factor(ugi_list$grna, levels = union(good_grna, sort(unique(ugi_list$grna))))))
well2ugi = well2ugi[,colSums(well2ugi) > 0]
grna_ord = order(colSums(well2ugi > 0), decreasing = T)
well2mouse = cell_stats$mouse.specimen.unique.ID; names(well2mouse) = rownames(cell_stats)
detected_ugi = well2ugi > 1
ugi_comb = rowSums(sweep(detected_ugi,2,2^(seq_len(ncol(detected_ugi)) - 1), "*"))
all_combs = unique(detected_ugi); rownames(all_combs) = rowSums(sweep(all_combs,2,2^(seq_len(ncol(all_combs)) - 1), "*"))
clones = as.vector(interaction(well2mouse[cells], ugi_comb[cells])); names(clones) = cells
clone_dist = table(clones, clusts)

ery_markers = c("Hba-a2", "Hbb-b1", "Beta-s"); neut_markers = c("S100a8", "S100a9", "Camp")
double_exp = c(colSums(sc_cl@clust_fp[ ery_markers, ] > 3) > 0) & c(colSums(sc_cl@clust_fp[ neut_markers, ] > 3) > 0)
good_clusts = setdiff(colnames(sc_cl@clust_fp), names(which(double_exp)))

tier3_annotations = read.table("annotations/tier3_annotations.txt", sep = "\t", stringsAsFactors = F, header = T)
markers = tier3_annotations[,1]; tier3_cols = c(tier3_annotations[,2], "darkgray"); thresh = tier3_annotations[,4]
above_exp = sc_cl@clust_fp[ markers, ] > thresh
clust_ass = max.col(t(above_exp), ties.method = "first"); clust_ass[colSums(above_exp) == 0] = length(markers) + 1
names(clust_ass) = colnames(sc_cl@clust_fp)

#######################

control_grna = c("CTRL_LacOp", "CTRL_Oris", "CTL_LacZ"); pos_grna = setdiff(colnames(well2ugi), c(control_grna, "Spi1_0"))
clone_names = sapply(strsplit(rownames(clone_dist), "\\."), "[", 2); names(clone_names) = rownames(clone_dist)
clone_comb = all_combs[clone_names,]; rownames(clone_comb) = rownames(clone_dist)
pos = c('Cebpa_1','Irf8-C')
neg = setdiff(colnames(clone_comb), c(pos, control_grna, "Spi1_0"))
#pos = c("Spi1_0")
clone_comb = all_combs[clone_names,]; rownames(clone_comb) = rownames(clone_dist)
clone2gene = rep("mix", nrow(clone_dist)); names(clone2gene) = rownames(clone_dist)
clone2gene[ rowSums(clone_comb[,pos]) == 1] = pos[ max.col(clone_comb[rowSums(clone_comb[,pos]) == 1,pos])]
clone2gene[ rowSums(clone_comb[,control_grna]) > 0 & rowSums(clone_comb[,setdiff(colnames(clone_comb), control_grna)]) == 0] = "CTL"
#clone2gene[ clone2gene == "CTL" & names(clone2gene) %in% grep("_7", names(clone2gene), value = T)] = "PU.1 CTL"
clone2gene[ rowSums(clone_comb[,setdiff(colnames(clone_comb), c(pos, control_grna))]) > 0] = "NEG"

good_clusts = setdiff(colnames(sc_cl@clust_fp), names(which(double_exp)))
cells = names(sc_cl@clusts)[cell_stats$virus.mix != "mix 4" & sc_cl@clusts %in% good_clusts]
umis_n = sweep(umis,2,colSums(umis),"/") * 1000
foc = log(1 + 7 * umis_n)

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

meta_order = c(16,2,11,3,14,15,13,12,4,6,10,5,7,8,1,9)
clust_ord = order(factor(clust_ass, levels = meta_order), secondary)

nms = choose_genes_from_clust(sc_cl, clust_ord, 10, 3, 60, bad_genes, c("Car1", "Mt2", "Klf1"), "max.col")

disp_cells = names(which(clones[cells] != "" & clone2gene[clones[cells]] != "NEG"))
message(length(disp_cells), " CRISP positive cells")
write.table(gsub(";.*", "", rev(nms)), row.names = F, col.names = F, quote = F, file = "supp_figures/figureS7/FigS7A_genes.txt")
cls = cumsum(table(factor(clust_ass[ sc_cl@clusts[disp_cells]], levels = meta_order))) / length(disp_cells)
png("supp_figures/figureS7/FigS7A.png", height=2000, width=2000)
par(mar = rep(0,4))
cell_ord = plot_sc_heatmap(sc_cl, nms, clust_ord, disp_cells, draw_cls = F)
abline(v = cls, lwd = 2, col = "gray40", lty = 2)
dev.off()

png("supp_figures/figureS7/FigS7A_bottom.png", height = 300, width = 1500)
par(mar = rep(0,4))
IM = t(detected_ugi[cell_ord,rev(c(control_grna, pos, neg))])
image(t(IM), col = c("white", "black"), axes = F)
abline(v = cls, lwd = 2, col = "gray40", lty = 2)
dev.off()

disp_cells = names(which(clones[cells] != "" & clone2gene[clones[cells]] == "NEG"))
message(length(disp_cells), " CRISP negative cells")
cls = cumsum(table(factor(clust_ass[ sc_cl@clusts[disp_cells]], levels = meta_order))) / length(disp_cells)

png("supp_figures/figureS7/FigS7B.png", height=2000, width=2000)
par(mar = rep(0,4))
cell_ord = plot_sc_heatmap(sc_cl, nms, clust_ord, disp_cells, draw_cls = F)
abline(v = cls, lwd = 2, col = "gray40", lty = 2)
dev.off()

png("supp_figures/figureS7/FigS7B_bottom.png", height = 300, width = 1500)
par(mar = rep(0,4))
IM = t(detected_ugi[cell_ord,rev(c(control_grna, pos, neg))])
image(t(IM), col = c("white", "black"), axes = F)
abline(v = cls, lwd = 2, col = "gray40", lty = 2)
dev.off()

##########################

marks = intersect(rownames(umis), rownames(tier3_cl@feat_mat))
d_old = t(log2(1+7*as.matrix(tier3_cl@scmat@mat[marks,])))
d_new = t(log2(1+7*umis[marks,]))
dcor = as.matrix(cor(t(d_new), t(d_old)))
m_knn = t(apply(dcor, 1, function(x) 1 - pmin(rank(-x)/50,1) ))
neigh = t(apply(m_knn, 1, function(x) names(tail(sort(x),n=50))))

old2new = data.frame(old = as.vector(neigh), new = rep(rownames(neigh), ncol(neigh)))
old2new$ass = clust_ass[ as.character(sc_cl@clusts[ old2new$new])]
X = tapply(old2new$ass, old2new$old, function(x) table(factor(x, levels = 1:15)))
Y = t(matrix(unlist(X), nrow = 15))
rownames(Y) = names(X)
Y = Y[rowSums(Y) > 50,]
old_ass = apply(Y,1,which.max)
pops = names(which(table(factor(old_ass, levels = 1:15)) > 50 & table(factor(clust_ass[as.character(sc_cl@clusts)], levels = 1:15)) > 50))
cells = c(names(old_ass)[ old_ass %in% pops], names(sc_cl@clusts)[ clust_ass[as.character(sc_cl@clusts)] %in% pops])
clusts = c(old_ass, clust_ass[ sc_cl@clusts]); names(clusts) = c(names(old_ass), names(sc_cl@clusts))
clusts = clusts[cells]
clusts = paste0(c(markers, "core")[ as.numeric(clusts)], ".", ifelse(names(clusts) %in% names(old_ass), "old", "new"))
names(clusts) = cells

sc_raw_mat = sc_pipe_build_mat(scdb_init(basedir="saved_work/all_cells"))
comb_mat = sc_raw_mat
comb_mat@mat = comb_mat@mat[,cells]
comb_cl = sc_cl; comb_cl@scmat = comb_mat
m = sc_to_bulk(comb_cl, clusts, bad_genes, cells, 50, choose_genes = F)
comb_cl@clusts = clusts
comb_cl@clust_fp = .calc_clusts_fp(comb_cl)
genes = setdiff(names(which(apply(comb_cl@clust_fp,1,max) > 2)), c(ribo_genes, cc_genes))
IM = log2(10 + m[genes,])
IM2 = t(scale(t(IM)))[,c(2,1,6,5,8,7,4,3,10,9,12,11)]
km = TGL.kmeans2(IM2, 20)
clusts = as.numeric(factor(km$cluster, levels = order(max.col(km$centers))))
png("supp_figures/figureS7/FigS7C.png", height = 1500, width = 1000)
par(mar = rep(0,4))
image(t(IM2[ order(clusts),]), col = colorRampPalette(c("blue", "white", "red"))(1000), axes = F)
cls = cumsum(table(clusts)) / length(clusts)
abline(h = cls)
dev.off()
zlim = quantile(IM2, c(0,1))
message("CRISPseq vs WT zlim: ", zlim[1], " - ", zlim[2])
####################

cell_stats = sc_cl@scmat@cell_metadata
cells = names(sc_cl@clusts)[cell_stats$virus.mix != "mix 4"]
clone_dist = table(clones[cells], sc_cl@clusts[cells])
clone_names = sapply(strsplit(rownames(clone_dist), "\\."), "[", 2); names(clone_names) = rownames(clone_dist)
big_clones = names(which(rowSums(clone_dist) > 20 & clone_names != "0"))

control_grna = c("CTRL_LacOp", "CTRL_Oris", "CTL_LacZ"); pos_grna = setdiff(colnames(well2ugi), c(control_grna, "Spi1_0"))
control_clones = intersect(big_clones, names(which(clone_names != "0" & rowSums(all_combs[clone_names, pos_grna]) == 0)))
#cells = names(which(clones != "" & clones %in% control_clones))

cond_dist = table(clones[cells],factor(clust_ass[sc_cl@clusts[cells]], levels = meta_order))
dist_n = cond_dist / rowSums(cond_dist)
em_skew = rowSums(dist_n[,1:4])
png("supp_figures/figureS7/FigS7D.png", height=1000, width=1000)
plot(rowSums(cond_dist), em_skew, log = "x", pch = 21, cex = 3, bg = ifelse(em_skew > 0.5, "red", "green"), axes = F, xlab = "", ylab = "")
axis(1); axis(2)
abline(h = 0.5, lty=2, lwd=2)
dev.off()

####################
pos = setdiff(names(which(colSums(detected_ugi) > 500)), control_grna)
clone_comb = all_combs[clone_names,]; rownames(clone_comb) = rownames(clone_dist)
clone2gene = rep("mix", nrow(clone_dist)); names(clone2gene) = rownames(clone_dist)
clone2gene[ rowSums(clone_comb[,pos]) == 1] = pos[ max.col(clone_comb[rowSums(clone_comb[,pos]) == 1,pos])]
clone2gene[ rowSums(clone_comb[,control_grna]) > 0 & rowSums(clone_comb[,setdiff(colnames(clone_comb), control_grna)]) == 0] = "CTL"
clone2gene[ rowSums(clone_comb[,setdiff(colnames(clone_comb), c(pos, control_grna))]) > 0] = "mix"
good_clones = names(which(clone2gene != "mix"))
cells = names(which(clones != "" & clones %in% good_clones & sc_cl@clusts %in% good_clusts))

lin_ord = c(11,3,14,15,13,7)
bad_groups = c(); good_groups = setdiff(seq_along(lin_ord), bad_groups)
cond_dist = table(clones[cells], factor(clust_ass[as.character(sc_cl@clusts[cells])], levels = lin_ord[ good_groups]))
n = rowSums(cond_dist)
clone_thresh = 2
dist_n = cond_dist / n
pvals = t(sapply(pos, function(y) {z = names(which(n >= clone_thresh & clone2gene[names(n)] %in% c("CTL",y)));
    sapply(as.character(lin_ord[good_groups]), function(x) wilcox.test(dist_n[z,x] ~ clone2gene[z])$p.value)}))
pvals_c = matrix((p.adjust(as.vector(pvals), "fdr")), nrow = nrow(pvals), ncol = ncol(pvals), dimnames = dimnames(pvals))
colnames(pvals_c) = c(markers, "none")[ as.numeric(colnames(pvals_c))]

IM = -log10(pvals_c)
zlim = c(0, max(IM))
message("all_grna_pval.png: ", zlim[1], " - ", zlim[2])
png("supp_figures/figureS7/FigS7E.png", height = 1000, width = 1200)
par(mar = c(5,0.5,0.5,0.5), fig = c(0,0.9,0,1))
image(t(IM), col = colorRampPalette(c("white", "tomato3"))(100), axes = F, zlim = zlim)
grid(nx = ncol(pvals_c), ny = nrow(pvals_c), lwd = 2.5)
par(fig = c(0.9,1,0,1), new = T, lwd = 6)
barplot(rowSums(table(clone2gene[clones[cells]], factor(clust_ass[as.character(sc_cl@clusts[cells])], levels = lin_ord[ good_groups]))[pos,]),
	names.arg = "", horiz = T, las = 2, yaxs = "i", col = "gray40", space = 0)
dev.off()

png("supp_figures/figureS7/FigS7E_cb.png", height = 100, width = 1000)
par(mar = rep(0,4))
image(matrix(1:100), col = colorRampPalette(c("white", "tomato3"))(100), axes = F)
dev.off()

pvals_melt = melt(pvals_c)
pvals_melt = pvals_melt[ pvals_melt$value < 0.05,]
write.table(pvals_melt, sep = "\t", quote = F, row.names = F)

###################

pos = c('Cebpa_1','Irf8-C')
clone_comb = all_combs[clone_names,]; rownames(clone_comb) = rownames(clone_dist)
clone2gene = rep("mix", nrow(clone_dist)); names(clone2gene) = rownames(clone_dist)
clone2gene[ rowSums(clone_comb[,pos]) == 1] = pos[ max.col(clone_comb[rowSums(clone_comb[,pos]) == 1,pos])]
clone2gene[ rowSums(clone_comb[,control_grna]) > 0 & rowSums(clone_comb[,setdiff(colnames(clone_comb), control_grna)]) == 0] = "CTL"
clone2gene[ rowSums(clone_comb[,setdiff(colnames(clone_comb), c(pos, control_grna))]) > 0] = "mix"
good_clones = names(which(clone2gene != "mix"))
cells = names(which(clones != "" & clones %in% good_clones & sc_cl@clusts %in% good_clusts))

cebpa_col = "#0b5ca8"; irf8_col = "#ab5617";
cond_dist = table(clones[cells], c(markers, "core")[clust_ass[as.character(sc_cl@clusts[cells])]])
dist_n = cond_dist / rowSums(cond_dist)

png("supp_figures/figureS7/FigS7F.png", height=1000, width=1000)
plot(dist_n[,"Ly86"], dist_n[,"Gstm1"], cex = 6, lwd=3, pch = 21, bg = c(cebpa_col, "gray80", irf8_col)[ as.numeric(factor(clone2gene[ rownames(dist_n)]))],
	axes = F, xlab = "", ylab = "")
axis(1); axis(2)
dev.off()