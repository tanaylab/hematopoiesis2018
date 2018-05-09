#############
# Figure 7
#############

message("generating figure 7")

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

dir.create("figures/figure7")

ery_markers = c("Hba-a2", "Hbb-b1", "Beta-s"); neut_markers = c("S100a8", "S100a9", "Camp")
double_exp = c(colSums(sc_cl@clust_fp[ ery_markers, ] > 3) > 0) & c(colSums(sc_cl@clust_fp[ neut_markers, ] > 3) > 0)
good_clusts = setdiff(colnames(sc_cl@clust_fp), names(which(double_exp)))

###############
# Project UGI on tier3 cluster model

ugi_coords = proj_ds_on_graph(tier3_2d, umis = umis,
            bg_cells = NULL, bw = 30, fn = "temp/proj_ugis.png", 
            reg = 10, bg_reg = 3, cex = 3, lwd = 2, markers = intersect(rownames(umis), rownames(sc_cl@feat_mat)))
save(ugi_coords, file = "saved_work/ugi_clusts/ugi_coords.Rda")

load("saved_work/ugi_clusts/ugi_coords.Rda")
x = ugi_coords$coords$x; y = ugi_coords$coords$y
x_cl = tapply(x[cells], sc_cl@clusts, mean)[good_clusts]
y_cl = tapply(y[cells], sc_cl@clusts, mean)[good_clusts]

tier3_annotations = read.table("annotations/tier3_annotations.txt", sep = "\t", stringsAsFactors = F, header = T)
markers = tier3_annotations[,1]; tier3_cols = c(tier3_annotations[,2], "darkgray"); thresh = tier3_annotations[,4]
above_exp = sc_cl@clust_fp[ markers, ] > thresh
clust_ass = max.col(t(above_exp), ties.method = "first"); clust_ass[colSums(above_exp) == 0] = length(markers) + 1
names(clust_ass) = colnames(sc_cl@clust_fp)

bad_cells= c()
P = .graph_con_comp(tier3_2d@clust_graph); large_comp = which.max(table(P))
tier3_clusts = names(which(P == large_comp))
tier3_cells = setdiff( names(which(tier3_cl@clusts > 0 & tier3_cl@clusts %in% tier3_clusts)), bad_cells)
outline = scr_find_outline(tier3_2d, reg = 0.75, cells = tier3_cells)

cells = names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% good_clusts))
above_exp = sc_cl@clust_fp[ markers, ] > thresh
clust_ass = max.col(t(above_exp), ties.method = "first"); clust_ass[colSums(above_exp) == 0] = length(markers) + 1
names(clust_ass) = colnames(sc_cl@clust_fp)
png("figures/figure7/Fig7B.png", height = 1500, width = 1200)
plot(tier3_2d@x, tier3_2d@y, axes = F, xlab = "", ylab = "", type = "n")
points(x[cells] ,y[cells], pch = 21, cex = 2, bg = tier3_cols[clust_ass[as.character(sc_cl@clusts[cells])]])
lines(outline[,1], outline[,2], lwd=4)
dev.off()

################
# Prove clonality

cell_stats = sc_cl@scmat@cell_metadata
cells = names(sc_cl@clusts)[cell_stats$virus.mix != "mix 4"]
clone_dist = table(clones[cells], sc_cl@clusts[cells])
clone_names = sapply(strsplit(rownames(clone_dist), "\\."), "[", 2); names(clone_names) = rownames(clone_dist)
big_clones = names(which(rowSums(clone_dist) > 20 & clone_names != "0"))

control_grna = c("CTRL_LacOp", "CTRL_Oris", "CTL_LacZ"); pos_grna = setdiff(colnames(well2ugi), c(control_grna, "Spi1_0"))
control_clones = intersect(big_clones, names(which(clone_names != "0" & rowSums(all_combs[clone_names, pos_grna]) == 0)))
cells = names(which(clones != "" & clones %in% control_clones))

lin_ord = c(2,11,16,3,14,15,13,12,4,6,10,5,7,8,1,9)
cond_dist = table(clones[cells],factor(clust_ass[sc_cl@clusts[cells]], levels = lin_ord))
#cond_dist = cond_dist[big_clones,]
dist_n = cond_dist / rowSums(cond_dist)
em_skew = log2((1 + rowSums(cond_dist[,1:4])) / (1 + rowSums(cond_dist[,-c(1:4)])))

png("figures/figure7/Fig7C.png", width = 800, height = 1000)
par(lwd = 3)
#clone_ord = order(rownames(dist_n) %in% control_clones, rowSums(dist_n[,1:4]))
clone_ord = order(rowSums(dist_n[,1:4]))
X = barplot(t(dist_n[ clone_ord,]) * 100, col = tier3_cols[lin_ord], ylim = c(-10,100), names.arg = rep("", nrow(cond_dist)))
dev.off()

cond_dist = table(clones[cells],factor(clust_ass[sc_cl@clusts[cells[sample(length(cells))]]], levels = lin_ord))
shuff_skew = log2((1 + rowSums(cond_dist[,1:4])) / (1 + rowSums(cond_dist[,-c(1:4)])))
dist_n = cond_dist / rowSums(cond_dist)
clone_ord = order(rowSums(dist_n[,1:4]))
png("figures/figure7/Fig7D.png", width = 800, height = 1000)
par(lwd=3)
barplot(t(dist_n[ clone_ord,]) * 100, col = tier3_cols[lin_ord], ylim = c(-10,100), names.arg = rep("", nrow(cond_dist)))
dev.off()

message("Clonality skew p-value: ", ks.test(em_skew, shuff_skew)$p.value)

##############
good_clusts = setdiff(colnames(sc_cl@clust_fp), names(which(double_exp)))
control_grna = c("CTRL_LacOp", "CTRL_Oris", "CTL_LacZ"); pos_grna = setdiff(colnames(well2ugi), c(control_grna))

cells = names(sc_cl@clusts)[cell_stats$virus.mix != "mix 4" & sc_cl@clusts %in% good_clusts]
clone_dist = table(clones[cells], clusts[cells])
clone_names = sapply(strsplit(rownames(clone_dist), "\\."), "[", 2); names(clone_names) = rownames(clone_dist)

pos = c('Cebpa_1','Irf8-C')
clone_comb = all_combs[clone_names,]; rownames(clone_comb) = rownames(clone_dist)
clone2gene = rep("mix", nrow(clone_dist)); names(clone2gene) = rownames(clone_dist)
if (length(pos) > 1) {
   clone2gene[ rowSums(clone_comb[,pos]) == 1] = pos[ max.col(clone_comb[rowSums(clone_comb[,pos]) == 1,pos])]
} else {
   clone2gene[ clone_comb[,pos]] = pos
}
clone2gene[ rowSums(clone_comb[,control_grna]) > 0 & rowSums(clone_comb[,setdiff(colnames(clone_comb), control_grna)]) == 0] = "CTL"
clone2gene[ rowSums(clone_comb[,setdiff(colnames(clone_comb), c(pos, control_grna))]) > 0] = "mix"
good_clones = names(which(clone2gene != "mix"))

#####################

bad_groups = c(); good_groups = setdiff(seq_along(lin_ord), bad_groups)
cond_dist = table(clones[cells], factor(clust_ass[as.character(sc_cl@clusts[cells])], levels = lin_ord[ good_groups]))
n = rowSums(cond_dist)
clone_thresh = 2
dist_n = cond_dist / n
pvals = t(sapply(pos, function(y) {z = names(which(n >= clone_thresh & clone2gene %in% c("CTL",y))); 
    sapply(as.character(lin_ord[good_groups]), function(x) wilcox.test(dist_n[z,x] ~ clone2gene[z])$p.value)}))
pvals_c = matrix((p.adjust(as.vector(pvals), "fdr")), nrow = nrow(pvals), ncol = ncol(pvals), dimnames = dimnames(pvals))
colnames(pvals_c) = c(markers, "none")[ as.numeric(colnames(pvals_c))]

y = pos
for (pop in as.character(c(11,14,15))) {
z = names(which(n > clone_thresh & clone2gene %in% c("CTL",y)))
heights = apply(dist_n[z,],2,tapply,factor(clone2gene[z], levels = c("CTL", y)), median)
dir.create("figures/figure7/Fig7E")
png(paste0("figures/figure7/Fig7E/", markers[as.numeric(pop)], ".png"), height = 1000, width = 1000)
plot(as.numeric(factor(clone2gene[z], levels = c("CTL", y))) + runif(length(z),-0.2,0.2), dist_n[z,pop] + 0.01, log = "y",
	xlab = "", ylab = "", xlim = c(0.5,3.5), axes = F, pch = 21, cex = 6, lwd = 3, ylim = c(0.01,1),
	bg = tier3_cols[as.numeric(pop)])
axis(2)
segments(c(0.8,1.8,2.8), heights[,pop] + 0.01, c(1.2,2.2,3.2), col = "black", lwd = 15)
dev.off()
}
