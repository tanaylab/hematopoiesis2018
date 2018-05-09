#############
# Figure 8
#############

message("generating figure 8")

dir.create("figures/figure8")
pu1_col = "#9d1a7d"

#############


sc_2d = sc_pipe_plots(scdb_init(basedir = "saved_work/pu1_clusts"))
#sc_2d = reposition_cc(sc_2d, coords = res$coords)$sc_2d
sc_cl = sc_2d@scl
pu1_cl = sc_cl;
pu1_umis = as.matrix(pu1_cl@scmat@mat)
pu1_n = sweep(pu1_umis,2,colSums(pu1_umis),"/") * 1000
pu1_foc = log(1 + 7 * pu1_n)
pu1_stats = pu1_cl@scmat@cell_metadata

bad_marks = c("Vpreb1", "Car2", "Hba-a2", "Beta-s", "Prss34", "Gstm1", "Cd74")
bad_clusts = names(which(apply( sc_cl@clust_fp[ bad_marks,],2,max) > 4))
good_clusts = setdiff(colnames(sc_cl@clust_fp), bad_clusts)
good_cells = names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% good_clusts))

outline = scr_find_outline(sc_2d, reg = 0.6, cells = good_cells)

png("figures/figure8/Fig8D.png", height = 1200, width = 1200)
plot(sc_2d@x[good_cells], sc_2d@y[good_cells], pch = 21, cex = 2, bg = ifelse(pu1_stats[good_cells, "treatment"] == "CRISPR_PU1", pu1_col, "gray60"),
	axes = F, xlab = "", ylab = "")
lines(outline[,1], outline[,2], lwd = 4)
dev.off()

nms = c("Itgam", "Mmp8", "Il1b", "Ccl6")

dir.create("figures/figure8/Fig8E")
filt_2d = sc_2d
filt_2d@scl@scmat@mat = filt_2d@scl@scmat@mat[,good_cells]
filt_2d@x = filt_2d@x[good_cells]; filt_2d@y = filt_2d@y[good_cells]
invisible(sapply(nms, function(x) plot_gene_2d(filt_2d, x, 1500, 1500, "figures/figure8/Fig8E", outline = outline,
        reg_factor = 10, positive_psize = 1, negative_psize = 1)))

###################
neut_genes = read.table("results/neut_genes.txt", stringsAsFactors=F, header = T)[[1]]
cells = names(sc_cl@clusts)
#pu1_comb = with(pu1_stats[cells,], paste0(treatment, ".", tier, ".", wells[cells, "Subject_ID"]))
pu1_comb = with(pu1_stats[cells,], paste0(treatment, ".", tier))
names(pu1_comb) = cells

lfp = log2(sc_cl@clust_fp)
wt_cells = rownames(pu1_stats)[ pu1_stats$treatment == "CRISPR_CTL"]
ga = "Ltf"; gb = "Ccl6"; a = lfp[ga,]; b = lfp[gb,]
g1 = intersect(wt_cells, names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% names(which(b > 2)))))
g2 = setdiff(intersect(wt_cells, names(which(sc_cl@clusts > 0 & sc_cl@clusts %in% good_clusts))), g1)
x = rowSums(pu1_n[,g2]) / length(g2) * min(length(g1), length(g2))
y = rowSums(pu1_n[,g1]) / length(g1) * min(length(g1), length(g2))
genes = setdiff(scr_chi_square_diff_genes(pu1_umis, g1 = g1,g2 = g2, pval = 1e-3, fdr = T), c(cc_genes, ribo_genes))
z = (y + 10) / (x + 10)
mature_genes = names(which(log2(z[genes]) > 1.5))
both = intersect(mature_genes, neut_genes)
mature_genes = setdiff(mature_genes, both)
neut_genes = setdiff(intersect(neut_genes, rownames(pu1_umis)), both)

sample_dist = table(pu1_comb, pu1_cl@clusts)
X = sample_dist/ rowSums(sample_dist)
ctl_pu1 = max.col(t(apply(X,2,tapply,c(1,1,2,2), sum))); names(ctl_pu1) = seq_along(ctl_pu1)
ly6g_3 = max.col(t(apply(X,2,tapply,c(1,2,1,2), sum))); names(ly6g_3) = seq_along(ly6g_3)
neut_clusts = tapply(colSums(pu1_umis[neut_genes,]), sc_cl@clusts, sum) / tapply(colSums(pu1_umis), sc_cl@clusts, sum) * 1000
mature_clusts = tapply(colSums(pu1_umis[mature_genes,]), sc_cl@clusts, sum) / tapply(colSums(pu1_umis), sc_cl@clusts, sum) * 1000
png("figures/figure8/Fig8F.png", height = 1000, width = 1000)
plot(mature_clusts, neut_clusts, pch = 21 + 2 * (ly6g_3 - 1), bg = c("gray60", pu1_col)[ctl_pu1], cex = 4, lwd = 3,
	axes = F, xlab = "", ylab = "", xlim = c(0,60), ylim = c(0,600))
axis(1); axis(2)
#text(mature_clusts, neut_clusts, good_clusts)
dev.off()

#################

iv_2d = sc_pipe_plots(scdb_init(basedir = "saved_work/iv_clusts"))
iv_cl = iv_2d@scl
iv_umis = as.matrix(iv_cl@scmat@mat)
iv_n = sweep(iv_umis,2,colSums(iv_umis),"/") * 1000
iv_foc = log(1 + 7 * iv_umis)
iv_stats = iv_cl@scmat@cell_metadata

bad_marks = c("Mcpt8", "H2-Aa", "Lgals3")
bad_clusts = names(which(apply( iv_cl@clust_fp[ bad_marks,],2,max) > 4))
iv_good_clusts = setdiff(colnames(iv_cl@clust_fp), bad_clusts)

cells = intersect(rownames(iv_stats)[ iv_stats$tier == "d9"], rownames(iv_stats)) #, names(which(iv_cl@clusts > 0 & iv_cl@clusts %in% good_clusts)))
iv_comb = with(iv_stats, paste0(treatment, ".", cytokines)); names(iv_comb) = rownames(iv_stats)

sample_dist = table(iv_comb[cells], iv_cl@clusts[cells])
X = sample_dist/ rowSums(sample_dist)
iv_ctl_pu1 = max.col(t(apply(X,2,tapply,c(1,1,2,2), sum))); names(iv_ctl_pu1) = seq_along(iv_ctl_pu1)
cytokines = max.col(t(apply(X,2,tapply,c(1,2,1,2), sum))); names(cytokines) = seq_along(cytokines)
iv_neut_clusts = tapply(colSums(iv_umis[intersect(rownames(iv_umis), neut_genes),cells]), iv_cl@clusts[cells], sum) /
        tapply(colSums(iv_umis[,cells]), iv_cl@clusts[cells], sum) * 1000
iv_mature_clusts = tapply(colSums(iv_umis[mature_genes,cells]), iv_cl@clusts[cells], sum) / tapply(colSums(iv_umis[,cells]), iv_cl@clusts[cells], sum) * 1000
png("figures/figure8/Fig8G.png", height = 1000, width = 1000)
plot(iv_mature_clusts, iv_neut_clusts, pch = 22 + 2 * (cytokines - 1), bg = c("gray60", pu1_col)[iv_ctl_pu1],
        cex = 4, lwd = 3, xlim = c(0,60), ylim = c(0,600), axes = F, xlab = "", ylab = "")
axis(1); axis(2)
#text(mature_clusts, neut_clusts, good_clusts)
dev.off()

###########

pu1_g1 = intersect(rownames(pu1_stats)[pu1_stats$tier == "Ly6g+" & pu1_stats$treatment == "CRISPR_PU1"], names(which(pu1_cl@clusts > 0 & pu1_cl@clusts %in% good_clusts)))
pu1_g2 = intersect(rownames(pu1_stats)[pu1_stats$tier == "Ly6g+" & pu1_stats$treatment == "CRISPR_CTL"], names(which(pu1_cl@clusts > 0 & pu1_cl@clusts %in% good_clusts)))

iv_g1 = intersect(rownames(iv_stats)[iv_stats$tier == "d9" & iv_stats$treatment == "CRISPR_PU1"], names(which(iv_cl@clusts > 0 & iv_cl@clusts %in% iv_good_clusts)))
iv_g2 = intersect(rownames(iv_stats)[iv_stats$tier == "d9" & iv_stats$treatment == "CRISPR_CTL"], names(which(iv_cl@clusts > 0 & iv_cl@clusts %in% iv_good_clusts)))

x1 = rowSums(pu1_n[,pu1_g2]) / length(pu1_g2) * min(length(pu1_g1), length(pu1_g2))
y1 = rowSums(pu1_n[,pu1_g1]) / length(pu1_g1) * min(length(pu1_g1), length(pu1_g2))

x2 = rowSums(iv_n[,iv_g2]) / length(iv_g2) * min(length(iv_g1), length(iv_g2))
y2 = rowSums(iv_n[,iv_g1]) / length(iv_g1) * min(length(iv_g1), length(iv_g2))

z1 = (y1 + 10) / (x1 + 10)
z2 = (y2 + 10) / (x2 + 10)

genes = intersect(rownames(pu1_n), rownames(iv_n))
disp_genes = names(which(apply(abs(log2(cbind(z1[genes], z2[genes]))),1,max) > 2))

png("figures/figure8/Fig8H.png", height = 1000, width = 1000)
plot(log2(z1[genes]), log2(z2[genes]), pch = 20, axes = F, xlab = "", ylab = "", col = ifelse(genes %in% disp_genes, "red", "navyblue"), cex = 1.5)
grid(lty = 2, col = "gray40", lwd = 2)
abline(h = 0, v = 0, col = "black", lwd = 4)
axis(1); axis(2)
dev.off()

message(length(disp_genes), " changed genes")
write.table(table( x = z1[disp_genes] >  1, y = z2[disp_genes] > 1) / length(disp_genes), sep = "\t", quote = F, col.names = NA)

png("figures/figure8/Fig8H_text.png", height = 1000, width = 1000)
plot(log2(z1[genes]), log2(z2[genes]), pch = 20, axes = F, xlab = "", ylab = "", col = "navyblue")
grid(); axis(1); axis(2)
text(log2(z1[disp_genes]), log2(z2[disp_genes]), disp_genes, col = "red")
dev.off()
