##########################
# Figure 1
##########################

message("generating figure 1")

###################
# Figure 1B-D - compare alternative sorts

dir.create("figures/figure1")
#load("tier3_clusts/tier3_clusts.Rda")
batch_stats = read.delim("annotations/batch_stats.txt", stringsAsFactors = F)
tier3_B = batch_stats[ batch_stats$tier == 3 & batch_stats$treatment == "-" & batch_stats$noise.estimation <= 0.05, "amplification.batch"]
tier1_B = batch_stats[ batch_stats$tier == 1 & batch_stats$noise.estimation <= 0.05, "amplification.batch"]
tier2_B = batch_stats[ batch_stats$tier == 2 & batch_stats$noise.estimation <= 0.05, "amplification.batch"]


tier1_cl = sc_pipe_mat_to_cmods(scdb_init(basedir = "saved_work/tier1_clusts"))
tier2_cl = sc_pipe_mat_to_cmods(scdb_init(basedir = "saved_work/tier2_clusts"))
tier3_cl = sc_pipe_mat_to_cmods(scdb_init(basedir = "saved_work/tier3_clusts"))
tier1_umis = as.matrix(tier1_cl@scmat@mat)
tier2_umis = as.matrix(tier2_cl@scmat@mat)
tier3_umis = as.matrix(tier3_cl@scmat@mat)

lin_ord = c(2,9,3,12,16,13,11,10,4,7,6,5,8,1,14,15)
markers = c("Ccl5", "Pf4", "Hba-a2","Prss34", "Siglech",  "Cd74", "Prg2", "Vpreb1", "Car1", "Ltf", "Fcnb", "Ly86", "Gstm1", "Myl4", "Fcrla", "Lyz1;Lyz2")[lin_ord]
mature_cols = c("darkorchid4", "darksalmon", "indianred4", "goldenrod1",
           "darkcyan", "cyan3", "goldenrod3", "navyblue", "indianred3", "darkolivegreen4", "darkgreen", "limegreen", "green4", "dodgerblue3", "darkslateblue", "palegreen4")[lin_ord]

tier1_cells = intersect(colnames(tier1_umis), rownames(wells)[ wells$Amp_batch_ID %in% tier1_B])
tier1_gated_cells = intersect(tier1_cells, rownames(wells)[ !lin_gate$gate])

marker.exp = as.data.frame(t(tier1_umis[ markers, tier1_gated_cells]))
wmax = max.col(marker.exp, ties.method = "first")
expressed = rowSums(marker.exp > 1)
identity = ifelse(expressed, markers[wmax], "none")
tier1_i = table(factor(identity, levels = c("none", markers)))

tier2_cells = intersect(colnames(tier2_umis), rownames(wells)[ wells$Amp_batch_ID %in% tier2_B])
tier2_gated_cells = intersect(tier2_cells, rownames(wells)[ !ckit_gate$gate])

marker.exp = as.data.frame(t(tier2_umis[ markers, tier2_gated_cells]))
wmax = max.col(marker.exp, ties.method = "first")
expressed = rowSums(marker.exp > 1)
identity = ifelse(expressed, markers[wmax], "none")
tier2_i = table(factor(identity, levels = c("none", markers)))

marker.exp = as.data.frame(t(tier3_umis[ markers, ]))
wmax = max.col(marker.exp, ties.method = "first")
expressed = rowSums(marker.exp > 1)
identity = ifelse(expressed, markers[wmax], "none")
tier3_i = table(factor(identity, levels = c("none", markers)))

all_i = cbind(tier1_i, tier2_i, tier3_i)
all_i  = t(t(all_i) / colSums(all_i))
all_i = t(t(all_i) * c(length(tier1_gated_cells) / length(tier1_cells), length(tier2_gated_cells) / length(tier2_cells), 1))
all_i = rbind( all_i, 1 - colSums(all_i))

png("figures/figure1/Fig1B.png", width = 1500, height = 200)
par(lwd = 2.5, mar = rep(1,4))
barplot(t(t(all_i[,1]) / sum(all_i[,1])), col = c("gray80", mature_cols, "white"), las = 2, horiz = T, axes = F)
axis(side = 3)
dev.off()

png("figures/figure1/Fig1C.png", width = 1500, height = 200)
par(lwd = 2.5, mar = rep(1,4))
barplot(t(t(all_i[,2]) / sum(all_i[,2])), col = c("gray80", mature_cols, "white"), las = 2, horiz = T, axes = F)
axis(side = 3)
dev.off()

png("figures/figure1/Fig1D.png", width = 1500, height = 200)
par(lwd = 2.5, mar = rep(1,4))
barplot(t(t(all_i[,3]) / sum(all_i[,3])), col = c("gray80", mature_cols, "white"), las = 2, horiz = T, axes = F)
axis(side = 3)
dev.off()

png("figures/figure1/Fig1_legend.png", width = 500, height = 1000)
plot(x = 1:10, y = 1:10, type = "n", axes = F, xlab = "", ylab = "")
legend("center", legend = rep("", length(markers) + 1), col = c("gray80", mature_cols), cex = 2.5, pch = 19, bty = "n")
dev.off()
############################
