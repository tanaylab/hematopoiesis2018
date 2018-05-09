################

dir.create("supp_figures/figureS1")
message("generation Supplementary Fig 1")

sc_raw_mat = sc_pipe_build_mat(scdb_init(basedir="saved_work/all_cells"))

B = unique(wells$Amp_batch_ID)

index_fn = "annotations/batch_stats.txt"
index = read.delim(index_fn, stringsAsFactors = F)
batch_stats = read.delim(index_fn, stringsAsFactors = F)
project_batches_stats = read.delim("results/TableS1.txt", stringsAsFactor = F)
batch_stats = merge(batch_stats[,c(1,2,3,4,5,6,11,12)], project_batches_stats[,c(2,3,6,7,8,9,10)], by.x = "amplification.batch", by.y = "amp_batch")
rownames(batch_stats) = batch_stats$amplification.batch
batch_stats = batch_stats[B,]
batch_ord = with(batch_stats, order(	factor(tier, levels = c(as.character(1:7), "Ly6g+", "d9", "d7")), 
                                        factor(treatment, levels = c("-", "Epo", "Epo (control)", "G-CSF", "G-CSF (control)", "CRISPR", "CRISPR_PU1", "CRISPR_CTL")),
                                        factor(virus.mix), factor(cytokines),
                                        mouse.specimen.unique.ID))

batch_stats = batch_stats[ batch_ord,]
batch_stats[ batch_stats == ""] = "-"
batch_stats$experiment = with(batch_stats, paste0(tier, ".",treatment, ".",virus.mix, ".", cytokines))
wells = with(wells, wells[ order(factor(Amp_batch_ID, levels = rownames(batch_stats))),])

cell_stats = sc_raw_mat@cell_metadata
wells$mouse = wells[rownames(wells), "Subject_ID"]
wells$experiment = batch_stats[wells$Amp_batch_ID, "experiment"]
wells$tier = cell_stats[rownames(wells), "tier"]

X = unique(wells[,c("Amp_batch_ID", "mouse")]); rownames(X) = X$Amp_batch_ID
batch_stats$mouse = X[ rownames(batch_stats), "mouse"]
tier_n = length(unique(wells$tier))
exp_n = length(unique(wells$experiment)); subject_n = length(unique(wells$mouse))
seq_n = length(unique(wells$Seq_batch_ID))

tier_cols = c("darkseagreen2", "darksalmon", "darkorchid1", "tomato", "rosybrown2", "orange", "khaki1", hsv(runif(tier_n - 7), runif(tier_n - 7, min = 0.5), runif(tier_n - 7, min = 0.5)))
subject_cols = hsv(runif(subject_n), runif(subject_n, min = 0.5), runif(subject_n, min = 0.5))
exp_cols = hsv(runif(exp_n), runif(exp_n, min = 0.5), runif(exp_n, min = 0.5))
seq_cols = hsv(runif(seq_n), runif(seq_n, min = 0.5), runif(seq_n, min = 0.5))

wells$col = tier_cols[ as.numeric(factor(wells$tier))]
batch_stats$col = tier_cols[ as.numeric(factor(batch_stats$tier))]
X = tapply(wells$filtered_in, factor(wells$Amp_batch_ID), mean)
batch_stats$ncells = X[ rownames(batch_stats)]
rownames(project_batches_stats) = project_batches_stats$amp_batch
batch_stats$noise.estimation = project_batches_stats[ rownames(batch_stats), "noise_estimation"]
ab_cls = with(wells, cumsum(table(factor(Amp_batch_ID, levels = unique(Amp_batch_ID)))))
exp_cls = with(wells, cumsum(table(factor(experiment, levels = unique(experiment)))))

png(paste0("supp_figures/figureS1/FigS1D-I.png"), width = 2000, height = 3000)
par(fig=c(0,1,0.6,0.8), new=TRUE, lwd = 2.5)
with(wells, plot(umicount, log = "y", xaxs = "i", axes = F, type = "n", xlab = "", ylab = ""))
abline(v = ab_cls, lty = 2, col = "darkgray")
abline(v = exp_cls, lty = 1, col = "black")
axis(2, las = 2)
with(wells, points(umicount, pch = 20, col = col, cex = 1.5))

par(fig=c(0,1,0.8,1), new=TRUE, lwd = 2.5)
with(wells, plot(readcount, log = "y", xaxs = "i", type = "n", xlab = "", ylab = "", axes = F))
abline(v = ab_cls, lty = 2, col = "darkgray")
axis(2, las = 2)
abline(v = exp_cls, lty = 1, col = "black")
with(wells, points(readcount, pch = 20, col = col, cex = 1.5))

par(fig=c(0,1,0.2,0.4), new=TRUE, lwd = 2.5)
with(batch_stats, plot(rowMeans(cbind(ab_cls,c(0,ab_cls[-length(ab_cls)]))), noise.estimation, pch = 18,
                                 xaxs = "i", xaxt = "n", type = "n", xlab = "", ylab = ""))
abline(v = ab_cls, lty = 2, col = "darkgray")
abline(v = exp_cls, lty = 1, col = "black")
with(batch_stats, points(rowMeans(cbind(ab_cls,c(0,ab_cls[-length(ab_cls)]))), noise.estimation, pch = 18,
                                   col = col, cex = 3))

par(fig=c(0,1,0.4,0.6), new=TRUE, lwd = 2.5)
with(batch_stats, plot(rowMeans(cbind(ab_cls,c(0,ab_cls[-length(ab_cls)]))), ncells, pch = 18,
                                 xaxs = "i", xaxt = "n", type = "n", xlab = "", ylab = ""))
abline(v = ab_cls, lty = 2, col = "darkgray")
abline(v = exp_cls, lty = 1, col = "black")
with(batch_stats, points(rowMeans(cbind(ab_cls,c(0,ab_cls[-length(ab_cls)]))), ncells, pch = 18,
                                   col = col, cex = 3))

par(fig=c(0,1,0.1,0.2), new=TRUE, lwd = 2.5)
with(wells, image(matrix(as.numeric(factor(Subject_ID, levels = unique(Subject_ID)))), col = subject_cols, axes = F))
par(fig=c(0,1,0,0.1), new=TRUE, lwd = 2.5)
with(wells, image(matrix(as.numeric(factor(Seq_batch_ID, levels = unique(Seq_batch_ID)))), col = seq_cols, axes = F))
dev.off()

Y = batch_stats[, c("tier", "treatment", "virus.mix", "cytokines", "mouse")]
X = merge(wells[,-c(14,15,16)], Y, by.x = "Amp_batch_ID", by.y = 0)
Y$id = with(Y, paste(tier, treatment, virus.mix, cytokines, sep = "."))
uY = unique(Y[,-5])
X$id = with(X, paste(tier, treatment, virus.mix, cytokines, sep = "."))
nbatches = table(Y$id)
nmice = rowSums(table(Y$id, Y$mouse) > 0)
Z = data.frame(ncells = table(X$id, X$filtered_in)[,2])
Z = merge(cbind(nbatches, nmice), Z, by.x = 0, by.y = 0)
rownames(Z) = Z[,1]; Z = Z[,-1]; colnames(Z)[1] = "nbatches"
Z = merge(uY, Z, by.x = "id", by.y = 0)
Z =  Z[order(  factor(Z$tier, levels = c(as.character(1:7), "Ly6g+", "d9", "d7")),
               factor(Z$treatment, levels = c("-", "Epo", "Epo (control)", "G-CSF", "G-CSF (control)", "CRISPR", "CRISPR_PU1", "CRISPR_CTL")),
               factor(Z$virus.mix), factor(Z$cytokines)),]
write.table(Z, sep = "\t", row.names = F, quote = F, file = "supp_figures/figureS1/figureS1_table.txt")
