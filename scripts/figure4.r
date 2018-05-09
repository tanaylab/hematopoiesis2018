###################
# Figure 4
##################
message("generating figure 4")
dir.create("figures/figure4")

igh_nm = grep("Igh", rownames(umis), value = T)
genes = c('Ifitm1','Zfpm1','Lmo4','Cd52','Cd34','Flt3','H2afy','Irf8','Sfpi1','Mpo','Serpina3f;Serpina3g','Hlf','Apoe','Car1','Gata2','Eltd1','Gpr56',
	'Klf1', 'Cebpa', 'Gfi1', 'Satb1', 'Myc', 'Dntt', 'Gata1', 'Prtn3', igh_nm)
b = 9
mat = umis_n
dir.create("figures/figure4/FigD-E")
palette = c("white", "cornsilk1", "orange","red3", "purple4", "midnightblue")
for (val in genes) {
	vals = mat[val,]
	norm_val = rep(1, length(vals)); names(norm_val) = names(vals)
	norm_val[ vals != 0] = as.numeric(cut(vals[ vals != 0], unique(quantile(vals[ vals != 0], (0:b)/b)), include.lowest = T)) + 1;
	cols = colorRampPalette(palette)(max(norm_val))
	png(paste0("figures/figure4/FigD-E/", gsub(";.*", "", val), ".png"), height = 1200, width = 1200)
	plot(sc_2d@x[good_cells], sc_2d@y[good_cells], cex = 1 + 0.7 * round((norm_val[good_cells] - 1) / max(norm_val) * 5),
        	pch = 21, bg = cols[norm_val[good_cells]], axes = F, xlab = "", ylab = "")
	lines(outline[,1], outline[,2], lwd = 4)
	dev.off()
}

png("figures/figure4/FigD-E/gene_cb.png", height = 200, width = 1000)
par(mar = rep(0,4))
image(matrix(1:100), col = colorRampPalette(palette)(100), axes = F)
dev.off()

###################

cell_stats = sc_cl@scmat@cell_metadata
bad_cells = rownames(cell_stats)[ cell_stats$tier == 3 | cell_stats$amplification.batch %in% c("AB3975", "AB3976")]
early_cells = setdiff(names(s_score)[ s_score > quantile(s_score, 0.6)], bad_cells)
early_sep = as.numeric(cut(s_score[ early_cells], quantile(s_score[ early_cells], c(0,4,7,8)/8), include.lowest = T)); names(early_sep) = early_cells
s_cols = rev(c("#e49a15", "#f7bb13","#f9e20a"))
igh_nm = grep("IGH", rownames(umis), value = T)
must_haves = c("Apoe", "Dntt", "Gata2", "Car2", igh_nm, "Car1", "Flt3", "Mpo", "Prtn3")
us = umis[,early_cells]
dus = .downsamp(us,1200) 
early_cells = colnames(dus)
rep_genes = setdiff(names(which(rowSums(us[,early_cells]) > 100)), c(ribo_genes, cc_genes))
C3 = cor(t(log2(1 + dus[unique(c(rep_genes, must_haves)), ])), s_score[early_cells], method = "spearman")
rep_genes = union(must_haves, rownames(C3)[ C3 < 0])
C1 = cor(log(1 + t(dus[rep_genes, names(which(early_sep[early_cells] == 1))])), method = "spearman"); diag(C1) = NA
genes = union(must_haves, names(head(sort(apply(C1,1,function(x) sort(x)[3])),45)))
C2 = C1[genes, genes]
zlim = max(max(C2,na.rm = T), -min(C2,na.rm = T))
#dimnames(C2) = list(gsub(";.*","", rownames(C2)), gsub(";.*","", colnames(C2)))
hc = hclust(dist(C2), method = "ward.D2")
ord_vec = C2["Dntt",] - C2["Apoe",]; ord_vec["Dntt"] = 1; ord_vec["Apoe"] = -1
hc = as.hclust(reorder(as.dendrogram(hc),
                             ord_vec,
                             agglo.FUN=max))
shades = colorRampPalette(c("midnightblue", "darkblue", "blue3", "cornflowerblue", "white", "darkgoldenrod1", "chocolate2", "darkorange3", "firebrick4"))(1000)
message("correlation zlim: ", zlim[1])

dir.create("figures/figure4/Fig4A-C")
write.table(gsub(";.*","", rev(genes[hc$order])), row.names = F, quote = F, col.names = F, file = "figures/figure4/Fig4A-C/markers.txt")
titles = c("C", "B", "A")
for (i in 1:3) {
  png(paste0("figures/figure4/Fig4A-C/", titles[i], "_cells.png"), width = 1000, height = 1000)
  plot(sc_2d@x[good_cells], sc_2d@y[good_cells], xlab = "", ylab = "", axes = F, type = "n")
  points(sc_2d@x[good_cells], sc_2d@y[good_cells], pch = 20, col = "gray40", cex = 3 * (good_cells %in% names(which(early_sep == i))))
  points(outline[,1], outline[,2], type = "l", lwd = 7)
  dev.off()
  C = cor(t(log2(1 + dus[genes, names(which(early_sep[early_cells] == i))])), method = "spearman"); diag(C) = NA
  png(paste0("figures/figure4/Fig4A-C/", titles[i], ".png"), width = 1500, height = 1500)
  par(mar = rep(5,4))
  image(pmin(pmax(C[hc$order, hc$order], -zlim), zlim), col = shades, axes = F, zlim = c(-zlim,zlim))
  mtext(gsub(";.*", "", rownames(C)[hc$order]), side = 2,las = 2, at = (1 - seq_len(nrow(C))) / (1 - nrow(C)))
  dev.off()
}

png("figures/figure4/Fig4A-C/colorbar.png", width = 1000, height = 500)
par(mar = rep(0,4))
image(matrix(1:1000), col = shades, axes = F)
dev.off()
