#srun -p zhao --pty apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis:/home --home /home /lustre/fs4/zhao_lab/scratch/ulee/docka/garnett.sif R

library(garnett)
library(gridExtra)
library(ggplot2)
library(dplyr)

#meant to be run in apptainer image where scratch_store/analysis is mounted at /home/
path_mel <- "/home/sc_liftoff/dmel_testis_fb"
path_yak <- "/home/sc_liftoff/dyak_testis_fb"
path_ana <- "/home/sc_liftoff/dana_testis_fb"

cds_mel <- load_cellranger_data(path_mel)
cds_yak <- load_cellranger_data(path_yak)
cds_ana <- load_cellranger_data(path_ana)

mel_dat <- exprs(cds_mel)
yak_dat <- exprs(cds_yak)
ana_dat <- exprs(cds_ana)

orthos_yak <- read.csv("/home/sc_liftoff/genomes_fb/dmel_orthos_yak.tsv", skip=5, header=F, sep="\t")
orthos_ana <- read.csv("/home/sc_liftoff/genomes_fb/dmel_orthos_ana.tsv", skip=5, header=F, sep="\t")

genes_yak <- unique(orthos_yak[,2])
genes_yak_uniq <- genes_yak[unlist(lapply(genes_yak, function(x) sum(orthos_yak[,2] == x) == 1))]
genes_yak_yak <- orthos_yak[which(orthos_yak[,2] %in% genes_yak_uniq), 7]
genes_yak_uniq_yak <- genes_yak_yak[unlist(lapply(genes_yak_yak, function(x) sum(orthos_yak[,7] == x) == 1))]
genes_yak_one2one <- intersect(genes_yak_uniq, orthos_yak[which(orthos_yak[,7] %in% genes_yak_uniq_yak), 2])

genes_ana <- unique(orthos_ana[,2])
genes_ana_uniq <- genes_ana[unlist(lapply(genes_ana, function(x) sum(orthos_ana[,2] == x) == 1))]
genes_ana_ana <- orthos_ana[which(orthos_ana[,2] %in% genes_ana_uniq), 7]
genes_ana_uniq_ana <- genes_ana_ana[unlist(lapply(genes_ana_ana, function(x) sum(orthos_ana[,7] == x) == 1))]
genes_ana_one2one <- intersect(genes_ana_uniq, orthos_ana[which(orthos_ana[,7] %in% genes_ana_uniq_ana), 2])

orthos_yak_uniq <- orthos_yak[which(orthos_yak[,2] %in% genes_yak_one2one),]
orthos_ana_uniq <- orthos_ana[which(orthos_ana[,2] %in% genes_ana_one2one),]

genes_shared <- intersect(genes_yak_one2one, genes_ana_one2one)
genes_shared <- genes_shared[!grepl("Rp[A-Z]", genes_shared)]
genes_shared <- genes_shared[!grepl("mt:", genes_shared)]
genes_shared <- genes_shared[!grepl("tRNA", genes_shared)]

genes_shared_fbgn_mel <- orthos_ana_uniq[which(orthos_ana_uniq[,2] %in% genes_shared), 1]
genes_shared_fbgn_yak <- orthos_yak_uniq[which(orthos_yak_uniq[,2] %in% genes_shared), 6]
genes_shared_fbgn_ana <- orthos_ana_uniq[which(orthos_ana_uniq[,2] %in% genes_shared), 6]

mel_dat_shared <- mel_dat[genes_shared_fbgn_mel,]
yak_dat_shared <- yak_dat[genes_shared_fbgn_yak,]
ana_dat_shared <- ana_dat[genes_shared_fbgn_ana,]

rownames(mel_dat_shared) <- genes_shared
rownames(yak_dat_shared) <- genes_shared
rownames(ana_dat_shared) <- genes_shared

rd_share <- DataFrame(row.names=genes_shared)
rd_share$id <- genes_shared
rd_share$gene_short_name <- genes_shared

cds_mel_share <- new_cell_data_set(mel_dat_shared, cell_metadata=colData(cds_mel), gene_metadata=rd_share)
cds_yak_share <- new_cell_data_set(yak_dat_shared, cell_metadata=colData(cds_yak), gene_metadata=rd_share)
cds_ana_share <- new_cell_data_set(ana_dat_shared, cell_metadata=colData(cds_ana), gene_metadata=rd_share)

cds_comb <- combine_cds(list(cds_mel_share, cds_yak_share, cds_ana_share))

cds_spec <- as.character(pData(cds_comb)[,"sample"])
cds_spec[cds_spec == "1"] <- "mel"
cds_spec[cds_spec == "2"] <- "yak"
cds_spec[cds_spec == "3"] <- "ana"

pData(cds_comb)$spec <- cds_spec

cds_comb_proc <- preprocess_cds(cds_comb)
cds_comb_proc_batch <- align_cds(cds_comb_proc, alignment_group="spec")
cds_comb_proc_batch <- preprocess_cds(cds_comb_proc_batch)

cds_comb_proc <- reduce_dimension(cds_comb_proc, reduction_method="UMAP")
cds_comb_proc <- cluster_cells(cds_comb_proc)
cds_comb_proc_batch <- reduce_dimension(cds_comb_proc_batch, reduction_method="UMAP")
cds_comb_proc_batch <- cluster_cells(cds_comb_proc_batch)

dat_evan <- readRDS("/home/sc_liftoff/190730_R517_wild_integrated_SCtransform.RDS")

ev_idents <- dat_evan@active.ident
mel_idents <- colnames(exprs(cds_comb_proc_batch))
mel_idents_under <- unlist(lapply(mel_idents, function(x) gsub("-", "_", x)))
mel_idents_under <- unlist(lapply(mel_idents_under, function(x) substr(x, 1, nchar(x)-2)))
mel_intersect <- which(mel_idents_under %in% names(ev_idents))

mel_intersect_labs <- ev_idents[mel_intersect]
names(mel_intersect_labs) <- unlist(lapply(names(mel_intersect_labs), function(x) gsub("_", "-1_", x)))
mel_intersect_labs <- unlist(lapply(mel_intersect_labs, function(x) gsub(",", "", x)))

pData(cds_comb_proc_batch)[,"known_type"] <- rep("unassigned", nrow(pData(cds_comb_proc_batch)))
cds_mel_intersect <- which(names(mel_intersect_labs) %in% rownames(pData(cds_comb_proc_batch)))
pData(cds_comb_proc_batch)[names(mel_intersect_labs[cds_mel_intersect]), "known_type"] <- as.character(mel_intersect_labs[cds_mel_intersect])

# png("sc_liftoff/temp.png", width=1200, height=1200, res=300, pointsize=3)
# plot_cells(cds_comb_proc_batch, label_groups_by_cluster=FALSE,  color_cells_by = "spec")
# dev.off()

# #mel pseudotime
# cds_mel_share <- preprocess_cds(cds_mel_share, num_dim=50)
# cds_mel_share <- reduce_dimension(cds_mel_share, reduction_method="tSNE")
# cds_mel_share <- reduce_dimension(cds_mel_share)

# pData(cds_mel_share)[,"known_type"] <- pData(cds_comb_proc_batch)[1:5000,"known_type"]
# colData(cds_mel_share)[,"known_type"] <- pData(cds_comb_proc_batch)[1:5000,"known_type"]

# png("sc_liftoff/temp.png", width=1200, height=1200, res=300, pointsize=3)
# plot_cells(cds_mel_share)
# dev.off()

# cds_mel_share <- cluster_cells(cds_mel_share)

# png("sc_liftoff/temp1.png", width=1200, height=1200, res=300, pointsize=3)
# plot_cells(cds_mel_share, color_cells_by="partition")
# dev.off()

# cont_par <- list(c(1), c(1/3))
# names(cont_par) <- c("euclidean_distance_ratio", "geodesic_distance_ratio")

# cds_mel_share <- learn_graph(cds_mel_share, use_partition=F, learn_graph_control=cont_par)

# png("sc_liftoff/temp2.png", width=1200, height=1200, res=300, pointsize=3)
# plot_cells(cds_mel_share, reduction_method="UMAP",
           # color_cells_by = "known_type", group_cells_by="cluster", label_principal_points=T)
# dev.off()

# cds_mel_share <- order_cells(cds_mel_share, root_pr_nodes="Y_99")

# png("sc_liftoff/temp3.png", width=1200, height=1200, res=300, pointsize=3)
# plot_cells(cds_mel_share, 
           # color_cells_by = "pseudotime",
           # label_groups_by_cluster=FALSE,
           # label_leaves=FALSE,
		   # label_roots = FALSE,
           # label_branch_points=FALSE,
		   # graph_label_size=1.5)
# dev.off()

# cds_yak_share <- preprocess_cds(cds_yak_share, num_dim=50)
# cds_yak_share <- reduce_dimension(cds_yak_share)

# png("sc_liftoff/temp.png", width=1200, height=1200, res=300, pointsize=3)
# plot_cells(cds_yak_share)
# dev.off()

# cds_yak_share <- cluster_cells(cds_yak_share)

# png("sc_liftoff/temp1.png", width=1200, height=1200, res=300, pointsize=3)
# plot_cells(cds_yak_share, color_cells_by="partition")
# dev.off()

# cont_par <- list(c(1), c(1/3))
# names(cont_par) <- c("euclidean_distance_ratio", "geodesic_distance_ratio")

# cds_yak_share <- learn_graph(cds_yak_share, use_partition=F, learn_graph_control=cont_par)

# png("sc_liftoff/temp2.png", width=1200, height=1200, res=300, pointsize=3)
# plot_cells(cds_yak_share, reduction_method="UMAP",
           # color_cells_by = "cluster", group_cells_by="cluster", label_principal_points=T)
# dev.off()

# cds_yak_share <- order_cells(cds_yak_share, root_pr_nodes="Y_2")

# png("sc_liftoff/temp3.png", width=1200, height=1200, res=300, pointsize=3)
# plot_cells(cds_yak_share, 
           # color_cells_by = "pseudotime",
           # label_groups_by_cluster=FALSE,
           # label_leaves=FALSE,
		   # label_roots = FALSE,
           # label_branch_points=FALSE,
		   # graph_label_size=1.5)
# dev.off()

# cds_ana_share <- preprocess_cds(cds_ana_share, num_dim=50)
# cds_ana_share <- reduce_dimension(cds_ana_share)

# png("sc_liftoff/temp.png", width=1200, height=1200, res=300, pointsize=3)
# plot_cells(cds_ana_share)
# dev.off()

# cds_ana_share <- cluster_cells(cds_ana_share)

# png("sc_liftoff/temp1.png", width=1200, height=1200, res=300, pointsize=3)
# plot_cells(cds_ana_share, color_cells_by="partition")
# dev.off()

# cont_par <- list(c(1), c(1/3))
# names(cont_par) <- c("euclidean_distance_ratio", "geodesic_distance_ratio")

# cds_ana_share <- learn_graph(cds_ana_share, use_partition=F, learn_graph_control=cont_par)

# png("sc_liftoff/temp2.png", width=1200, height=1200, res=300, pointsize=3)
# plot_cells(cds_ana_share, reduction_method="UMAP",
           # color_cells_by = "cluster", group_cells_by="cluster", label_principal_points=T)
# dev.off()

# cds_ana_share <- order_cells(cds_ana_share, root_pr_nodes="Y_75")

# png("sc_liftoff/temp3.png", width=1200, height=1200, res=300, pointsize=3)
# plot_cells(cds_ana_share, 
           # color_cells_by = "pseudotime",
           # label_groups_by_cluster=FALSE,
           # label_leaves=FALSE,
		   # label_roots = FALSE,
           # label_branch_points=FALSE,
		   # graph_label_size=1.5)
# dev.off()

#find genes with low between species var, high between cluster var in mel
spec_dat1 <- pData(cds_comb_proc)[,"spec"]
type_dat1 <- pData(cds_comb_proc_batch)[,"known_type"]

expr_mat1 <- exprs(cds_comb_proc)

# expr_mat_mel <- expr_mat[,which(spec_dat == "mel")]
# expr_mat_mel <- rbind(as.data.frame(as.matrix(expr_mat_mel)), type_dat[which(spec_dat == "mel")])
# rownames(expr_mat_mel) <- c(rownames(expr_mat), "known_type")

# expr_mat <- rbind(as.data.frame(as.matrix(expr_mat)), spec_dat)
# rownames(expr_mat) <- c(rownames(exprs(cds_comb_proc)), "spec")

# anova_stat <- matrix(0, nrow=nrow(exprs(cds_comb_proc)), ncol=2)
# rownames(anova_stat) <- rownames(exprs(cds_comb_proc))
# colnames(anova_stat) <- c("spec", "known_type")

# for(i in 1:nrow(exprs(cds_comb_proc))){
	# print(i)
	# anova_stat[i, 1] <- oneway.test(as.numeric(expr_mat[i,]) ~ as.character(expr_mat["spec",]), var.equal=F)$statistic
	# anova_stat[i, 2] <- oneway.test(as.numeric(expr_mat_mel[i,]) ~ as.character(expr_mat_mel["known_type",]), var.equal=F)$statistic
# }

#as a function
anova_function <- function(cds_in, spec_data, type_data, expr_matrix, target_spec_label = "mel", type_label = "known_type", spec_label="spec", verbose=FALSE, shuffle=FALSE){
	expr_matrix_mel <- expr_matrix[,which(spec_data == target_spec_label)]
	expr_matrix_mel <- rbind(as.data.frame(as.matrix(expr_matrix_mel)), type_data[which(spec_data == target_spec_label)])
	rownames(expr_matrix_mel) <- c(rownames(expr_matrix), type_label)

	expr_matrix <- rbind(as.data.frame(as.matrix(expr_matrix)), spec_data)
	rownames(expr_matrix) <- c(rownames(exprs(cds_in)), spec_label)

	anova_statistics <- matrix(0, nrow=nrow(exprs(cds_in)), ncol=2)
	rownames(anova_statistics) <- rownames(exprs(cds_in))
	colnames(anova_statistics) <- c(spec_label, type_label)

	if(shuffle){
		expr_matrix[spec_label,] <- sample(expr_matrix[spec_label,], replace=F)
		expr_matrix_mel[type_label,] <- sample(expr_matrix_mel[type_label,], replace=F)
	}

	for(i in 1:nrow(exprs(cds_in))){
		if(verbose) print(i)
		anova_statistics[i, 1] <- oneway.test(as.numeric(expr_matrix[i,]) ~ as.character(expr_matrix[spec_label,]), var.equal=F)$statistic
		anova_statistics[i, 2] <- oneway.test(as.numeric(expr_matrix_mel[i,]) ~ as.character(expr_matrix_mel[type_label,]), var.equal=F)$statistic
	}

	return(anova_statistics)
}

anova_stat <- anova_function(cds_comb_proc, spec_dat1, type_dat1, expr_mat1, "mel", "known_type", "spec", TRUE, FALSE)

write.csv(anova_stat, file="/home/sc_liftoff/anova_stat.csv")
#anova_stat <- read.csv("/home/sc_liftoff/anova_stat.csv", row.names="X")


idx <- intersect(which(!is.na(anova_stat[,1])), which(!is.na(anova_stat[,2])))
anova_cor <- cor(log(anova_stat[idx,1]), log(anova_stat[idx,2]))
anova_cor <- round(anova_cor, digits=3)

gene_list <- intersect(which(anova_stat[,1] < median(anova_stat[,1], na.rm=T)), which(anova_stat[,2] > median(anova_stat[,2], na.rm=T)))
gene_list <- rownames(anova_stat)[gene_list]

png("sc_liftoff/anova_plot.png", width=1200, height=1200, res=300, pointsize=5)
anova_main <- bquote("ANOVA log(F-Statistic) by Gene, " ~ rho == .(anova_cor))
plot(log(anova_stat[idx,2]), log(anova_stat[idx,1]), main=anova_main, xlab="", ylab="", cex.axis=1.75, cex.main=2)
mtext("Predictive of Cell Type", side=1, line=3.2, cex=2)
mtext("Predictive of Species", side=2, line=2.7, cex=2)
abline(v=log(median(anova_stat[,2], na.rm=T)))
abline(h=log(median(anova_stat[,1], na.rm=T)))
dev.off()

png("sc_liftoff/anova_plot_scram.png", width=1200, height=1200, res=300, pointsize=5)
anova_main <- "Scrambled ANOVA log(F-Statistic) by Gene"
plot(log(sample(anova_stat[idx,2], replace=F)), log(sample(anova_stat[idx,1], replace=F)), main=anova_main, xlab="", ylab="", cex.axis=1.75, cex.main=2)
mtext("Predictive of Cell Type", side=1, line=3.2, cex=2)
mtext("Predictive of Species", side=2, line=2.7, cex=2)
abline(v=median(log(anova_stat[,2]), na.rm=T))
abline(h=median(log(anova_stat[,1]), na.rm=T))
dev.off()


svg("sc_liftoff/anova_plot.svg", width=4, height=4, pointsize=5)
anova_main <- bquote("ANOVA log(F-Statistic) by Gene, " ~ rho == .(anova_cor))
plot(log(anova_stat[idx,2]), log(anova_stat[idx,1]), main=anova_main, xlab="", ylab="", cex.axis=1.75, cex.main=2)
mtext("Predictive of Cell Type", side=1, line=3.2, cex=2)
mtext("Predictive of Species", side=2, line=2.7, cex=2)
abline(v=log(median(anova_stat[,2], na.rm=T)))
abline(h=log(median(anova_stat[,1], na.rm=T)))
dev.off()

svg("sc_liftoff/anova_plot_scram.svg", width=4, height=4, pointsize=5)
anova_main <- "Scrambled ANOVA log(F-Statistic) by Gene"
plot(log(sample(anova_stat[idx,2], replace=F)), log(sample(anova_stat[idx,1], replace=F)), main=anova_main, xlab="", ylab="", cex.axis=1.75, cex.main=2)
mtext("Predictive of Cell Type", side=1, line=3.2, cex=2)
mtext("Predictive of Species", side=2, line=2.7, cex=2)
abline(v=median(log(anova_stat[,2]), na.rm=T))
abline(h=median(log(anova_stat[,1]), na.rm=T))
dev.off()


write.table(gene_list, "sc_liftoff/reduced_gene_list_monocle.txt", row.names=F, col.names=F, quote=F)
write.table(genes_shared_fbgn_mel[genes_shared %in% gene_list], "sc_liftoff/reduced_gene_list_monocle_fbgn.txt", row.names=F, col.names=F, quote=F)


cds_comb_diffgene <- new_cell_data_set(exprs(cds_comb)[gene_list,], cell_metadata=colData(cds_comb_proc_batch), gene_metadata=rd_share[gene_list,])
cds_comb_diffgene <- preprocess_cds(cds_comb_diffgene)
cds_comb_diffgene <- align_cds(cds_comb_diffgene, alignment_group="spec")
cds_comb_diffgene <- reduce_dimension(cds_comb_diffgene)
cds_comb_diffgene <- cluster_cells(cds_comb_diffgene, resolution=3e-4)

clu_umap <- cds_comb_diffgene@clusters$UMAP$clusters
pData(cds_comb_diffgene)[,"clu_umap"] <- as.character(clu_umap)
colData(cds_comb_diffgene)[,"clu_umap"] <- as.character(clu_umap)

spec_lab <- colData(cds_comb_diffgene)$spec
type_dat <- colData(cds_comb_diffgene)$known_type

cds_comb_diffgene_mel <- cds_comb_diffgene[,!type_dat == "unassigned"]
cds_comb_diffgene_yak <- cds_comb_diffgene[,spec_lab == "yak"]
cds_comb_diffgene_ana <- cds_comb_diffgene[,spec_lab == "ana"]

png("sc_liftoff/genelist_clusters_1.png", width=2400*2, height=2400*2, res=600, pointsize=3)
p1 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="known_type", group_cells_by="known_type", group_label_size = 3) + ggtitle("Known Annotation (D. melanogaster)")
p2 <- plot_cells(cds_comb_diffgene_mel, group_label_size = 3) + ggtitle("D. melanogaster")
p3 <- plot_cells(cds_comb_diffgene_yak, group_label_size = 3) + ggtitle("D. yakuba")
p4 <- plot_cells(cds_comb_diffgene_ana, group_label_size = 3) + ggtitle("D. ananassae")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

png("sc_liftoff/genelist_clusters_pca_1.png", width=2400*2, height=2400*2, res=600, pointsize=3)
p1 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="known_type", group_cells_by="known_type", reduction_method="PCA", group_label_size = 3) + ggtitle("Known Annotation (D. melanogaster)")
p2 <- plot_cells(cds_comb_diffgene_mel, reduction_method="PCA", color_cells_by="clu_umap", group_label_size = 3) + ggtitle("D. melanogaster")
p3 <- plot_cells(cds_comb_diffgene_yak, reduction_method="PCA", color_cells_by="clu_umap", group_label_size = 3) + ggtitle("D. yakuba")
p4 <- plot_cells(cds_comb_diffgene_ana, reduction_method="PCA", color_cells_by="clu_umap", group_label_size = 3) + ggtitle("D. ananassae")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

png("sc_liftoff/genelist_clusters_1_nolab.png", width=2400*2, height=2400*2, res=600, pointsize=3)
p1 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="known_type", group_cells_by="known_type", group_label_size = 0) + ggtitle("Known Annotation (D. melanogaster)")
p2 <- plot_cells(cds_comb_diffgene_mel, group_label_size = 0) + ggtitle("D. melanogaster")
p3 <- plot_cells(cds_comb_diffgene_yak, group_label_size = 0) + ggtitle("D. yakuba")
p4 <- plot_cells(cds_comb_diffgene_ana, group_label_size = 0) + ggtitle("D. ananassae")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

png("sc_liftoff/genelist_clusters_pca_1_nolab.png", width=2400*2, height=2400*2, res=600, pointsize=3)
p1 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="known_type", group_cells_by="known_type", reduction_method="PCA", group_label_size = 0) + ggtitle("Known Annotation (D. melanogaster)")
p2 <- plot_cells(cds_comb_diffgene_mel, reduction_method="PCA", color_cells_by="clu_umap", group_label_size = 0) + ggtitle("D. melanogaster")
p3 <- plot_cells(cds_comb_diffgene_yak, reduction_method="PCA", color_cells_by="clu_umap", group_label_size = 0) + ggtitle("D. yakuba")
p4 <- plot_cells(cds_comb_diffgene_ana, reduction_method="PCA", color_cells_by="clu_umap", group_label_size = 0) + ggtitle("D. ananassae")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()


svg("sc_liftoff/genelist_clusters_1.svg", width=8, height=8, pointsize=3)
p1 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="known_type", group_cells_by="known_type", group_label_size = 3) + ggtitle("Known Annotation (D. melanogaster)")
p2 <- plot_cells(cds_comb_diffgene_mel, group_label_size = 3) + ggtitle("D. melanogaster")
p3 <- plot_cells(cds_comb_diffgene_yak, group_label_size = 3) + ggtitle("D. yakuba")
p4 <- plot_cells(cds_comb_diffgene_ana, group_label_size = 3) + ggtitle("D. ananassae")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

svg("sc_liftoff/genelist_clusters_pca_1.svg", width=8, height=8, pointsize=3)
p1 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="known_type", group_cells_by="known_type", reduction_method="PCA", group_label_size = 3) + ggtitle("Known Annotation (D. melanogaster)")
p2 <- plot_cells(cds_comb_diffgene_mel, reduction_method="PCA", color_cells_by="clu_umap", group_label_size = 3) + ggtitle("D. melanogaster")
p3 <- plot_cells(cds_comb_diffgene_yak, reduction_method="PCA", color_cells_by="clu_umap", group_label_size = 3) + ggtitle("D. yakuba")
p4 <- plot_cells(cds_comb_diffgene_ana, reduction_method="PCA", color_cells_by="clu_umap", group_label_size = 3) + ggtitle("D. ananassae")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

svg("sc_liftoff/genelist_clusters_1_nolab.svg", width=8, height=8, pointsize=3)
p1 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="known_type", group_cells_by="known_type", group_label_size = 0) + ggtitle("Known Annotation (D. melanogaster)")
p2 <- plot_cells(cds_comb_diffgene_mel, group_label_size = 0) + ggtitle("D. melanogaster")
p3 <- plot_cells(cds_comb_diffgene_yak, group_label_size = 0) + ggtitle("D. yakuba")
p4 <- plot_cells(cds_comb_diffgene_ana, group_label_size = 0) + ggtitle("D. ananassae")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

svg("sc_liftoff/genelist_clusters_pca_1_nolab.svg", width=8, height=8, pointsize=3)
p1 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="known_type", group_cells_by="known_type", reduction_method="PCA", group_label_size = 0) + ggtitle("Known Annotation (D. melanogaster)")
p2 <- plot_cells(cds_comb_diffgene_mel, reduction_method="PCA", color_cells_by="clu_umap", group_label_size = 0) + ggtitle("D. melanogaster")
p3 <- plot_cells(cds_comb_diffgene_yak, reduction_method="PCA", color_cells_by="clu_umap", group_label_size = 0) + ggtitle("D. yakuba")
p4 <- plot_cells(cds_comb_diffgene_ana, reduction_method="PCA", color_cells_by="clu_umap", group_label_size = 0) + ggtitle("D. ananassae")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

#cluster analysis
n_reclu <- 10
reclu_lab_fine <- c("Late spermatocyte", "Early spermatocyte I", "Early spermatid II", "Early spermatocyte II", "Early spermatid I", "Late spermatogonia", "Somatic", "Late spermatid", "ananassae spermatocyte", "GSC/Early spermatogonia")
reclu_lab_coarse <- c("Late spermatocyte", "Early spermatocyte", "Early spermatid", "Early spermatocyte", "Early spermatid", "Late spermatogonia", "Somatic", "Late spermatid", "ananassae spermatocyte", "GSC/Early spermatogonia") 

reclu_fine <- as.character(clu_umap)
reclu_coarse <- as.character(clu_umap)

for(i in 1:n_reclu){
	reclu_fine[reclu_fine == as.character(i)] <- reclu_lab_fine[i]
	reclu_coarse[reclu_coarse == as.character(i)] <- reclu_lab_coarse[i]
}

pData(cds_comb_diffgene)[,"reclu_fine"] <- reclu_fine
colData(cds_comb_diffgene)[,"reclu_fine"] <- reclu_fine

pData(cds_comb_diffgene)[,"reclu_coarse"] <- reclu_coarse
colData(cds_comb_diffgene)[,"reclu_coarse"] <- reclu_coarse

cds_comb_diffgene_mel <- cds_comb_diffgene[,!type_dat == "unassigned"]
cds_comb_diffgene_yak <- cds_comb_diffgene[,spec_lab == "yak"]
cds_comb_diffgene_ana <- cds_comb_diffgene[,spec_lab == "ana"]


png("sc_liftoff/genelist_clusters_reclu.png", width=2400*2, height=2400*2, res=600, pointsize=3)
p1 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="known_type", group_cells_by="known_type", group_label_size = 3) + ggtitle("Known Annotation (D. melanogaster)")
p2 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="reclu_coarse", group_cells_by="reclu_coarse", group_label_size = 3) + ggtitle("D. melanogaster")
p3 <- plot_cells(cds_comb_diffgene_yak, color_cells_by="reclu_coarse", group_cells_by="reclu_coarse", group_label_size = 3) + ggtitle("D. yakuba")
p4 <- plot_cells(cds_comb_diffgene_ana, color_cells_by="reclu_coarse", group_cells_by="reclu_coarse", group_label_size = 3) + ggtitle("D. ananassae")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

png("sc_liftoff/genelist_clusters_pca_reclu.png", width=2400*2, height=2400*2, res=600, pointsize=3)
p1 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="known_type", group_cells_by="known_type", reduction_method="PCA", group_label_size = 3) + ggtitle("Known Annotation (D. melanogaster)")
p2 <- plot_cells(cds_comb_diffgene_mel, reduction_method="PCA", color_cells_by="reclu_coarse", group_label_size = 3) + ggtitle("D. melanogaster")
p3 <- plot_cells(cds_comb_diffgene_yak, reduction_method="PCA", color_cells_by="reclu_coarse", group_label_size = 3) + ggtitle("D. yakuba")
p4 <- plot_cells(cds_comb_diffgene_ana, reduction_method="PCA", color_cells_by="reclu_coarse", group_label_size = 3) + ggtitle("D. ananassae")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

png("sc_liftoff/genelist_clusters_reclu_nolab.png", width=2400*2, height=2400*2, res=600, pointsize=3)
p1 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="known_type", group_cells_by="known_type", group_label_size = 0) + ggtitle("Known Annotation (D. melanogaster)")
p2 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="reclu_coarse", group_cells_by="reclu_coarse", group_label_size = 0) + ggtitle("D. melanogaster")
p3 <- plot_cells(cds_comb_diffgene_yak, color_cells_by="reclu_coarse", group_cells_by="reclu_coarse", group_label_size = 0) + ggtitle("D. yakuba")
p4 <- plot_cells(cds_comb_diffgene_ana, color_cells_by="reclu_coarse", group_cells_by="reclu_coarse", group_label_size = 0) + ggtitle("D. ananassae")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

png("sc_liftoff/genelist_clusters_pca_reclu_nolab.png", width=2400*2, height=2400*2, res=600, pointsize=3)
p1 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="known_type", group_cells_by="known_type", reduction_method="PCA", group_label_size = 0) + ggtitle("Known Annotation (D. melanogaster)")
p2 <- plot_cells(cds_comb_diffgene_mel, reduction_method="PCA", color_cells_by="reclu_coarse", group_label_size = 0) + ggtitle("D. melanogaster")
p3 <- plot_cells(cds_comb_diffgene_yak, reduction_method="PCA", color_cells_by="reclu_coarse", group_label_size = 0) + ggtitle("D. yakuba")
p4 <- plot_cells(cds_comb_diffgene_ana, reduction_method="PCA", color_cells_by="reclu_coarse", group_label_size = 0) + ggtitle("D. ananassae")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()



svg("sc_liftoff/genelist_clusters_reclu.svg", width=8, height=8, pointsize=3)
p1 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="known_type", group_cells_by="known_type", group_label_size = 3) + ggtitle("Known Annotation (D. melanogaster)")
p2 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="reclu_coarse", group_cells_by="reclu_coarse", group_label_size = 3) + ggtitle("D. melanogaster")
p3 <- plot_cells(cds_comb_diffgene_yak, color_cells_by="reclu_coarse", group_cells_by="reclu_coarse", group_label_size = 3) + ggtitle("D. yakuba")
p4 <- plot_cells(cds_comb_diffgene_ana, color_cells_by="reclu_coarse", group_cells_by="reclu_coarse", group_label_size = 3) + ggtitle("D. ananassae")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

svg("sc_liftoff/genelist_clusters_pca_reclu.svg", width=8, height=8, pointsize=3)
p1 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="known_type", group_cells_by="known_type", reduction_method="PCA", group_label_size = 3) + ggtitle("Known Annotation (D. melanogaster)")
p2 <- plot_cells(cds_comb_diffgene_mel, reduction_method="PCA", color_cells_by="reclu_coarse", group_label_size = 3) + ggtitle("D. melanogaster")
p3 <- plot_cells(cds_comb_diffgene_yak, reduction_method="PCA", color_cells_by="reclu_coarse", group_label_size = 3) + ggtitle("D. yakuba")
p4 <- plot_cells(cds_comb_diffgene_ana, reduction_method="PCA", color_cells_by="reclu_coarse", group_label_size = 3) + ggtitle("D. ananassae")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

svg("sc_liftoff/genelist_clusters_reclu_nolab.svg", width=8, height=8, pointsize=3)
p1 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="known_type", group_cells_by="known_type", group_label_size = 0) + ggtitle("Known Annotation (D. melanogaster)")
p2 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="reclu_coarse", group_cells_by="reclu_coarse", group_label_size = 0) + ggtitle("D. melanogaster")
p3 <- plot_cells(cds_comb_diffgene_yak, color_cells_by="reclu_coarse", group_cells_by="reclu_coarse", group_label_size = 0) + ggtitle("D. yakuba")
p4 <- plot_cells(cds_comb_diffgene_ana, color_cells_by="reclu_coarse", group_cells_by="reclu_coarse", group_label_size = 0) + ggtitle("D. ananassae")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

svg("sc_liftoff/genelist_clusters_pca_reclu_nolab.svg", width=8, height=8, pointsize=3)
p1 <- plot_cells(cds_comb_diffgene_mel, color_cells_by="known_type", group_cells_by="known_type", reduction_method="PCA", group_label_size = 0) + ggtitle("Known Annotation (D. melanogaster)")
p2 <- plot_cells(cds_comb_diffgene_mel, reduction_method="PCA", color_cells_by="reclu_coarse", group_label_size = 0) + ggtitle("D. melanogaster")
p3 <- plot_cells(cds_comb_diffgene_yak, reduction_method="PCA", color_cells_by="reclu_coarse", group_label_size = 0) + ggtitle("D. yakuba")
p4 <- plot_cells(cds_comb_diffgene_ana, reduction_method="PCA", color_cells_by="reclu_coarse", group_label_size = 0) + ggtitle("D. ananassae")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()



pData(cds_comb_proc)[,"reclu_fine"] <- reclu_fine
colData(cds_comb_proc)[,"reclu_fine"] <- reclu_fine

pData(cds_comb_proc)[,"reclu_coarse"] <- reclu_coarse
colData(cds_comb_proc)[,"reclu_coarse"] <- reclu_coarse

marker_test_diffgene <- top_markers(cds_comb_diffgene, group_cells_by="reclu_coarse")
marker_test_proc <- top_markers(cds_comb_proc, group_cells_by="reclu_coarse")

top_specific_markers_diffgene <- marker_test_diffgene %>%
							filter(fraction_expressing >= 0.10) %>%
							group_by(cell_group) %>%
							top_n(2, pseudo_R2)

top_specific_marker_ids_diffgene <- unique(top_specific_markers_diffgene %>% pull(gene_id))

top_specific_markers_proc <- marker_test_proc %>%
							filter(fraction_expressing >= 0.10) %>%
							group_by(cell_group) %>%
							top_n(2, pseudo_R2)

top_specific_marker_ids_proc <- unique(top_specific_markers_proc %>% pull(gene_id))


png("sc_liftoff/markers_comb_genelist.png", width=1200*3, height=1200*1.75, res=300)
p <- plot_genes_by_group(cds_comb_diffgene, top_specific_marker_ids_diffgene, group_cells_by="reclu_coarse", ordering_type="maximal_on_diag", max.size=10)
p$data$Group <- factor(p$data$Group, levels=c("Somatic", "GSC/Early spermatogonia", "Late spermatogonia", "ananassae spermatocyte", "Early spermatocyte", "Late spermatocyte", "Early spermatid",  "Late spermatid"))
p
dev.off()

png("sc_liftoff/markers_comb_all.png", width=1200*3, height=1200*1.75, res=300)
p <- plot_genes_by_group(cds_comb_proc, top_specific_marker_ids_proc, group_cells_by="reclu_coarse", ordering_type="maximal_on_diag", max.size=10)
p$data$Group <- factor(p$data$Group, levels=c("Somatic", "GSC/Early spermatogonia", "Late spermatogonia", "ananassae spermatocyte", "Early spermatocyte", "Late spermatocyte", "Early spermatid",  "Late spermatid"))
p
dev.off()


svg("sc_liftoff/markers_comb_genelist.svg", width=4*3, height=4*1.75)
p <- plot_genes_by_group(cds_comb_diffgene, top_specific_marker_ids_diffgene, group_cells_by="reclu_coarse", ordering_type="maximal_on_diag", max.size=10)
p$data$Group <- factor(p$data$Group, levels=c("Somatic", "GSC/Early spermatogonia", "Late spermatogonia", "ananassae spermatocyte", "Early spermatocyte", "Late spermatocyte", "Early spermatid",  "Late spermatid"))
p
dev.off()

svg("sc_liftoff/markers_comb_all.svg", width=4*3, height=4*1.75)
p <- plot_genes_by_group(cds_comb_proc, top_specific_marker_ids_proc, group_cells_by="reclu_coarse", ordering_type="maximal_on_diag", max.size=10)
p$data$Group <- factor(p$data$Group, levels=c("Somatic", "GSC/Early spermatogonia", "Late spermatogonia", "ananassae spermatocyte", "Early spermatocyte", "Late spermatocyte", "Early spermatid",  "Late spermatid"))
p
dev.off()


# png("sc_liftoff/markers_comb_genelist_vln.png", width=1200*2.25, height=1200*2.25, res=300)
# p <- plot_genes_violin(cds_comb_diffgene[top_specific_marker_ids_diffgene,], group_cells_by="reclu_coarse", ncol=2) +
      # theme(axis.text.x=element_text(angle=45, hjust=1))
# p$data$Group <- factor(p$data$Group, levels=c("Somatic", "GSC/Early spermatogonia", "Late spermatogonia", "Early spermatocyte", "Late spermatocyte", "Early spermatid",  "Late spermatid", "ananassae spermatocyte"))
# p
# dev.off()

# png("sc_liftoff/markers_comb_all_vln.png", width=1200*2.25, height=1200*2.25, res=300)
# p <- plot_genes_violin(cds_comb_proc[top_specific_marker_ids_proc,], group_cells_by="reclu_coarse", ncol=2) +
      # theme(axis.text.x=element_text(angle=45, hjust=1))
# p$data$Group <- factor(p$data$Group, levels=c("Somatic", "GSC/Early spermatogonia", "Late spermatogonia", "Early spermatocyte", "Late spermatocyte", "Early spermatid",  "Late spermatid", "ananassae spermatocyte"))
# p
# dev.off()


#write out by cell type for txburst
for(i in levels(as.factor(reclu_lab_coarse))){
	for(j in levels(as.factor(cds_spec))){
		out_file <- paste("sc_liftoff/bursting/", j, "_", gsub("[/ ]", "_", i), ".csv", sep="")
		out_idx <- intersect(which(colData(cds_comb_diffgene)[,"reclu_coarse"] == i), which(colData(cds_comb_diffgene)[,"spec"] == j))
		
		bar_codes <- colData(cds_comb_diffgene)[out_idx,"barcode"]
		
		if(j == "ana") write.csv(as.matrix(exprs(cds_ana[,bar_codes])), file=out_file)
		if(j == "mel") write.csv(as.matrix(exprs(cds_mel[,bar_codes])), file=out_file)
		if(j == "yak") write.csv(as.matrix(exprs(cds_yak[,bar_codes])), file=out_file)
		
	}
}

write.csv(colData(cds_comb_diffgene)[,c("spec", "known_type", "reclu_fine", "reclu_coarse")], file="sc_liftoff/ct_assignments.csv")

save(list=ls(), file="sc_liftoff/wrkspce_run_integrated.rdata")

#analyze txburst

#diff exprs by cluster
#use cell type idents into seurat

cds_comb_diffgene <- learn_graph(cds_comb_diffgene)
cds_comb_diffgene <- order_cells(cds_comb_diffgene)

png("sc_liftoff/genelist_pseudotime.png", width=2400, height=2400, res=600, pointsize=3)
plot_cells(cds_comb_diffgene, color_cells_by="reclu_coarse", label_groups_by_cluster=F, label_leaves=F, label_branch_points=F)
dev.off()

#use old umap and full data from cds_comb_proc
cds_comb_diffgene_all <- cds_comb_diffgene

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

col_wheel <- gg_color_hue(8)

ct_cols <- col_wheel[c(7, 3, 6, 2, 5, 1, 4, 8)]

plot_single_gene <- function(in_cds, gene_name, y_limit=2, m_lab){
	gg_color_hue <- function(n) {
	  hues = seq(15, 375, length = n + 1)
	  hcl(h = hues, l = 65, c = 100)[1:n]
	}

	col_wheel <- gg_color_hue(8)

	ct_cols <- col_wheel[c(7, 3, 6, 2, 5, 1, 4, 8)]

	p1 <- plot_cells(in_cds, genes=c(gene_name), label_cell_groups=F, alpha=0.7) + ggtitle(paste0(m_lab, ": ", gene_name))
	p1 <- p1 + theme(text=element_text(size=7.5), plot.title=element_text(size=12)) + scale_color_gradient2(low="skyblue1", mid="navyblue", high="red3", midpoint=65)
	p2 <- plot_genes_violin(in_cds[gene_name,], group_cells_by="reclu_coarse", ncol=8, log_scale=F) 
	p2$data$reclu_coarse <- factor(p2$data$reclu_coarse, levels=c("Somatic", "GSC/Early spermatogonia", "Late spermatogonia", "ananassae spermatocyte", "Early spermatocyte", "Late spermatocyte", "Early spermatid",  "Late spermatid"))
	# p2 <- p2 + scale_fill_manual(values=ct_cols) + coord_cartesian(ylim=c(0, y_limit)) + theme(text=element_text(size=8), axis.text.x=element_text(angle=25, hjust=1))
	p2 <- p2 + scale_fill_manual(values=ct_cols) + theme(text=element_text(size=8), axis.text.x=element_text(angle=25, hjust=1))


	return(list(p1, p2))
}

plot_single_gene_spec <- function(all_cds, gene_name, y_limits = c(2, 2, 2, 1)){
	sp_lab <- colData(all_cds)$spec
	ty_dat <- colData(all_cds)$known_type

	all_cds_mel <- all_cds[,!ty_dat == "unassigned"]
	all_cds_yak <- all_cds[,sp_lab == "yak"]
	all_cds_ana <- all_cds[,sp_lab == "ana"]

	pl_list <- list()
	pl_list[[1]] <- plot_single_gene(all_cds, gene_name, y_limits[1], "Combined")
	pl_list[[2]] <- plot_single_gene(all_cds_mel, gene_name, y_limits[2], "melanogaster")
	pl_list[[3]] <- plot_single_gene(all_cds_yak, gene_name, y_limits[3], "yakuba")
	pl_list[[4]] <- plot_single_gene(all_cds_ana, gene_name, y_limits[4], "ananassae")

	names(pl_list) <- c("all", "mel", "yak", "ana")

	return(pl_list)
}


write_out_fig <- function(in_cds, out_lab, gene_nam){
	out_label <- paste0(out_lab, gene_nam, ".png") 
	png(out_label, width=2400*4, height=2400*1.4, res=600, pointsize=3)
	p_list_all <- plot_single_gene_spec(in_cds, gene_nam, c(2, 2, 2.25, 1))

	grid.arrange(p_list_all$all[[1]], p_list_all$mel[[1]], p_list_all$yak[[1]], p_list_all$ana[[1]],
		p_list_all$all[[2]], p_list_all$mel[[2]], p_list_all$yak[[2]], p_list_all$ana[[2]],
		ncol=4, heights=c(1.75,1))

	dev.off()
}

write_out_fig(cds_comb_diffgene_all, "sc_liftoff/all_figs/fig_", "aly")

lapply(rownames(cds_comb_diffgene_all), function(x) write_out_fig(cds_comb_diffgene_all, "sc_liftoff/all_figs/fig_", x))


plot_single_gene_add <- function(in_cds, full_cds, gene_name, y_limit=2, m_lab){
	gg_color_hue <- function(n) {
	  hues = seq(15, 375, length = n + 1)
	  hcl(h = hues, l = 65, c = 100)[1:n]
	}

	col_wheel <- gg_color_hue(8)

	ct_cols <- col_wheel[c(7, 3, 6, 2, 5, 1, 4, 8)]
	
	in_cds_local <- in_cds

	if(!gene_name %in% rownames(exprs(in_cds))){
		in_cds_local@assays@data@listData$counts[1,] <- exprs(full_cds)[gene_name,]
		rownames(in_cds_local@assays@data@listData$counts)[1] <- gene_name
		names(in_cds_local@rowRanges)[1] <- gene_name
		in_cds_local@rowRanges@elementMetadata@listData$gene_short_name[1] <- gene_name
	}

	p1 <- plot_cells(in_cds_local, genes=c(gene_name), label_cell_groups=F, alpha=0.7) + ggtitle(paste0(m_lab, ": ", gene_name))
	p1 <- p1 + theme(text=element_text(size=7.5), plot.title=element_text(size=12)) + scale_color_gradient2(low="skyblue1", mid="navyblue", high="red3", midpoint=65)
	p2 <- plot_genes_violin(in_cds_local[gene_name,], group_cells_by="reclu_coarse", ncol=8, log_scale=F) 
	p2$data$reclu_coarse <- factor(p2$data$reclu_coarse, levels=c("Somatic", "GSC/Early spermatogonia", "Late spermatogonia", "ananassae spermatocyte", "Early spermatocyte", "Late spermatocyte", "Early spermatid",  "Late spermatid"))
	# p2 <- p2 + scale_fill_manual(values=ct_cols) + coord_cartesian(ylim=c(0, y_limit)) + theme(text=element_text(size=8), axis.text.x=element_text(angle=25, hjust=1))
	p2 <- p2 + scale_fill_manual(values=ct_cols) + theme(text=element_text(size=8), axis.text.x=element_text(angle=25, hjust=1))


	return(list(p1, p2))
}

plot_single_gene_spec_add <- function(all_cds, full_cds, gene_name, y_limits = c(2, 2, 2, 1)){
	sp_lab <- colData(all_cds)$spec
	ty_dat <- colData(all_cds)$known_type

	all_cds_mel <- all_cds[,!ty_dat == "unassigned"]
	all_cds_yak <- all_cds[,sp_lab == "yak"]
	all_cds_ana <- all_cds[,sp_lab == "ana"]
	
	full_cds_mel <- full_cds[,!ty_dat == "unassigned"]
	full_cds_yak <- full_cds[,sp_lab == "yak"]
	full_cds_ana <- full_cds[,sp_lab == "ana"]	

	pl_list <- list()
	pl_list[[1]] <- plot_single_gene_add(all_cds, full_cds, gene_name, y_limits[1], "Combined")
	pl_list[[2]] <- plot_single_gene_add(all_cds_mel, full_cds_mel, gene_name, y_limits[2], "melanogaster")
	pl_list[[3]] <- plot_single_gene_add(all_cds_yak, full_cds_yak, gene_name, y_limits[3], "yakuba")
	pl_list[[4]] <- plot_single_gene_add(all_cds_ana, full_cds_ana, gene_name, y_limits[4], "ananassae")

	names(pl_list) <- c("all", "mel", "yak", "ana")

	return(pl_list)
}

write_out_fig_add <- function(in_cds, full_cds, out_lab, gene_nam){
	out_label <- paste0(out_lab, gene_nam, ".png")
	png(out_label, width=2400*4, height=2400*1.4, res=600, pointsize=3)
	p_list_all <- plot_single_gene_spec_add(in_cds, full_cds, gene_nam, c(2, 2, 2.25, 1))

	grid.arrange(p_list_all$all[[1]], p_list_all$mel[[1]], p_list_all$yak[[1]], p_list_all$ana[[1]],
		p_list_all$all[[2]], p_list_all$mel[[2]], p_list_all$yak[[2]], p_list_all$ana[[2]],
		ncol=4, heights=c(1.75,1))

	dev.off()
}


raz_genes <- c("aly", "aub", "bam", "CadN", "can", "cher", "CycB", "Dic61B", "dlg1", "esg", "eya", "f-cup", "Fas3", "fzo", "hh", "Hml", "kl-2", "kl-3", "kl-5", "kmg", "Mst77F", "Mst84Db", "Mst84Dc", "Mst87F", "MtnA", "Nep5", "p53", "p-cup", "piwi", "Rbp4", "sa", "shg", "so", "soti", "stg", "Syt1", "tj", "tomboy20", "upd1", "vas", "wa-cup", "zpg")

#kl-3, Mst84Db, Mst84Dc, Mst87F, not in orthologs list entirely
#kl-5 not in yakuba
#fzo, kl-5, vas not in ana
#sa duplicated in yakuba

raz_genes <- subset(raz_genes, raz_genes %in% rownames(exprs(cds_comb)))

#aly and wa-cup already in 198 genes

raz_genes <- subset(raz_genes, !raz_genes %in% rownames(exprs(cds_comb_diffgene_all)))

symbs_mel <- read.csv("/home/sc_liftoff/genomes_fb/dmel_orthologs_in_drosophila_species_fb_2022_01.tsv", skip=5, header=F, sep="\t")


lapply(raz_genes, function(x) write_out_fig_add(cds_comb_diffgene_all, cds_comb, "sc_liftoff/all_figs_raz/fig_", x))


pData(cds_comb) <- cbind(pData(cds_comb), pData(cds_comb_diffgene_all)[,"reclu_coarse"])
colnames(pData(cds_comb))[5] <- "reclu_coarse"

genes_ana_clu <- read.csv("sc_liftoff/an_top_10.tsv", sep="\t", header=T)

genes_ana_clu <- subset(genes_ana_clu, gene  %in% rownames(exprs(cds_comb))) 

lapply(genes_ana_clu[,"gene"], function(x) write_out_fig_add(cds_comb_diffgene_all, cds_comb, "sc_liftoff/all_figs_an_top20/fig_", x))

hcr_genes <- c("Rbp4", "B52", "Pkd2", "soti")

lapply(hcr_genes, function(x) write_out_fig_add(cds_comb_diffgene_all, cds_comb, "sc_liftoff/hcr_marks/fig_", x))

redo_genes <- c("Obp99a", "eEF1alpha1", "CG34434", "CG15876", "CAH13", "schuy", "whip", "Pzl")

lapply(redo_genes, function(x) write_out_fig_add(cds_comb_diffgene_all, cds_comb, "sc_liftoff/all_figs_redo/fig_", x))


save(list=ls(), file="sc_liftoff/wrkspce_run_integrated.rdata")


#similarity matrix
#first make averages by gene 
#pairwise 

reclu_coarse_fac <- factor(as.factor(reclu_coarse), levels=c("Somatic", "GSC/Early spermatogonia", "Late spermatogonia", "ananassae spermatocyte", "Early spermatocyte", "Late spermatocyte", "Early spermatid",  "Late spermatid"))
pData(cds_comb)[,"reclu_coarse"] <- reclu_coarse_fac
colData(cds_comb)[,"reclu_coarse"] <- reclu_coarse_fac
colData(cds_comb)$known_type <- colData(cds_comb_diffgene_all)$known_type

spec <- colData(cds_comb)$spec
type1 <- colData(cds_comb)$known_type

cor_analysis_mel <- cds_comb[,!type1 == "unassigned"]
cor_analysis_yak <- cds_comb[,spec == "yak"]
cor_analysis_ana <- cds_comb[,spec == "ana"]

ct_avgs <- matrix(0, nrow=nrow(exprs(cds_comb)), ncol=8*3)

ct_avgs_list <- lapply(list(cor_analysis_mel, cor_analysis_yak, cor_analysis_ana), function(x){
		return(lapply(levels(reclu_coarse_fac), function(y){
			return(rowMeans(exprs(x)[,which(colData(x)$reclu_coarse == y)]))
		}))
	})

ct_avgs <- do.call("cbind", ct_avgs_list)
ct_avgs <- do.call("cbind", ct_avgs)

colnames(ct_avgs) <- unlist(lapply(c("mel", "yak", "ana"), function(x) unlist(lapply(levels(reclu_coarse_fac), function(y) paste0(x, "_", y)))))
rownames(ct_avgs) <- rownames(exprs(cds_comb))

ct_cors <- matrix(0, nrow=8*3, ncol=8*3)
colnames(ct_cors) <- colnames(ct_avgs)
rownames(ct_cors) <- colnames(ct_cors)

for(i in 1:24){
	for(j in 1:24){
		ct_cors[i,j] <- cor(ct_avgs[,rownames(ct_cors)[i]], ct_avgs[,colnames(ct_cors)[j]], method="spearman")
	}
}

library(RColorBrewer)

png("sc_liftoff/ct_pair_cors.png", width=3600, height=3600, res=300)
heatmap(ct_cors[24:1,], Rowv=NA, Colv=NA, margins=c(13,13), scale="none", col= colorRampPalette(brewer.pal(8, "Oranges"))(25))
legend(x="bottomright", legend=c(paste0("cor=", as.character(signif(min(as.numeric(ct_cors)), digits=3))), paste("cor=", as.character(signif(median(as.numeric(ct_cors)), digits=2))), "cor=1.00"), 
     fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))
dev.off()



png("test1.png", width=2400*2, height=2400*1.4, res=600, pointsize=3)
temp1 <- plot_single_gene_add(cds_comb_diffgene_all, cds_comb, "stg", m_lab="Combined: stg")
grid.arrange(temp1[[1]], temp4[[1]], temp1[[2]], temp4[[2]], ncol=2, heights=c(1.75, 1))
dev.off()

pr_graph_test_res <- graph_test(cds_comb_proc, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))





# # cluster/garnett classifier
# ev_class <- train_cell_classifier(cds=cds_comb_diffgene, marker_file="/home/sc_liftoff/ev_markfree.txt", db="none", num_unknown=1000)
# feat_genes_comb_proc <- get_feature_genes(ev_class, node="root", db="none")

# cds_comb_diffgene_pred <- classify_cells(cds_comb_diffgene, ev_class, db='none', cluster_extend=T)

# assgnd <- which(!pData(cds_comb_diffgene_pred)[,"cell_type"] == "Unknown")
# sum(pData(cds_comb_diffgene_pred)[assgnd,"cluster_ext_type"] == pData(cds_comb_diffgene_pred)[assgnd,"cell_type"])

# temp1 <- pData(cds_comb_diffgene)[,"sample"]
# temp1[which(temp1 == 1)] <- "mel"
# temp1[which(temp1 == 2)] <- "yak"
# temp1[which(temp1 == 3)] <- "ana"

# temp <- cds_comb_diffgene
# pData(temp)[,"spec"] <- as.factor(temp1)

# temp <- align_cds(temp, alignment_group="spec")
# temp <- reduce_dimension(temp)
# png("sc_liftoff/clu_fb.png", width=1200, height=1200, res=300, pointsize=3)
# plot_cells(temp, color_cells_by="spec")
# dev.off()

