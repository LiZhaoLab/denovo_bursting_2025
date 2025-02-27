#srun -p zhao --pty apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis:/home --home /home /lustre/fs4/zhao_lab/scratch/ulee/docka/garnett.sif R
library(garnett)
library(ggplot2)

setwd("sc_liftoff")

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

gen_bootstrap <- function(n_gene, idx_in, idx_out, dist_in, dist_out){
	res_list <- lapply(1:n_gene, function(x){
			res <- rep(0, length(c(idx_in, idx_out)))
			if(length(idx_in) > 0) res[idx_in] <- sample(dist_in, length(idx_in), replace=T)
			res[idx_out] <- sample(dist_out, length(idx_out), replace=T)
			return(res)
		})
		
	res_mtx <- do.call("rbind", res_list)
	return(res_mtx)
}

gen_ref_mtx <- function(lab_vec, gene_names, dist_in, dist_out){
	exp_mat <- matrix(sample(dist_out, length(gene_names)*length(lab_vec), replace=T), nrow=length(gene_names), ncol=length(lab_vec))
	rownames(exp_mat) <- gene_names

	label_names <- levels(as.factor(as.character(lab_vec)))

	idx_genes <- 1:round(length(gene_names)/length(label_names), digits=0)

	#make_ref_mtx_sp
	for(i in label_names){
		idx_ct_in <- which(lab_vec == i)
		idx_ct_out <- unique(sort(c(which(!lab_vec == i), which(is.na(lab_vec)))))
		if(is.na(i)) {
			idx_ct_in <- c()
			idx_ct_out <- 1:length(idx_genes)
		}

		j <- which(label_names == i) - 1
		idx_genes_this <- idx_genes + j*length(idx_genes)
		print(idx_genes_this)
		exp_mat[idx_genes_this,] <- gen_bootstrap(length(idx_genes), idx_ct_in, idx_ct_out, dist_in, dist_out)
	}
	
	return(exp_mat)
}

gen_ref_mtx_hack <- function(lab_vec, gene_names, dist_in, dist_out){
	exp_mat <- matrix(sample(dist_out, length(gene_names)*length(lab_vec), replace=T), nrow=length(gene_names), ncol=length(lab_vec))
	rownames(exp_mat) <- gene_names

	label_names <- levels(as.factor(as.character(lab_vec)))

	idx_genes <- 1:floor(length(gene_names)/length(label_names))

	idx_ct_in_sp1 <- 1:1800
	idx_ct_out_sp1 <- 1801:2700
	idx_genes_sp1 <- 1:6666
	
	idx_ct_in_sp3 <- 1801:2700
	idx_ct_out_sp3 <- 1:1800
	idx_genes_sp3 <- 6667:9999
	
	exp_mat[idx_genes_sp1,] <- gen_bootstrap(6666, idx_ct_in_sp1, idx_ct_out_sp1, dist_in, dist_out)
	exp_mat[idx_genes_sp3,] <- gen_bootstrap(3333, idx_ct_in_sp3, idx_ct_out_sp3, dist_in, dist_out)
	
	return(exp_mat)
}


calc_cor <- function(anova_stat){
	idx <- intersect(which(!is.na(anova_stat[,1])), which(!is.na(anova_stat[,2])))
	idx <- intersect(idx, which(anova_stat[,1] > 0))
	idx <- intersect(idx, which(anova_stat[,2] > 0))
	anova_cor <- cor(log(anova_stat[idx,1]), log(anova_stat[idx,2]))
	anova_cor <- round(anova_cor, digits=3)
	return(anova_cor)
}

gen_genelist <- function(anov_dat, n_top, dir_x, dir_y){
	dat_ord <- cbind(rank(dir_x*anov_dat[,1]), rank(dir_y*anov_dat[,2]))
	list_max <- apply(dat_ord, 1, function(x) return(min(c(x[1], x[2]))))
	list_ord <- order(list_max, decreasing=TRUE)
	
	tops <- rownames(anov_dat)[list_ord[1:n_top]]
	
	print(paste("Median overall x value: ", median(anov_dat[,1]), sep=""))
	print(paste("Median overall y value: ", median(anov_dat[,2]), sep=""))
	
	print(paste("Median subset x value: ", median(anov_dat[tops,1]), sep=""))
	print(paste("Median subset y value: ", median(anov_dat[tops,2]), sep=""))
	
	return(tops)
}



#------------------------main-------------------

#load in real data
path_mel <- "/home/sc_liftoff/dmel_testis_fb"
cds_mel <- load_cellranger_data(path_mel)
mel_dat <- exprs(cds_mel)

exp_means <- rowMeans2(mel_dat)
names(exp_means) <- rownames(mel_dat)

exp_means <- exp_means[order(exp_means, decreasing=T)]


percentile_top <- 90
percentile_bot <- 80
percentile_evenhigher <- 95
idx_top <- round(length(exp_means)*(1-percentile_top/100), digits=0)
idx_bot <- round(length(exp_means)*(1-percentile_bot/100), digits=0)
idx_evenhigher <- round(length(exp_means)*(1-percentile_evenhigher/100), digits=0)

high_exprs <- mel_dat[names(exp_means)[idx_top],]

low_exprs <- mel_dat[names(exp_means)[idx_bot],]

evenhigher_exprs <- mel_dat[names(exp_means)[idx_evenhigher],]

#make simulated data set, one with tissue specific genes being same, one with idents changed
species <- c("sp1", "sp2", "sp3")
cell_types <- unlist(lapply(1:9, function(x) paste("ct_", x, sep="")))

n_genes <- 9999
n_cell_per_spec <- 900

gene_labs <- unlist(lapply(1:n_genes, function(x) paste("gene_", x, sep="")))
cell_labs <- unlist(lapply(1:(n_cell_per_spec * length(species)), function(x) paste("cell_", x, sep="")))

n_types <- length(cell_types)
n_cell_per_type <- n_cell_per_spec/n_types
n_spec <- length(species)
n_cell_tot <- n_spec * n_cell_per_spec

species_lab <- unlist(lapply(species, function(x) rep(x, n_cell_per_spec)))
ct_lab <- rep(unlist(lapply(cell_types, function(x) rep(x, n_cell_per_type))), n_spec)
ct_lab_asym <- ct_lab

species_lab_asym <- species_lab
species_lab_asym[c(201:300, 501:600, 201:300 + 900, 501:600 + 900)] <- NA

species_lab_asym2 <- species_lab
species_lab_asym2[which(species_lab == "sp2")] <- "sp1"

#make reference and mixed expression matrices
ref_mat_sp <- gen_ref_mtx(species_lab, gene_labs, high_exprs, low_exprs)
ref_mat_sp_asym <- gen_ref_mtx(species_lab_asym, gene_labs, high_exprs, low_exprs)
# ref_mat_sp_asym[2222:3332,] <- gen_bootstrap(1111, c(), 1:2700, high_exprs, low_exprs)
# ref_mat_sp_asym[5555:6665,] <- gen_bootstrap(1111, c(), 1:2700, high_exprs, low_exprs)

ref_mat_sp_asym2 <- gen_ref_mtx_hack(species_lab_asym2, gene_labs, high_exprs, low_exprs)
ref_mat_ct <- gen_ref_mtx(ct_lab, gene_labs, evenhigher_exprs, low_exprs)
ref_mat_ct_asym <- gen_ref_mtx(ct_lab_asym, gene_labs, evenhigher_exprs, low_exprs)

# ct_lab1 <- unlist(lapply(1:9, function(x) rep(x, n_cell_per_spec/3)))
# ref_mat_ct <- gen_ref_mtx(ct_lab1, gene_labs, high_exprs, low_exprs)

p_mix <- 0.5
exp_mat_all <- (1-p_mix)*ref_mat_sp + p_mix*ref_mat_ct

beta_shape <- 7
beta_vals <- rbeta(nrow(ref_mat_sp), beta_shape, beta_shape)

exp_mat_all <- lapply(1:nrow(ref_mat_sp), function(x){
		# pr <- runif(1)
		pr <- beta_vals[x]
		return(pr*ref_mat_sp[x,] + (1-pr)*ref_mat_ct[x,])
	})

exp_mat_all <- do.call("rbind", exp_mat_all)
colnames(exp_mat_all) <- cell_labs
rownames(exp_mat_all) <- rownames(ref_mat_sp)	


exp_mat_all_asym <- lapply(1:nrow(ref_mat_sp), function(x){
		pr <- beta_vals[x]
		return(pr*ref_mat_sp_asym[x,] + (1-pr)*ref_mat_ct_asym[x,])
	})

exp_mat_all_asym <- do.call("rbind", exp_mat_all_asym)
colnames(exp_mat_all_asym) <- cell_labs
rownames(exp_mat_all_asym) <- rownames(ref_mat_sp)	

exp_mat_all_asym2 <- lapply(1:nrow(ref_mat_sp), function(x){
		pr <- beta_vals[x]
		return(pr*ref_mat_sp_asym2[x,] + (1-pr)*ref_mat_ct[x,])
	})

exp_mat_all_asym2 <- do.call("rbind", exp_mat_all_asym2)
colnames(exp_mat_all_asym2) <- cell_labs
rownames(exp_mat_all_asym2) <- rownames(ref_mat_sp)	


ref_mat_noise <- matrix(sample(c(high_exprs, low_exprs), nrow(ref_mat_sp)*ncol(ref_mat_sp), replace=T), nrow=nrow(ref_mat_sp))
colnames(ref_mat_noise) <- cell_labs
rownames(ref_mat_noise) <- rownames(ref_mat_sp)

exp_mat_cor <- lapply(1:nrow(exp_mat_all), function(x){
		pr <- runif(1)
		return(pr*exp_mat_all[x,] + (1-pr)*ref_mat_noise[x,])
	})

exp_mat_cor <- do.call("rbind", exp_mat_cor)
colnames(exp_mat_cor) <- colnames(exp_mat_all)
rownames(exp_mat_cor) <- rownames(exp_mat_all)	
	
exp_mat_all_asym <- lapply(1:nrow(exp_mat_all), function(x){
		pr <- runif(1)
		return(pr*exp_mat_all_asym[x,] + (1-pr)*ref_mat_noise[x,])
	})
	
exp_mat_all_asym <- do.call("rbind", exp_mat_all_asym)
colnames(exp_mat_all_asym) <- colnames(exp_mat_all)
rownames(exp_mat_all_asym) <- rownames(exp_mat_all)

exp_mat_all_asym2 <- lapply(1:nrow(exp_mat_all), function(x){
		pr <- runif(1)
		return(pr*exp_mat_all_asym2[x,] + (1-pr)*ref_mat_noise[x,])
	})
	
exp_mat_all_asym2 <- do.call("rbind", exp_mat_all_asym2)
colnames(exp_mat_all_asym2) <- colnames(exp_mat_all)
rownames(exp_mat_all_asym2) <- rownames(exp_mat_all)


cdat <- cbind(species_lab, ct_lab)
rownames(cdat) <- 1:n_cell_tot

sim_dat_sp <- new_cell_data_set(ref_mat_sp, cell_metadata=cdat)
sim_dat_ct <- new_cell_data_set(ref_mat_ct, cell_metadata=cdat)
sim_dat_comb <- new_cell_data_set(exp_mat_all, cell_metadata=cdat)
sim_dat_cor <- new_cell_data_set(exp_mat_cor, cell_metadata=cdat)
sim_dat_cor_asym <- new_cell_data_set(exp_mat_all_asym, cell_metadata=cdat)
sim_dat_cor_asym2 <- new_cell_data_set(exp_mat_all_asym2, cell_metadata=cdat)

sim_dat_sp <- cluster_cells(reduce_dimension(preprocess_cds(sim_dat_sp), reduction_method="UMAP"))
sim_dat_ct <- cluster_cells(reduce_dimension(preprocess_cds(sim_dat_ct), reduction_method="UMAP"))
sim_dat_comb <- cluster_cells(reduce_dimension(preprocess_cds(sim_dat_comb), reduction_method="UMAP"))
sim_dat_cor <- cluster_cells(reduce_dimension(preprocess_cds(sim_dat_cor), reduction_method="UMAP"))
sim_dat_cor_asym <- cluster_cells(reduce_dimension(preprocess_cds(sim_dat_cor_asym), reduction_method="UMAP"))
sim_dat_cor_asym2 <- cluster_cells(reduce_dimension(preprocess_cds(sim_dat_cor_asym2), reduction_method="UMAP"))


sim_dat_cor_batch <- cluster_cells(reduce_dimension(align_cds(preprocess_cds(sim_dat_cor), alignment_group="species_lab"), reduction_method="UMAP"))


anova_stat_sp <- anova_function(sim_dat_sp, species_lab, ct_lab, ref_mat_sp, "sp1", "species_lab", "ct_lab", TRUE, FALSE)
anova_stat_ct <- anova_function(sim_dat_ct, species_lab, ct_lab, ref_mat_ct, "sp1", "species_lab", "ct_lab", TRUE, FALSE)
anova_stat_comb <- anova_function(sim_dat_comb, species_lab, ct_lab, exp_mat_all, "sp1", "species_lab", "ct_lab", TRUE, FALSE)
anova_stat_cor <- anova_function(sim_dat_cor, species_lab, ct_lab, exp_mat_cor, "sp1", "species_lab", "ct_lab", TRUE, FALSE)
anova_stat_cor_asym <- anova_function(sim_dat_cor_asym, species_lab, ct_lab, exp_mat_all_asym, "sp1", "species_lab", "ct_lab", TRUE, FALSE)
anova_stat_cor_asym2 <- anova_function(sim_dat_cor_asym2, species_lab, ct_lab, exp_mat_all_asym2, "sp1", "species_lab", "ct_lab", TRUE, FALSE)


an_cor_sp <- calc_cor(anova_stat_sp)
an_cor_ct <- calc_cor(anova_stat_ct)
an_cor_comb <- calc_cor(anova_stat_comb)
an_cor_cor <- calc_cor(anova_stat_cor)
an_cor_cor_asym <- calc_cor(anova_stat_cor_asym)
an_cor_cor_asym2 <- calc_cor(anova_stat_cor_asym2)

top_list <- gen_genelist(cbind(anova_stat_cor[,2], anova_stat_cor[,1]), 200, 1, -1)
top_list_asym <- gen_genelist(cbind(anova_stat_cor_asym[,2], anova_stat_cor_asym[,1]), 200, 1, -1)
top_list_asym2 <- gen_genelist(cbind(anova_stat_cor_asym2[,2], anova_stat_cor_asym2[,1]), 200, 1, -1)



sim_dat_cor_fix <- new_cell_data_set(exp_mat_cor[top_list,], cell_metadata=cdat)
sim_dat_cor_fix <- cluster_cells(reduce_dimension(preprocess_cds(sim_dat_cor_fix), reduction_method="UMAP"))

sim_dat_asym_fix <- new_cell_data_set(exp_mat_all_asym[top_list_asym,], cell_metadata=cdat) 
sim_dat_asym_fix <- cluster_cells(reduce_dimension(preprocess_cds(sim_dat_asym_fix), reduction_method="UMAP"))

sim_dat_asym2_fix <- new_cell_data_set(exp_mat_all_asym2[top_list_asym2,], cell_metadata=cdat) 
sim_dat_asym2_fix <- cluster_cells(reduce_dimension(preprocess_cds(sim_dat_asym2_fix), reduction_method="UMAP"))

idx_sp <- intersect(!is.na(log(anova_stat_sp[,1])), !is.na(log(anova_stat_sp[,2]))) 
idx_ct <- intersect(!is.na(log(anova_stat_ct[,1])), !is.na(log(anova_stat_ct[,2])))
idx_comb <- intersect(!is.na(log(anova_stat_comb[,1])), !is.na(log(anova_stat_comb[,2])))
idx_cor <- intersect(!is.na(log(anova_stat_cor[,1])), !is.na(log(anova_stat_cor[,2])))
idx_asym <- intersect(!is.na(log(anova_stat_cor_asym[,1])), !is.na(log(anova_stat_cor_asym[,2])))
idx_asym2 <- intersect(!is.na(log(anova_stat_cor_asym2[,1])), !is.na(log(anova_stat_cor_asym2[,2])))



#plot simulated set
png("sim_anova_sp.png", width=1200, height=1200, res=300, pointsize=5)
anova_main <- bquote("Species: ANOVA log(F-Statistic) by Gene, " ~ rho == .(an_cor_sp))
plot(log(anova_stat_sp[idx_sp,2]), log(anova_stat_sp[idx_sp,1]), main=anova_main, xlab="", ylab="", cex.axis=1.75, cex.main=2)
mtext("Predictive of Cell Type", side=1, line=3.2, cex=2)
mtext("Predictive of Species", side=2, line=2.7, cex=2)
abline(v=median(log(anova_stat_sp[idx_sp,2]), na.rm=T))
abline(h=median(log(anova_stat_sp[idx_sp,1]), na.rm=T))
dev.off()

png("sim_anova_ct.png", width=1200, height=1200, res=300, pointsize=5)
anova_main <- bquote("Cell Type: ANOVA log(F-Statistic) by Gene, " ~ rho == .(an_cor_ct))
plot(log(anova_stat_ct[idx_ct,2]), log(anova_stat_ct[idx_ct,1]), main=anova_main, xlab="", ylab="", cex.axis=1.75, cex.main=2)
mtext("Predictive of Cell Type", side=1, line=3.2, cex=2)
mtext("Predictive of Species", side=2, line=2.7, cex=2)
abline(v=median(log(anova_stat_ct[idx_ct,2]), na.rm=T))
abline(h=median(log(anova_stat_ct[idx_ct,1]), na.rm=T))
dev.off()

png("sim_anova_comb.png", width=1200, height=1200, res=300, pointsize=5)
anova_main <- bquote("Combined: ANOVA log(F-Statistic) by Gene, " ~ rho == .(an_cor_comb))
plot(log(anova_stat_comb[idx_comb,2]), log(anova_stat_comb[idx_comb,1]), main=anova_main, xlab="", ylab="", cex.axis=1.75, cex.main=2)
mtext("Predictive of Cell Type", side=1, line=3.2, cex=2)
mtext("Predictive of Species", side=2, line=2.7, cex=2)
abline(v=median(log(anova_stat_comb[idx_comb,2]), na.rm=T))
abline(h=median(log(anova_stat_comb[idx_comb,1]), na.rm=T))
dev.off()

png("sim_anova_cor.png", width=1200, height=1200, res=300, pointsize=5)
anova_main <- bquote("Simulated: ANOVA log(F-Statistic) by Gene, " ~ rho == .(an_cor_cor))
plot(log(anova_stat_cor[idx_cor,2]), log(anova_stat_cor[idx_cor,1]), main=anova_main, xlab="", ylab="", cex.axis=1.75, cex.main=2)
mtext("Predictive of Cell Type", side=1, line=3.2, cex=2)
mtext("Predictive of Species", side=2, line=2.7, cex=2)
abline(v=log(median(anova_stat_cor[idx_cor,2], na.rm=T)))
abline(h=log(median(anova_stat_cor[idx_cor,1], na.rm=T)))
abline(v=log(min(anova_stat_cor[top_list, 2])), lty=2)
abline(h=log(max(anova_stat_cor[top_list, 1])), lty=2)
dev.off()

png("sim_anova_asym.png", width=1200, height=1200, res=300, pointsize=5)
anova_main <- bquote("Simulated: ANOVA log(F-Statistic) by Gene, " ~ rho == .(an_cor_cor_asym))
plot(log(anova_stat_cor_asym[idx_asym,2]), log(anova_stat_cor_asym[idx_asym,1]), main=anova_main, xlab="", ylab="", cex.axis=1.75, cex.main=2)
mtext("Predictive of Cell Type", side=1, line=3.2, cex=2)
mtext("Predictive of Species", side=2, line=2.7, cex=2)
abline(v=log(median(anova_stat_cor_asym[idx_asym,2], na.rm=T)))
abline(h=log(median(anova_stat_cor_asym[idx_asym,1], na.rm=T)))
abline(v=log(min(anova_stat_cor_asym[top_list_asym, 2])), lty=2)
abline(h=log(max(anova_stat_cor_asym[top_list_asym, 1])), lty=2)
dev.off()

png("sim_anova_asym2.png", width=1200, height=1200, res=300, pointsize=5)
anova_main <- bquote("Simulated: ANOVA log(F-Statistic) by Gene, " ~ rho == .(an_cor_cor_asym2))
plot(log(anova_stat_cor_asym2[idx_asym2,2]), log(anova_stat_cor_asym2[idx_asym2,1]), main=anova_main, xlab="", ylab="", cex.axis=1.75, cex.main=2)
mtext("Predictive of Cell Type", side=1, line=3.2, cex=2)
mtext("Predictive of Species", side=2, line=2.7, cex=2)
abline(v=log(median(anova_stat_cor_asym2[idx_asym2,2], na.rm=T)))
abline(h=log(median(anova_stat_cor_asym2[idx_asym2,1], na.rm=T)))
abline(v=log(min(anova_stat_cor_asym2[top_list_asym2, 2])), lty=2)
abline(h=log(max(anova_stat_cor_asym2[top_list_asym2, 1])), lty=2)
dev.off()

png("sim_heat_sp.png", width=1200, height=1200, res=300, pointsize=5)
heatmap(ref_mat_sp, Rowv=NA, Colv=NA)
dev.off()

png("sim_heat_ct.png", width=1200, height=1200, res=300, pointsize=5)
heatmap(ref_mat_ct, Rowv=NA, Colv=NA)
dev.off()

png("sim_heat_comb.png", width=1200, height=1200, res=300, pointsize=5)
heatmap(exp_mat_all, Rowv=NA, Colv=NA)
dev.off()

png("sim_heat_cor.png", width=7, height=7, res=1200, units="in", pointsize=5)
heatmap(log(exp_mat_cor+1), Rowv=NA, Colv=NA)
dev.off()

png("sim_heat_asym.png", width=7, height=7, res=1200, units="in", pointsize=5)
heatmap(log(exp_mat_all_asym+1), Rowv=NA, Colv=NA, scale="none")
dev.off()

png("sim_heat_asym2.png", width=7, height=7, res=1200, units="in", pointsize=5)
heatmap(log(exp_mat_all_asym2+1), Rowv=NA, Colv=NA, scale="none")
dev.off()


png("sim_umap_sp.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_sp, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data (Species)")
dev.off()

png("sim_umap_ct.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_ct, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data (Cell Type)")
dev.off()

png("sim_umap_comb.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_comb, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data, All Genes")
dev.off()

png("sim_umap_cor.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_cor, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data, All Genes")
dev.off()

png("sim_umap_all_asym.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_cor_asym, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data, All Genes")
dev.off()

png("sim_umap_all_asym2.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_cor_asym2, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data, All Genes")
dev.off()


png("sim_umap_cor_fix.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_cor_fix, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data, Top 200 Genes")
dev.off()

png("sim_umap_cor_fix_ct.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_cor_fix, color_cells_by="ct_lab", group_cells_by="ct_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data, Top 200 Genes")
dev.off()

png("sim_umap_asym.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_asym_fix, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data, Top 200 Genes")
dev.off()

png("sim_umap_asym_ct.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_asym_fix, color_cells_by="ct_lab", group_cells_by="ct_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data, Top 200 Genes")
dev.off()

png("sim_umap_asym2.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_asym2_fix, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data, Top 200 Genes")
dev.off()

png("sim_umap_asym2_ct.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_asym2_fix, color_cells_by="ct_lab", group_cells_by="ct_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data, Top 200 Genes")
dev.off()


png("sim_umap_cor_batch.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_cor_batch, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data, Batch Corrected")
dev.off()

png("sim_pca_sp.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_sp, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F, reduction_method="PCA") +
	ggtitle("Simulated Data (Species)")
dev.off()

png("sim_pca_ct.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_ct, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F, reduction_method="PCA") +
	ggtitle("Simulated Data (Cell Type)")
dev.off()

png("sim_pca_comb.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_comb, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F, reduction_method="PCA") +
	ggtitle(paste("Simulated Data, p_mix=", p_mix, sep=""))
dev.off()

png("sim_pca_cor.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_cor, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F, reduction_method="PCA") +
	ggtitle(paste("Simulated Data, p_mix=", p_mix, sep=""))
dev.off()

png("sim_pca_cor_fix.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_cor_fix, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F, reduction_method="PCA") +
	ggtitle("Simulated Data, p_mix=", p_mix, sep=""))
dev.off()

png("sim_pca_cor_fix_ct.png", width=1200, height=1200, res=300, pointsize=5)
plot_cells(sim_dat_cor_fix, color_cells_by="ct_lab", group_cells_by="ct_lab", group_label_size = 3, label_cell_groups = F, reduction_method="PCA") +
	ggtitle(paste("Simulated Data, p_mix=", p_mix, sep=""))
dev.off()

png("sim_beta_vals.png", width=1200, height=1200, res=300, pointsize=5)
plot(density(beta_vals[which(rownames(anova_stat_cor) %in% top_list)]), xlim=c(0, 1), lwd=2, col="red", main="Cell Type/Species Bias")
lines(density(beta_vals), lwd=2, col="black")
legend(0.75, 6, legend=c("All Genes", "Top 200 Genes"), col=c("black", "red"), lwd=3)
dev.off()




# png("test.png")
# exp_mat_cor_df <- expand.grid(rowz=rownames(exp_mat_cor), colz=colnames(exp_mat_cor))
# exp_mat_cor_df$vals <- as.numeric(exp_mat_cor)
# ggplot(exp_mat_cor_df, aes(rowz, colz, fill=vals)) + geom_tile()
# dev.off()

svg("sim_heat_cor.svg", width=4, height=4, pointsize=5)
#no render

dev.off()

svg("sim_umap_sp.svg", width=4, height=4, pointsize=5)
plot_cells(sim_dat_sp, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data (Species)")
dev.off()

svg("sim_umap_ct.svg", width=4, height=4, pointsize=5)
plot_cells(sim_dat_ct, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data (Cell Type)")
dev.off()


svg("sim_umap_cor.svg", width=4, height=4, pointsize=5)
plot_cells(sim_dat_cor, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data, All Genes")
dev.off()

svg("sim_umap_cor_fix.svg", width=4, height=4, pointsize=5)
plot_cells(sim_dat_cor_fix, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data, Top 200 Genes")
dev.off()

svg("sim_umap_cor_fix_ct.svg", width=4, height=4, pointsize=5)
plot_cells(sim_dat_cor_fix, color_cells_by="ct_lab", group_cells_by="ct_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data, Top 200 Genes")
dev.off()

svg("sim_umap_cor_batch.svg", width=4, height=4, pointsize=5)
plot_cells(sim_dat_cor_batch, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F) +
	ggtitle("Simulated Data, Batch Corrected")
dev.off()

svg("sim_pca_sp.svg", width=4, height=4, pointsize=5)
plot_cells(sim_dat_sp, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F, reduction_method="PCA") +
	ggtitle("Simulated Data (Species)")
dev.off()

svg("sim_pca_ct.svg", width=4, height=4, pointsize=5)
plot_cells(sim_dat_ct, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F, reduction_method="PCA") +
	ggtitle("Simulated Data (Cell Type)")
dev.off()

svg("sim_pca_comb.svg", width=4, height=4, pointsize=5)
plot_cells(sim_dat_comb, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F, reduction_method="PCA") +
	ggtitle(paste("Simulated Data, p_mix=", p_mix, sep=""))
dev.off()

svg("sim_pca_cor.svg", width=4, height=4, pointsize=5)
plot_cells(sim_dat_cor, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F, reduction_method="PCA") +
	ggtitle(paste("Simulated Data, p_mix=", p_mix, sep=""))
dev.off()

svg("sim_pca_cor_fix.svg", width=4, height=4, pointsize=5)
plot_cells(sim_dat_cor_fix, color_cells_by="species_lab", group_cells_by="species_lab", group_label_size = 3, label_cell_groups = F, reduction_method="PCA") +
	ggtitle("Simulated Data, p_mix=", p_mix, sep=""))
dev.off()

svg("sim_pca_cor_fix_ct.svg", width=4, height=4, pointsize=5)
plot_cells(sim_dat_cor_fix, color_cells_by="ct_lab", group_cells_by="ct_lab", group_label_size = 3, label_cell_groups = F, reduction_method="PCA") +
	ggtitle(paste("Simulated Data, p_mix=", p_mix, sep=""))
dev.off()

svg("sim_beta_vals.svg", width=4, height=4, pointsize=5)
plot(density(beta_vals[which(rownames(anova_stat_cor) %in% top_list)]), xlim=c(0, 1), lwd=2, col="red", main="Cell Type/Species Bias")
lines(density(beta_vals), lwd=2, col="black")
legend(0.75, 6, legend=c("All Genes", "Top 200 Genes"), col=c("black", "red"), lwd=3)
dev.off()

svg("sim_anova_cor.svg", width=4, height=4,pointsize=5)
anova_main <- bquote("Simulated: ANOVA log(F-Statistic) by Gene, " ~ rho == .(an_cor_cor))
plot(log(anova_stat_cor[idx_cor,2]), log(anova_stat_cor[idx_cor,1]), main=anova_main, xlab="", ylab="", cex.axis=1.75, cex.main=2)
mtext("Predictive of Cell Type", side=1, line=3.2, cex=2)
mtext("Predictive of Species", side=2, line=2.7, cex=2)
abline(v=log(median(anova_stat_cor[idx_cor,2], na.rm=T)))
abline(h=log(median(anova_stat_cor[idx_cor,1], na.rm=T)))
abline(v=log(min(anova_stat_cor[top_list, 2])), lty=2)
abline(h=log(max(anova_stat_cor[top_list, 1])), lty=2)
dev.off()


save(list=ls(), file="sims.rdata")