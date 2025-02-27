#srun -p zhao --pty apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis:/home --home /home /lustre/fs4/zhao_lab/scratch/ulee/docka/samap_fig.sif R
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(extrafont)


mel_dat <- Read10X(data.dir="/home/sc_liftoff/dmel_testis_fb/outs/filtered_feature_bc_matrix")
write.csv(mel_dat, file="/home/sc_liftoff/bursting/dmel_testis_fb_raw.csv")
yak_dat <- Read10X(data.dir="/home/sc_liftoff/dyak_testis_fb/outs/filtered_feature_bc_matrix")
write.csv(yak_dat, file="/home/sc_liftoff/bursting/dyak_testis_fb_raw.csv")
ana_dat <- Read10X(data.dir="/home/sc_liftoff/dana_testis_fb/outs/filtered_feature_bc_matrix")
write.csv(ana_dat, file="/home/sc_liftoff/bursting/dana_testis_fb_raw.csv")

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

mel_seu <- CreateSeuratObject(counts=mel_dat_shared, project="melyakana_fb", min.cells=3, min.features=200)
yak_seu <- CreateSeuratObject(counts=yak_dat_shared, project="melyakana_fb", min.cells=3, min.features=200)
ana_seu <- CreateSeuratObject(counts=ana_dat_shared, project="melyakana_fb", min.cells=3, min.features=200)

# mel_seu <- subset(mel_seu, subset=nFeature_RNA > 200 & nFeature_RNA < 2500)
# yak_seu <- subset(yak_seu, subset=nFeature_RNA > 200 & nFeature_RNA < 2500)
# ana_seu <- subset(ana_seu, subset=nFeature_RNA > 200 & nFeature_RNA < 2500)

mel_seu <- NormalizeData(mel_seu)
yak_seu <- NormalizeData(yak_seu)
ana_seu <- NormalizeData(ana_seu)

mel_seu <- FindVariableFeatures(mel_seu, selection.method = "vst", nfeatures = 2000)
yak_seu <- FindVariableFeatures(yak_seu, selection.method = "vst", nfeatures = 2000)
ana_seu <- FindVariableFeatures(ana_seu, selection.method = "vst", nfeatures = 2000)

mel_seu <- ScaleData(mel_seu, features=genes_shared)
yak_seu <- ScaleData(yak_seu, features=genes_shared)
ana_seu <- ScaleData(ana_seu, features=genes_shared)

seu_list <- list(mel_seu, yak_seu, ana_seu)



features <- SelectIntegrationFeatures(object.list=seu_list)
seu_anchors <- FindIntegrationAnchors(object.list=seu_list, anchor.features=features)

seu_anchors_all <- FindIntegrationAnchors(object.list=seu_list, anchor.features=genes_shared)

seu_comb <- IntegrateData(seu_anchors)

seu_comb_all_shared <- IntegrateData(seu_anchors_all, features.to.integrate=genes_shared)

seu_comb <- FindVariableFeatures(object = seu_comb, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

DefaultAssay(seu_comb) <- "integrated"
seu_comb <- ScaleData(seu_comb, verbose=F)
seu_comb_all_shared <- ScaleData(seu_comb_all_shared, verbose=F)


ct_assign <- read.csv(file="sc_liftoff/ct_assignments.csv", row.names="X")
seu_comb <- AddMetaData(seu_comb, ct_assign$reclu_coarse, "reclu_coarse")
seu_comb <- AddMetaData(seu_comb, ct_assign$spec, "spec")
Idents(seu_comb) <- "reclu_coarse"

seu_comb_all_shared <- AddMetaData(seu_comb_all_shared, ct_assign$reclu_coarse, "reclu_coarse")
seu_comb_all_shared <- AddMetaData(seu_comb_all_shared, ct_assign$spec, "spec")
Idents(seu_comb_all_shared) <- "reclu_coarse"

seu_comb_byspec <- SplitObject(seu_comb, split.by="spec")

summary_mel <- table(Idents(seu_comb_byspec$mel))
summary_yak <- table(Idents(seu_comb_byspec$yak))
summary_ana <- table(Idents(seu_comb_byspec$ana))

seu_comb@active.ident <- factor(seu_comb@active.ident, levels=c("Somatic", "GSC/Early spermatogonia", "Late spermatogonia", "ananassae spermatocyte", "Early spermatocyte", "Late spermatocyte", "Early spermatid",  "Late spermatid"))

seu_comb_markers_allspec <- FindAllMarkers(seu_comb)
seu_comb_markers_mel <- FindAllMarkers(seu_comb_byspec$mel)
seu_comb_markers_yak <- FindAllMarkers(seu_comb_byspec$yak)
seu_comb_markers_ana <- FindAllMarkers(seu_comb_byspec$ana)

seu_comb_markers_allspec$cluster <- factor(seu_comb_markers_allspec$cluster, levels=c("Somatic", "GSC/Early spermatogonia", "Late spermatogonia", "ananassae spermatocyte", "Early spermatocyte", "Late spermatocyte", "Early spermatid",  "Late spermatid"))

marks_top <- seu_comb_markers_allspec %>% group_by(cluster) %>% slice_max(order_by=avg_log2FC, n=5)

marks_top_10 <- seu_comb_markers_allspec %>% group_by(cluster) %>% slice_max(order_by=avg_log2FC, n=50)

marks_top_all <- seu_comb_markers_allspec %>% group_by(cluster) %>% slice_max(order_by=avg_log2FC, n=999999)
temp <- unlist(lapply(1:8, function(x) sum(as.numeric(marks_top_all$cluster) == x)))
temp <- unlist(lapply(temp, (function(x) 1:x)), recursive=T)
marks_top_all <- data.frame(rank=temp, marks_top_all)

write.table(marks_top_10[1:50 + 50*7,], file="sc_liftoff/an_top_10.tsv", sep="\t", quote=F)
write.table(marks_top_all, file="sc_liftoff/all_top_all.tsv", sep="\t", quote=F, row.names=F)


marks_top_1 <- seu_comb_markers_allspec %>% group_by(cluster) %>% slice_max(order_by=avg_log2FC, n=1)

genes_converter <- genes_shared_fbgn_mel
names(genes_converter) <- genes_shared

marks_top_names_fbgn <- genes_converter[as.character(unlist(marks_top[,"gene"]))]

write.table(marks_top_names_fbgn, file="sc_liftoff/top_marks_seurat_afterct.txt", quote=F, row.names=F, col.names=F)

loadfonts()

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

col_wheel <- gg_color_hue(8)

ct_cols <- col_wheel[c(7, 3, 6, 2, 5, 1, 4, 8)]

png("sc_liftoff/vln_marks_afterct.png", width=2400*3.2, height=1200*8, res=300, pointsize=2)
p <- VlnPlot(seu_comb, features=as.character(unlist(marks_top[,"gene"])), ncol=5, pt.size=0.001, cols=ct_cols)
p <- p & theme(text=element_text(family="Roboto"))
p
dev.off()

png("sc_liftoff/vln_marks_afterct_nopts.png", width=2400*3.2, height=1200*8, res=300, pointsize=2)
p <- VlnPlot(seu_comb, features=as.character(unlist(marks_top[,"gene"])), ncol=5, pt.size=0, cols=ct_cols)
p <- p & theme(text=element_text(family="Roboto"))
p
dev.off()

png("sc_liftoff/vln_marks_afterct_nopts_tops.png", width=2400*8.2, height=1200*6, res=300, pointsize=2)
p <- VlnPlot(seu_comb, features=as.character(unlist(marks_top_1[,"gene"])), ncol=4, pt.size=0, cols=ct_cols)
p <- p & theme(text=element_text(family="Roboto", size=50), axis.text.x=element_text(size=40), plot.margin=margin(l=2, unit="cm"))
p <- p & labs(x="")
p
dev.off()


png("sc_liftoff/dots_afterct_1.png", width=8000, height=1200, res=300, pointsize=3)
p <- DotPlot(seu_comb, features=as.character(unlist(marks_top[1:20,"gene"])), cluster.idents=F)
p <- p & theme(text=element_text(family="Roboto"))
p
dev.off()

png("sc_liftoff/dots_afterct_2.png", width=8000, height=1200, res=300, pointsize=3)
p <- DotPlot(seu_comb, features=as.character(unlist(marks_top[21:40,"gene"])), cluster.idents=F)
p <- p & theme(text=element_text(family="Roboto"))
p
dev.off()


svg("sc_liftoff/vln_marks_afterct.svg", width=12.8, height=32, pointsize=2)
p <- VlnPlot(seu_comb, features=as.character(unlist(marks_top[,"gene"])), ncol=5, pt.size=0.001, cols=ct_cols)
p <- p & theme(text=element_text(family="Roboto"))
p
dev.off()

svg("sc_liftoff/vln_marks_afterct_nopts.svg", width=25.6, height=32, pointsize=2)
p <- VlnPlot(seu_comb, features=as.character(unlist(marks_top[,"gene"])), ncol=5, pt.size=0, cols=ct_cols)
p <- p & theme(text=element_text(family="Roboto"))
p
dev.off()

svg("sc_liftoff/vln_marks_afterct_nopts_tops.svg", width=65.6, height=24, pointsize=2)
p <- VlnPlot(seu_comb, features=as.character(unlist(marks_top_1[,"gene"])), ncol=4, pt.size=0, cols=ct_cols)
p <- p & theme(text=element_text(family="Roboto", size=50), axis.text.x=element_text(size=40), plot.margin=margin(l=2, unit="cm"))
p <- p & labs(x="")
p
dev.off()


svg("sc_liftoff/dots_afterct_1.svg", width=26.6666, height=4, pointsize=3)
p <- DotPlot(seu_comb, features=as.character(unlist(marks_top[1:20,"gene"])), cluster.idents=F)
p <- p & theme(text=element_text(family="Roboto"))
p
dev.off()

svg("sc_liftoff/dots_afterct_2.svg", width=26.6666, height=4, pointsize=3)
p <- DotPlot(seu_comb, features=as.character(unlist(marks_top[21:40,"gene"])), cluster.idents=F)
p <- p & theme(text=element_text(family="Roboto"))
p
dev.off()


dnt_dat <- read.csv("sc_liftoff/dnt_data.csv", header=T, colClasses=c("spec_code"="character"))

all_burst_size_mean <- read.csv("sc_liftoff/bursting/all_burst_size_mean.csv", header=T, row.names="X")
all_burst_freq_mean <- read.csv("sc_liftoff/bursting/all_burst_freq_mean.csv", header=T, row.names="X")

dnt_L2_size_mean <- read.csv("sc_liftoff/bursting/dnt_L2_size_mean.csv", header=T, row.names="X")
dnt_L1_size_mean <- read.csv("sc_liftoff/bursting/dnt_L1_size_mean.csv", header=T, row.names="X")
dnt_R1_size_mean <- read.csv("sc_liftoff/bursting/dnt_R1_size_mean.csv", header=T, row.names="X")
dnt_R2_size_mean <- read.csv("sc_liftoff/bursting/dnt_R2_size_mean.csv", header=T, row.names="X")
dnt_L2_freq_mean <- read.csv("sc_liftoff/bursting/dnt_L2_freq_mean.csv", header=T, row.names="X")
dnt_L1_freq_mean <- read.csv("sc_liftoff/bursting/dnt_L1_freq_mean.csv", header=T, row.names="X")
dnt_R1_freq_mean <- read.csv("sc_liftoff/bursting/dnt_R1_freq_mean.csv", header=T, row.names="X")
dnt_R2_freq_mean <- read.csv("sc_liftoff/bursting/dnt_R2_freq_mean.csv", header=T, row.names="X")

#fixing headers after renaming ananssae cell type after HCR! run only once, before txburst is re-run
reorder_an_cytes <- function(in_df){
	temp <- in_df
	temp <- data.frame(temp[,1:5], temp[,10], temp[,6:9], temp[,11])
	names(temp) <- c(names(in_df)[1:5], "ananassae_spermatocyte", names(in_df)[6:9], names(in_df)[11])
	return(temp)
}

all_burst_size_mean <- reorder_an_cytes(all_burst_size_mean)
all_burst_freq_mean <- reorder_an_cytes(all_burst_freq_mean)

dnt_L2_size_mean <- reorder_an_cytes(dnt_L2_size_mean)
dnt_L1_size_mean <- reorder_an_cytes(dnt_L1_size_mean)
dnt_R1_size_mean <- reorder_an_cytes(dnt_R1_size_mean)
dnt_R2_size_mean <- reorder_an_cytes(dnt_R2_size_mean)
dnt_L2_freq_mean <- reorder_an_cytes(dnt_L2_freq_mean)
dnt_L1_freq_mean <- reorder_an_cytes(dnt_L1_freq_mean)
dnt_R1_freq_mean <- reorder_an_cytes(dnt_R1_freq_mean)
dnt_R2_freq_mean <- reorder_an_cytes(dnt_R2_freq_mean)



#make a feature that is mean exprs level for l2, l1, r1, r2 and all (no-zero) genes
#vlnplot for those features by species and cell type

seu_comb_all_shared@active.ident <- factor(seu_comb_all_shared@active.ident, levels=c("Somatic", "GSC/Early spermatogonia", "Late spermatogonia", "ananassae spermatocyte", "Early spermatocyte", "Late spermatocyte", "Early spermatid",  "Late spermatid"))

gene_list_indat <- function(whi_table, whi_spec, whi_ct){
	return(whi_table[which(whi_table$species == whi_spec & whi_table[, whi_ct] > 0), "Symbol"])
}

calc_mean_by_cell <- function(in_seu, cats_list){
	exp_mtx <- in_seu@assays$RNA@counts
	in_specs <- in_seu@meta.data$spec
	in_clusts <- in_seu@meta.data$reclu_coarse
	
	in_clusts <- unlist(lapply(in_clusts, function(x) gsub(" ", "_", x)))
	in_clusts <- unlist(lapply(in_clusts, function(x) gsub("/", "_", x)))

	
	out_meta_dat <- matrix(0, nrow=length(cats_list), ncol=ncol(exp_mtx))
	rownames(out_meta_dat) <- names(cats_list)
	colnames(out_meta_dat) <- colnames(exp_mtx)

	for(i in 1:ncol(exp_mtx)){	
		props_list <- lapply(cats_list, function(x){
				print(c(i, in_specs[i], in_clusts[i]))
				gen_list <- gene_list_indat(x, in_specs[i], in_clusts[i])
				gen_list <- intersect(gen_list, rownames(exp_mtx))

				return(mean(exp_mtx[gen_list, i], na.rm=T))
			})
		
		names(props_list) <- names(cats_list)
		
		out_meta_dat[,i] <- unlist(props_list)
	}
	
	return(out_meta_dat)
}

cat_size_list <- list(all_burst_size_mean, dnt_L2_size_mean, dnt_L1_size_mean, dnt_R1_size_mean, dnt_R2_size_mean)

names(cat_size_list) <- c("size_all", "size_L2", "size_L1", "size_R1", "size_R2")

meta_dat_size <- calc_mean_by_cell(seu_comb_all_shared, cat_size_list)

seu_comb_all_shared[["size_all"]] <- meta_dat_size["size_all",]
seu_comb_all_shared[["size_L2"]] <- meta_dat_size["size_L2",]
seu_comb_all_shared[["size_L1"]] <- meta_dat_size["size_L1",]
seu_comb_all_shared[["size_R1"]] <- meta_dat_size["size_R1",]
seu_comb_all_shared[["size_R2"]] <- meta_dat_size["size_R2",]


seu_comb_all_shared_copy <- seu_comb_all_shared
seu_comb_all_shared_copy <- RenameCells(seu_comb_all_shared_copy, add.cell.id="copy")

seu_comb_all_shared_copy[["size_L2"]] <- meta_dat_size["size_all",]
seu_comb_all_shared_copy[["size_L1"]] <- meta_dat_size["size_all",]
seu_comb_all_shared_copy[["size_R1"]] <- meta_dat_size["size_all",]
seu_comb_all_shared_copy[["size_R2"]] <- meta_dat_size["size_all",]

seu_comb_all_copy_merge <- merge(seu_comb_all_shared, seu_comb_all_shared_copy)

seu_comb_all_copy_merge[["dist"]] <- c(rep("group", 15000), rep("all", 15000))

seu_comb_all_copy_merge@active.ident <- factor(seu_comb_all_copy_merge@active.ident, levels=c("Somatic", "GSC/Early spermatogonia", "Late spermatogonia", "ananassae spermatocyte", "Early spermatocyte", "Late spermatocyte", "Early spermatid",  "Late spermatid"))


seu_comb_all_copy_merge_spec <- SplitObject(seu_comb_all_copy_merge, split.by="spec")


svg("sc_liftoff/vln_neis.svg", width=12.8*2, height=12.8*0.625*2, pointsize=2)
p1 <- VlnPlot(seu_comb_all_copy_merge_spec[["mel"]], features=c("size_L2", "size_L1", "size_R1", "size_R2"), split.by="dist", split.plot=TRUE, same.y.lims=TRUE, y.max=25, ncol=4, pt.size=0, cols=ct_cols)
p2 <- VlnPlot(seu_comb_all_copy_merge_spec[["yak"]], features=c("size_L2", "size_L1", "size_R1", "size_R2"), split.by="dist", split.plot=TRUE, same.y.lims=TRUE, y.max=45, ncol=4, pt.size=0, cols=ct_cols)
p3 <- VlnPlot(seu_comb_all_copy_merge_spec[["ana"]], features=c("size_L2", "size_L1", "size_R1", "size_R2"), split.by="dist", split.plot=TRUE, same.y.lims=TRUE, y.max=10, ncol=4, pt.size=0, cols=ct_cols)
p1/p2/p3
dev.off()

png("sc_liftoff/vln_neis.png", width=2400*3.2, height=1200*8/2, res=300, pointsize=2)
p1/p2/p3

dev.off()

pv_counts_mtx_L2 <- matrix(0, nrow=3, ncol=8)
rownames(pv_counts_mtx_L2) <- names(seu_comb_all_copy_merge_spec)
colnames(pv_counts_mtx_L2) <- levels(seu_comb_all_copy_merge_spec[[1]]@active.ident)

pv_counts_mtx_L1 <- pv_counts_mtx_L2
pv_counts_mtx_R1 <- pv_counts_mtx_L2
pv_counts_mtx_R2 <- pv_counts_mtx_L2

rownames(pv_counts_mtx_L2) <- unlist(lapply(rownames(pv_counts_mtx_L2), function(x) paste(x, "_L2", sep="")))
rownames(pv_counts_mtx_L1) <- unlist(lapply(rownames(pv_counts_mtx_L1), function(x) paste(x, "_L1", sep="")))
rownames(pv_counts_mtx_R1) <- unlist(lapply(rownames(pv_counts_mtx_R1), function(x) paste(x, "_R1", sep="")))
rownames(pv_counts_mtx_R2) <- unlist(lapply(rownames(pv_counts_mtx_R2), function(x) paste(x, "_R2", sep="")))


for(i in 1:3){
	this_seu <- seu_comb_all_copy_merge_spec[[i]]
	
	for(j in 1:8){
	
		this_idx_group <- intersect(which(this_seu@active.ident == levels(this_seu@active.ident)[j]), which(this_seu[["dist"]] == "group"))
		this_idx_all <- intersect(which(this_seu@active.ident == levels(this_seu@active.ident)[j]), which(this_seu[["dist"]] == "all"))
		
		
		pv_counts_mtx_L2[i, j] <- t.test(this_seu[["size_L2"]][this_idx_group,], this_seu[["size_L2"]][this_idx_all,])$p.value
		pv_counts_mtx_L1[i, j] <- t.test(this_seu[["size_L1"]][this_idx_group,], this_seu[["size_L1"]][this_idx_all,])$p.value
		pv_counts_mtx_R1[i, j] <- t.test(this_seu[["size_R1"]][this_idx_group,], this_seu[["size_R1"]][this_idx_all,])$p.value
		pv_counts_mtx_R2[i, j] <- t.test(this_seu[["size_R2"]][this_idx_group,], this_seu[["size_R2"]][this_idx_all,])$p.value

	}
}

pv_counts_mtx_all <- do.call("rbind", list(pv_counts_mtx_L2, pv_counts_mtx_L1, pv_counts_mtx_R1, pv_counts_mtx_R2))

benj_hoch <- function(in_mtx, in_fdr){
	this_mtx <- in_mtx
	this_fdr <- in_fdr
	temp <- as.numeric(this_mtx)
	temp1 <- temp[order(temp)]
	temp2 <- rbind(temp1, this_fdr*(1:length(unlist(in_mtx))/length(unlist(in_mtx))))
	temp3 <- temp2[1,] - temp2[2,]
	if(sum(temp3 <= 0) == 0) return(rep(0, 3))
	
	temp4 <- which(temp3 <= 0)
	
	temp5 <- rep(0, 3)
	if(length(temp4) > 0){
		temp5 <- lapply(temp1[temp4], function(x) which(this_mtx == x, arr.ind=TRUE))
		temp5 <- do.call("rbind", temp5)
		temp5 <- cbind(temp5, temp1[temp4])
		temp6 <- data.frame(temp5)
		rownames(temp6) <- unlist(lapply(1:nrow(temp5), function(x) paste("result_", x, sep="")))
		temp6[,1] <- rownames(in_mtx)[temp5[,1]]
		temp6[,2] <- colnames(in_mtx)[temp5[,2]]
		
		colnames(temp6)[3] <- "pval"
	}

	return(temp6)
}

bh_counts_mtx_all <- benj_hoch(pv_counts_mtx_all, 0.05)

write.csv(pv_counts_mtx_all, file="sc_liftoff/pv_counts_mtx_all.csv")

write.csv(pv_counts_mtx_all < 0.05/96, file="sc_liftoff/pv_counts_mtx_all_bool.csv")

save(list=ls(), file="sc_liftoff/wrkspce_run_seurat_after_ctassign.rdata")
