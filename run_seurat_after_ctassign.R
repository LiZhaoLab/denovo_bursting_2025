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

seu_comb <- IntegrateData(seu_anchors)

seu_comb <- FindVariableFeatures(object = seu_comb, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

DefaultAssay(seu_comb) <- "integrated"
seu_comb <- ScaleData(seu_comb, verbose=F)

ct_assign <- read.csv(file="sc_liftoff/ct_assignments.csv", row.names="X")
seu_comb <- AddMetaData(seu_comb, ct_assign$reclu_coarse, "reclu_coarse")
seu_comb <- AddMetaData(seu_comb, ct_assign$spec, "spec")
Idents(seu_comb) <- "reclu_coarse"

seu_comb_byspec <- SplitObject(seu_comb, split.by="spec")

summary_mel <- table(Idents(seu_comb_byspec$mel))
summary_yak <- table(Idents(seu_comb_byspec$yak))
summary_ana <- table(Idents(seu_comb_byspec$ana))

seu_comb@active.ident <- factor(seu_comb@active.ident, levels=c("Somatic", "GSC/Early spermatogonia", "Late spermatogonia", "Early spermatocyte", "Late spermatocyte", "Early spermatid",  "Late spermatid", "ananassae spermatid"))

seu_comb_markers_allspec <- FindAllMarkers(seu_comb)
seu_comb_markers_mel <- FindAllMarkers(seu_comb_byspec$mel)
seu_comb_markers_yak <- FindAllMarkers(seu_comb_byspec$yak)
seu_comb_markers_ana <- FindAllMarkers(seu_comb_byspec$ana)

seu_comb_markers_allspec@cluster <- factor(seu_comb_markers_allspec@cluster, levels=c("Somatic", "GSC/Early spermatogonia", "Late spermatogonia", "Early spermatocyte", "Late spermatocyte", "Early spermatid",  "Late spermatid", "ananassae spermatid"))

marks_top <- seu_comb_markers_allspec %>% group_by(cluster) %>% slice_max(order_by=avg_log2FC, n=5)

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

svg("sc_liftoff/vln_marks_afterct_nopts.svg", width=12.8, height=32, pointsize=2)
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



save(list=ls(), file="sc_liftoff/wrkspce_run_seurat_after_ctassign.rdata")