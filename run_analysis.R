#srun -p zhao --pty apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis:/home --home /home /lustre/fs4/zhao_lab/scratch/ulee/docka/garnett.sif R
library(stringr)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(ggsignif)

#bursting analysis
load("sc_liftoff/bursting/allpickslist.rdata")
dnt_dat <- read.csv("sc_liftoff/dnt_data.csv", header=T, colClasses=c("spec_code"="character"))

plot_new <- FALSE

extract_geneid <- function(a) return(str_match(a, "gene_id\\s*(.*?)\\s*;")[1,2])
extract_genesymbol <- function(a) return(str_match(a, "gene_symbol\\s*(.*?)\\s*;")[1,2])

orthos_yak <- read.csv("/home/sc_liftoff/genomes_fb/dmel_orthos_yak.tsv", skip=5, header=F, sep="\t")
orthos_ana <- read.csv("/home/sc_liftoff/genomes_fb/dmel_orthos_ana.tsv", skip=5, header=F, sep="\t")

mel_gtf <- read.csv("/home/supacoil/dmel_r_6_50_geneonly_sorted.gtf", header=F, sep="\t")
mel_gtf_older <- read.csv("/home/supacoil/dmel_r_6_15_geneonly_sorted.gtf", header=F, sep="\t")

mel_symb_all <- unlist(lapply(mel_gtf[,9], function(x) extract_genesymbol(x)))
mel_geneid_all <- unlist(lapply(mel_gtf[,9], function(x) extract_geneid(x)))

mel_symb_all_old <- unlist(lapply(mel_gtf_older[,9], function(x) extract_genesymbol(x)))
mel_geneid_all_old <- unlist(lapply(mel_gtf_older[,9], function(x) extract_geneid(x)))

idx_ncrna <- union(which(mel_gtf[,3] %in% c("miRNA", "ncRNA", "pre_miRNA", "rRNA", "snRNA", "snoRNA", "tRNA", "pseudogene")), which(grepl("RNA:", mel_symb_all)))
idx_ncrna <- union(idx_ncrna, which(grepl("CR", mel_symb_all)))
idx_ncrna <- union(idx_ncrna, which(grepl("mir", mel_symb_all)))

idx_ncrna <- sort(idx_ncrna, decreasing=FALSE)

mel_symb_ncrna <- mel_symb_all[idx_ncrna]
mel_symb_genes <- setdiff(mel_symb_all, mel_symb_ncrna)

mel_geneid_ncrna <- mel_geneid_all[idx_ncrna]
mel_geneid_genes <- setdiff(mel_geneid_all, mel_geneid_ncrna)

mel_onlygenes <- setdiff(mel_symb_genes, mel_symb_ncrna)

symb_id_conv <- mel_geneid_all
names(symb_id_conv) <- mel_symb_all

symb_id_conv_rev <- names(symb_id_conv)
names(symb_id_conv_rev) <- symb_id_conv

symb_id_conv_old <- mel_geneid_all_old
names(symb_id_conv_old) <- mel_symb_all_old


#check for presence/absence in cell types and species
ct_labs <- c("GSC_Early_spermatogonia", "Late_spermatogonia", "Early_spermatocyte", "Late_spermatocyte", "Early_spermatid", "Late_spermatid", "ananassae_spermatid", "Somatic")
spec_labs <- c("mel", "yak", "ana")
spec_ct_labs <- lapply(spec_labs, function(x){
		unlist(lapply(ct_labs, function(y){
			return(paste(x, "_", y, "_ML", sep=""))
		}))
	})
	
spec_ct_labs <- do.call("rbind", spec_ct_labs)

#if given a gene symbol, identify yak and ana orthos, return two matrices, one for with QC one without
get_stat_mtx <- function(gen_symb, converter, species_labs, type_labs, ort_yak_mtx, ort_ana_mtx, picks_list, freq_size="burst_freq", qc=FALSE){
	null_res <- data.frame(matrix(0, nrow=length(species_labs), ncol=length(type_labs)))
	null_res <- cbind(species_labs, null_res)
	null_res <- cbind(rep(gen_symb, length(species_labs)), null_res)
	colnames(null_res) <- c("symbol", "species", type_labs)
	
	#lookup gene symbol in mel
	mel_fbgn <- NA
	
	if(!freq_size %in% c("burst_freq", "burst_size")){
		print("Invalid freq_size, try burst_freq or burst_size")
		return(null_res)
	}
	
	if(gen_symb %in% names(converter)) mel_fbgn <- converter[gen_symb]
	else{
		print("Missing symbol info")
		return(null_res)
	}
		
	orts_yak <- c()
	orts_yak_symb <- c()
	orts_ana <- c()
	orts_ana_symb <- c()
	
	if(mel_fbgn %in% ort_yak_mtx[,1]){
		orts_yak <- ort_yak_mtx[ort_yak_mtx[,1] == mel_fbgn, 6]
		orts_yak_symb <- unlist(lapply(ort_yak_mtx[ort_yak_mtx[,1] == mel_fbgn, 7], function(x) return(substr(x, 1, nchar(x)))))
	}
	
	if(mel_fbgn %in% ort_ana_mtx[,1]){
		orts_ana <- ort_ana_mtx[ort_ana_mtx[,1] == mel_fbgn, 6]
		orts_ana_symb <- unlist(lapply(ort_ana_mtx[ort_ana_mtx[,1] == mel_fbgn, 7], function(x) return(substr(x, 1, nchar(x)))))
	}
	
	print(paste("Number of yakuba orthologs found: ", length(orts_yak), sep=""))
	print(paste("Number of ananassae orthologs found: ", length(orts_ana), sep=""))
	
	if(sum(length(orts_yak) + length(orts_ana) == 0)){
		print("No orthologs found")
		return(null_res)
	}
	
	result <- as.data.frame(matrix(0, nrow=length(c(orts_yak_symb, orts_ana_symb)) + 1, ncol=length(type_labs)))
	
	iterator_symbs <- c(gen_symb, orts_yak_symb, orts_ana_symb)
	
	colnames(result) <- type_labs
	
	iterator <- c(mel_fbgn, orts_yak, orts_ana)
	iterator_specs <- c("mel", rep("yak", length(orts_yak)), rep("ana", length(orts_ana)))
	
	for(i in 1:length(iterator)){
		result[i, ] <- unlist(lapply(type_labs, function(x){
			this_list_lab <- paste(iterator_specs[i], "_", x, "_ML", sep="")
			this_pick <- picks_list[[this_list_lab]]
			
			if(!iterator[i] %in% rownames(this_pick)){
				print("Ortholog not detected")
				return(null_res)
			}
			else{
				this_return <- this_pick[which(rownames(this_pick) == iterator[i]), freq_size]
				this_qc <- !this_pick[which(rownames(this_pick) == iterator[i]), "keep"]
				if(qc && this_qc){
					this_return <- 0
				}
				
				return(this_return)
			}
			
		}))
	}
	
	result <- cbind(iterator_specs, result)
	result <- cbind(iterator_symbs, result)
	colnames(result)[1:2] <- c("symbol", "species")
	
	return(result)
}

test1 <- get_stat_mtx("CG9570", symb_id_conv, spec_labs, ct_labs, orthos_yak, orthos_ana, all_picks_list, "burst_freq", TRUE)
test2 <- get_stat_mtx("Adh", symb_id_conv, spec_labs, ct_labs, orthos_yak, orthos_ana, all_picks_list, "burst_freq", TRUE)
test3 <- get_stat_mtx("GstD4", symb_id_conv, spec_labs, ct_labs, orthos_yak, orthos_ana, all_picks_list, "burst_size", TRUE)


#characterize burst freq and burst size for all genes

all_freq_list <- lapply(names(symb_id_conv), function(x) get_stat_mtx(x, symb_id_conv, spec_labs, ct_labs, orthos_yak, orthos_ana, all_picks_list, "burst_freq", TRUE))
all_size_list <- lapply(names(symb_id_conv), function(x) get_stat_mtx(x, symb_id_conv, spec_labs, ct_labs, orthos_yak, orthos_ana, all_picks_list, "burst_size", TRUE))

names(all_freq_list) <- names(symb_id_conv)
names(all_size_list) <- names(symb_id_conv)

all_freq_list <- all_freq_list[!grepl("ana", names(all_freq_list))]
all_size_list <- all_size_list[!grepl("ana", names(all_size_list))]

colMeans_nozero <- function(in_mat){
	result <- apply(in_mat, 2, function(x){
		n_zero <- sum(x == 0)
		
		if(n_zero == length(x)) return(0)
		
		sumz <- sum(x)
		return(sumz/(length(x) - n_zero))
	})
	
	return(result)
}

collapse_stat_mtx <- function(stat_mtx, by.mean=TRUE){
	if(!by.mean) print("Using first row by species")
	
	spec_all <- c("mel", "yak", "ana")
	spec_levels <- levels(as.factor(stat_mtx[,"species"]))
	
	if(length(spec_levels) < length(spec_all)) return(stat_mtx[,-1])
	if(nrow(stat_mtx) == length(spec_all)) return(stat_mtx[,-1])
	
	stat_mtx <- as.data.frame(stat_mtx)

	result <- lapply(spec_all, function(x){
		if(sum(stat_mtx[,"species"] == x) == 1) return(stat_mtx[which(stat_mtx[,"species"] == x),-1:-2])
		
		submtx <- stat_mtx[which(stat_mtx[,"species"] == x),-1:-2]
		
		# submtx <- matrix(as.numeric(submtx), nrow=nrow(submtx))
		if(by.mean) return(colMeans_nozero(submtx))
		else return(submtx[1,])
	})
		
	result <- cbind(spec_all, as.data.frame(do.call("rbind", result)))
	colnames(result) <- colnames(stat_mtx)[-1]
	
	return(result)
}

test1_1 <- collapse_stat_mtx(test1)
test2_1 <- collapse_stat_mtx(test2)
test3_1 <- collapse_stat_mtx(test3, FALSE)

summary_all_freq_mean <- lapply(all_freq_list, function(x){collapse_stat_mtx(x, TRUE)})
summary_all_freq_first <- lapply(all_freq_list, function(x) collapse_stat_mtx(x, FALSE))

summary_all_size_mean <- lapply(all_size_list, function(x) collapse_stat_mtx(x, TRUE))
summary_all_size_first <- lapply(all_size_list, function(x) collapse_stat_mtx(x, FALSE))

save(list=ls(), file="sc_liftoff/bursting/wrkspce.rdata")

#collate as matrix and print
collate_as_mtx <- function(list_o_statmtx, var_lab){
	name_labs <- names(list_o_statmtx)
	
	list_withlabs <- list()
	
	for(i in 1:length(name_labs)){
		this_statmtx <- list_o_statmtx[[i]]
		this_res <- cbind(rep(name_labs[i], nrow(this_statmtx)), rep(var_lab, nrow(this_statmtx)))
		this_res <- cbind(this_res, this_statmtx)
		colnames(this_res) <- c("Symbol", "Property", colnames(this_statmtx))
		list_withlabs[[i]] <- this_res
	}
	
	result <- do.call("rbind", list_withlabs)
	return(result)
}

collate_two <- function(mtx1, mtx2){
	result <- rbind(mtx1, mtx2)
	
	result[seq(1, 2*nrow(mtx1) - 1, by=2),] <- mtx1
	result[seq(2, 2*nrow(mtx1), by=2),] <- mtx2
	
	return(result)
}

search_annots <- function(in_list, annots){
	result <- unlist(lapply(in_list, function(x){
		sear <- grepl(x, annots)
		if(sum(sear) >= 1) return(annots[sear][1])
		else return(x)
	}))

	return(result)
}

dnt_dat_annots <- search_annots(dnt_dat[,"symbol"], names(summary_all_freq_mean))

dnt_freq_mean <- summary_all_freq_mean[dnt_dat_annots]
dnt_size_mean <- summary_all_size_mean[dnt_dat_annots]

dnt_idx <- intersect(which(!is.na(names(dnt_freq_mean))), which(!is.na(names(dnt_size_mean))))

dnt_freq_mean <- dnt_freq_mean[dnt_idx]
dnt_size_mean <- dnt_size_mean[dnt_idx]

dnt_freq_table <- collate_as_mtx(dnt_freq_mean, "freq")

dnt_size_table <- collate_as_mtx(dnt_size_mean, "size")

dnt_burst_comb <- collate_two(dnt_freq_table, dnt_size_table)

write.csv(dnt_burst_comb, file="sc_liftoff/bursting/dnt_bursting_comb.csv")

gen_burst_table <- function(gene_list, freq_list, size_list, out_file){
	freq_mean <- freq_list[gene_list]
	size_mean <- size_list[gene_list]

	this_idx <- intersect(which(!is.na(names(freq_mean))), which(!is.na(names(size_mean))))

	freq_mean <- freq_mean[this_idx]
	size_mean <- size_mean[this_idx]

	freq_table <- collate_as_mtx(freq_mean, "freq")
	size_table <- collate_as_mtx(size_mean, "size")

	burst_comb <- collate_two(freq_table, size_table)

	write.csv(burst_comb, file=out_file)

	return(list(freq_table, size_table))
}

dnt_L2 <- gen_burst_table(symb_id_conv_rev[dnt_dat[,"nei_L2_fbgn_mel"]], summary_all_freq_mean, summary_all_size_mean, "sc_liftoff/bursting/dnt_L2_bursting_comb.csv")
dnt_L2_freq_mean <- dnt_L2[[1]]
dnt_L2_size_mean <- dnt_L2[[2]]

dnt_L1 <- gen_burst_table(symb_id_conv_rev[dnt_dat[,"nei_L1_fbgn_mel"]], summary_all_freq_mean, summary_all_size_mean, "sc_liftoff/bursting/dnt_L1_bursting_comb.csv")
dnt_L1_freq_mean <- dnt_L1[[1]]
dnt_L1_size_mean <- dnt_L1[[2]]

dnt_R1 <- gen_burst_table(symb_id_conv_rev[dnt_dat[,"nei_R1_fbgn_mel"]], summary_all_freq_mean, summary_all_size_mean, "sc_liftoff/bursting/dnt_R1_bursting_comb.csv")
dnt_R1_freq_mean <- dnt_R1[[1]]
dnt_R1_size_mean <- dnt_R1[[2]]

dnt_R2 <- gen_burst_table(symb_id_conv_rev[dnt_dat[,"nei_R2_fbgn_mel"]], summary_all_freq_mean, summary_all_size_mean, "sc_liftoff/bursting/dnt_R2_bursting_comb.csv")
dnt_R2_freq_mean <- dnt_R2[[1]]
dnt_R2_size_mean <- dnt_R2[[2]]

all_burst <- gen_burst_table(names(summary_all_freq_mean), summary_all_freq_mean, summary_all_size_mean, "sc_liftoff/bursting/all_bursting_comb.csv")
all_burst_freq_mean <- all_burst[[1]]
all_burst_size_mean <- all_burst[[2]]

#perform analysis
#plots of data by cell type, each point is species
#with or without 0s
rm_zero <- function(in_vec){
	idx_nonzero <- which(!in_vec == 0)
	return(in_vec[idx_nonzero])
}

format_df <- function(in_mat1, in_mat2, lab1, lab2){
	
	summz <- function(x){
		prop <- x[1,"Property"]
		ct_labby <- colnames(x)[4:11]
		
		idx_mel <- which(x[,"species"] == "mel")
		idx_yak <- which(x[,"species"] == "yak")
		idx_ana <- which(x[,"species"] == "ana")
	
		ct_dat_list_mel <- apply(x[idx_mel,4:11], 2, function(y) return(rm_zero(y)))
		ct_dat_list_yak <- apply(x[idx_yak,4:11], 2, function(y) return(rm_zero(y)))
		ct_dat_list_ana <- apply(x[idx_ana,4:11], 2, function(y) return(rm_zero(y)))
		
		names(ct_dat_list_mel) <- ct_labby
		names(ct_dat_list_yak) <- ct_labby
		names(ct_dat_list_ana) <- ct_labby
		
		result <- data.frame(spec="", cell_type="", parameter=0)
		
		for(z in ct_labby){
			this_vec_mel <- ct_dat_list_mel[[z]]
			this_vec_yak <- ct_dat_list_yak[[z]]
			this_vec_ana <- ct_dat_list_ana[[z]]
			
			result <- rbind(result, data.frame(spec=rep("mel", length(this_vec_mel)), cell_type=rep(z, length(this_vec_mel)), parameter=unlist(this_vec_mel)))
			result <- rbind(result, data.frame(spec=rep("yak", length(this_vec_yak)), cell_type=rep(z, length(this_vec_yak)), parameter=unlist(this_vec_yak)))
			result <- rbind(result, data.frame(spec=rep("ana", length(this_vec_ana)), cell_type=rep(z, length(this_vec_ana)), parameter=unlist(this_vec_ana)))
		}
		
		return(result[-1,])
	}
	
	big_result1 <- summz(in_mat1)
	big_result1 <- cbind(status=rep(lab1, nrow(big_result1)), big_result1)
	big_result2 <- summz(in_mat2)
	big_result2 <- cbind(status=rep(lab2, nrow(big_result2)), big_result2)
	
	return(rbind(big_result1, big_result2))
}

#plot for all genes vs L2, L1, R1, R2... for burst size and freq

plot_ind <- function(in_mat, exp_mat, in_mat_lab, exp_mat_lab, whi_spec, y_low, y_high, ylabel, in_title, fix_ana=FALSE){
	this_dat <- format_df(in_mat, exp_mat, in_mat_lab, exp_mat_lab)	
	this_dat <- this_dat[which(this_dat[,"spec"] == whi_spec),]
	
	ct <- this_dat$cell_type
	ct[which(ct == "GSC_Early_spermatogonia")] <- "GSC/Early Spermatogonia"
	ct[which(ct == "Late_spermatogonia")] <- "Late Spermatogonia"
	ct[which(ct == "Early_spermatocyte")] <- "Early Spermatocyte"
	ct[which(ct == "Late_spermatocyte")] <- "Late Spermatocyte"
	ct[which(ct == "Early_spermatid")] <- "Early Spermatid"
	ct[which(ct == "Late_spermatid")] <- "Late Spermatid"
	ct[which(ct == "ananassae_spermatid")] <- "ananassae Spermatocyte"
	
	stat <- this_dat$status
	stat[which(stat == "All_Genes")] <- "All Genes"
	stat[which(stat == "Neighbor_L2")] <- "Neighbor L2"
	stat[which(stat == "Neighbor_L1")] <- "Neighbor L1"
	stat[which(stat == "Neighbor_R1")] <- "Neighbor R1"
	stat[which(stat == "Neighbor_R2")] <- "Neighbor R2"
	stat[which(stat == "X_linked")] <- "X Linked"
	stat[which(stat == "Autosomal")] <- "Autosomal"
	
	this_dat$cell_type <- factor(ct, levels=c("Somatic", "GSC/Early Spermatogonia", "Late Spermatogonia", "ananassae Spermatocyte", "Early Spermatocyte", "Late Spermatocyte", "Early Spermatid",  "Late Spermatid"))
	
	this_dat$status <- factor(stat, levels=c("Neighbor L2", "Neighbor L1", "Neighbor R1", "Neighbor R2", "All Genes", "X Linked", "Autosomal"))
	
	if(fix_ana){
		if(whi_spec %in% c("mel", "yak")){
			this_dat <- subset(this_dat, cell_type %in% c("Somatic", "GSC/Early Spermatogonia", "Late Spermatogonia", "Early Spermatocyte", "Late Spermatocyte", "Early Spermatid",  "Late Spermatid"))
		}

		if(whi_spec == "ana"){
			this_dat <- subset(this_dat, cell_type %in% c("Somatic", "GSC/Early Spermatogonia", "Late Spermatogonia", "ananassae Spermatocyte", "Late Spermatocyte", "Early Spermatid",  "Late Spermatid"))
		}
		# ct <- subset(ct, 
	
	}
	
	ggplot(this_dat, aes(x=cell_type, y=parameter, fill=status), outlier.shape = NA) + geom_boxplot() + coord_cartesian(ylim=c(y_low, y_high)) + ylab(ylabel) + xlab("Cell Type") + labs(fill="Gene Set") + ggtitle(in_title)
}


plot_ind1 <- function(in_mat, exp_mat, in_mat_lab, exp_mat_lab, whi_spec, y_low, y_high, ylabel, in_title, fix_ana=FALSE){
	this_dat <- format_df(in_mat, exp_mat, in_mat_lab, exp_mat_lab)	
	this_dat <- this_dat[which(this_dat[,"spec"] == whi_spec),]
	
	ct <- this_dat$cell_type
	ct[which(ct == "GSC_Early_spermatogonia")] <- "GSC/Early Spermatogonia"
	ct[which(ct == "Late_spermatogonia")] <- "Late Spermatogonia"
	ct[which(ct == "Early_spermatocyte")] <- "Early Spermatocyte"
	ct[which(ct == "Late_spermatocyte")] <- "Late Spermatocyte"
	ct[which(ct == "Early_spermatid")] <- "Early Spermatid"
	ct[which(ct == "Late_spermatid")] <- "Late Spermatid"
	ct[which(ct == "ananassae_spermatid")] <- "ananassae Spermatocyte"
	
	stat <- this_dat$status
	stat[which(stat == "All_Genes")] <- "All Genes"
	stat[which(stat == "Neighbor_L2")] <- "Neighbor L2"
	stat[which(stat == "Neighbor_L1")] <- "Neighbor L1"
	stat[which(stat == "Neighbor_R1")] <- "Neighbor R1"
	stat[which(stat == "Neighbor_R2")] <- "Neighbor R2"
	stat[which(stat == "X_linked")] <- "X Linked"
	stat[which(stat == "Autosomal")] <- "Autosomal"
	
	this_dat$cell_type <- factor(ct, levels=c("Somatic", "GSC/Early Spermatogonia", "Late Spermatogonia", "ananassae Spermatocyte", "Early Spermatocyte", "Late Spermatocyte", "Early Spermatid",  "Late Spermatid"))
	
	this_dat$status <- factor(stat, levels=c("Neighbor L2", "Neighbor L1", "Neighbor R1", "Neighbor R2", "All Genes", "X Linked", "Autosomal"))
	
	if(fix_ana){
		if(whi_spec %in% c("mel", "yak")){
			this_dat <- subset(this_dat, cell_type %in% c("Somatic", "GSC/Early Spermatogonia", "Late Spermatogonia", "Early Spermatocyte", "Late Spermatocyte", "Early Spermatid", "Late Spermatid"))
		}

		if(whi_spec == "ana"){
			this_dat <- subset(this_dat, cell_type %in% c("Somatic", "GSC/Early Spermatogonia", "Late Spermatogonia", "ananassae Spermatocyte", "Late Spermatocyte", "Early Spermatid", "Late Spermatid"))
			this_dat$cell_type <- factor(this_dat$cell_type, levels=c("Somatic", "GSC/Early Spermatogonia", "Late Spermatogonia", "ananassae Spermatocyte", "Late Spermatocyte", "Early Spermatid", "Late Spermatid"))

			this_labs <- factor(c("Somatic", "GSC/Early Spermatogonia", "Late Spermatogonia", "ananassae Spermatocyte", "Late Spermatocyte", "Early Spermatid", "Late Spermatid"))
			
			return(ggplot(this_dat, aes(x=cell_type, y=parameter, fill=status), outlier.shape = NA) + geom_boxplot() + scale_x_discrete("cell_type", breaks=this_labs, drop=FALSE) + coord_cartesian(ylim=c(y_low, y_high)) + ylab(ylabel) + xlab("Cell Type") + labs(fill="Gene Set") + ggtitle(in_title))
		}
		# ct <- subset(ct, 
	
	}
	
	ggplot(this_dat, aes(x=cell_type, y=parameter, fill=status), outlier.shape = NA) + geom_boxplot() + coord_cartesian(ylim=c(y_low, y_high)) + ylab(ylabel) + xlab("Cell Type") + labs(fill="Gene Set") + ggtitle(in_title)
}



plot_onlyone <- function(in_mat, in_mat_lab, whi_spec, y_low, y_high, ylabel, in_title, fix_ana=FALSE){
	this_dat <- format_df(in_mat, in_mat, in_mat_lab, in_mat_lab)	
	this_dat <- this_dat[which(this_dat[,"spec"] == whi_spec),]

	gg_color_hue <- function(n) {
	  hues = seq(15, 375, length = n + 1)
	  hcl(h = hues, l = 65, c = 100)[1:n]
	}

	col_wheel <- gg_color_hue(8)

	ct_cols <- col_wheel[c(7, 3, 6, 2, 5, 1, 4, 8)]
	
	ct <- this_dat$cell_type
	ct[which(ct == "GSC_Early_spermatogonia")] <- "GSC/Early Spermatogonia"
	ct[which(ct == "Late_spermatogonia")] <- "Late Spermatogonia"
	ct[which(ct == "Early_spermatocyte")] <- "Early Spermatocyte"
	ct[which(ct == "Late_spermatocyte")] <- "Late Spermatocyte"
	ct[which(ct == "Early_spermatid")] <- "Early Spermatid"
	ct[which(ct == "Late_spermatid")] <- "Late Spermatid"
	ct[which(ct == "ananassae_spermatid")] <- "ananassae Spermatocyte"
	
	stat <- this_dat$status
	stat[which(stat == "All_Genes")] <- "All Genes"
	stat[which(stat == "Neighbor_L2")] <- "Neighbor L2"
	stat[which(stat == "Neighbor_L1")] <- "Neighbor L1"
	stat[which(stat == "Neighbor_R1")] <- "Neighbor R1"
	stat[which(stat == "Neighbor_R2")] <- "Neighbor R2"
	stat[which(stat == "X_linked")] <- "X Linked"
	stat[which(stat == "Autosomal")] <- "Autosomal"
	
	this_dat$cell_type <- factor(ct, levels=c("Somatic", "GSC/Early Spermatogonia", "Late Spermatogonia", "ananassae Spermatocyte", "Early Spermatocyte", "Late Spermatocyte", "Early Spermatid",  "Late Spermatid"))
	
	this_dat$status <- factor(stat, levels=c("Neighbor L2", "Neighbor L1", "Neighbor R1", "Neighbor R2", "All Genes", "X Linked", "Autosomal"))
	
	this_breaks <- factor(c("Somatic", "GSC/Early Spermatogonia", "Late Spermatogonia", "ananassae Spermatocyte", "Early Spermatocyte", "Late Spermatocyte", "Early Spermatid",  "Late Spermatid"))
	
	this_limits <- 1:8
	
	if(fix_ana){
		if(whi_spec %in% c("mel", "yak")){
			this_dat <- subset(this_dat, cell_type %in% c("Somatic", "GSC/Early Spermatogonia", "Late Spermatogonia", "Early Spermatocyte", "Late Spermatocyte", "Early Spermatid", "Late Spermatid"))
			this_dat$cell_type <- factor(this_dat$cell_type, levels=c("Somatic", "GSC/Early Spermatogonia", "Late Spermatogonia", "Early Spermatocyte", "Late Spermatocyte", "Early Spermatid", "Late Spermatid"))
			
			ct_cols <- ct_cols[1:7]
			this_limits <- 1:7
			
			this_breaks <- factor(c("Somatic", "GSC/Early Spermatogonia", "Late Spermatogonia", "Early Spermatocyte", "Late Spermatocyte", "Early Spermatid", "Late Spermatid"))
			
			return(ggplot(this_dat, aes(x=cell_type, y=parameter), outlier.shape = NA) + geom_boxplot(fill=ct_cols) + scale_x_discrete("cell_type", breaks=this_breaks, drop=FALSE) + coord_cartesian(ylim=c(y_low, y_high)) + ylab(ylabel) + xlab("Cell Type") + labs(fill="Gene Set") + ggtitle(in_title))
		}

		if(whi_spec == "ana"){
			this_dat <- subset(this_dat, cell_type %in% c("Somatic", "GSC/Early Spermatogonia", "Late Spermatogonia", "ananassae Spermatocyte", "Late Spermatocyte", "Early Spermatid", "Late Spermatid"))
			this_dat$cell_type <- factor(this_dat$cell_type, levels=c("Somatic", "GSC/Early Spermatogonia", "Late Spermatogonia", "ananassae Spermatocyte", "Late Spermatocyte", "Early Spermatid", "Late Spermatid"))
			
			ct_cols <- ct_cols[c(1:7)]			
			this_labs <- factor(c("Somatic", "GSC/Early Spermatogonia", "Late Spermatogonia", "ananassae Spermatocyte", "Late Spermatocyte", "Early Spermatid", "Late Spermatid"))
						
			return(ggplot(this_dat, aes(x=cell_type, y=parameter), outlier.shape = NA) + geom_boxplot(fill=ct_cols) + scale_x_discrete("cell_type", breaks=this_labs, drop=FALSE) + coord_cartesian(ylim=c(y_low, y_high)) + ylab(ylabel) + xlab("Cell Type") + labs(fill="Gene Set") + ggtitle(in_title))
		}
		# ct <- subset(ct, 
	
	}
	
	ggplot(this_dat, aes(x=cell_type, y=parameter), outlier.shape = NA) + geom_boxplot(fill=ct_cols) + scale_x_discrete("cell_type", breaks=this_breaks, drop=FALSE) + coord_cartesian(ylim=c(y_low, y_high)) + ylab(ylabel) + xlab("Cell Type") + labs(fill="Gene Set") + ggtitle(in_title)
}



strand_flipper <- function(dnt_table, sym_conv, gtf_table){
	result <- dnt_table
	
	for(i in 1:nrow(dnt_table)){
		if(dnt_table[i, "fbgn_mel"] %in% symb_id_conv){
			stnd <- gtf_table[which(symb_id_conv == dnt_table[i, "fbgn_mel"])[1], 7]
			if(stnd == "-"){
				#flippem
				result[i,c("nei_L2_sym_mel", "nei_L2_fbgn_mel", "nei_L2_fbgn_yak", "nei_L2_fbgn_ana", "nei_L1_sym_mel", "nei_L1_fbgn_mel", "nei_L1_fbgn_yak", "nei_L1_fbgn_ana")] <- dnt_table[i, c("nei_R2_sym_mel", "nei_R2_fbgn_mel", "nei_R2_fbgn_yak", "nei_R2_fbgn_ana", "nei_R1_sym_mel", "nei_R1_fbgn_mel", "nei_R1_fbgn_yak", "nei_R1_fbgn_ana")]
				result[i,c("nei_R2_sym_mel", "nei_R2_fbgn_mel", "nei_R2_fbgn_yak", "nei_R2_fbgn_ana", "nei_R1_sym_mel", "nei_R1_fbgn_mel", "nei_R1_fbgn_yak", "nei_R1_fbgn_ana")] <- dnt_table[i, c("nei_L2_sym_mel", "nei_L2_fbgn_mel", "nei_L2_fbgn_yak", "nei_L2_fbgn_ana", "nei_L1_sym_mel", "nei_L1_fbgn_mel", "nei_L1_fbgn_yak", "nei_L1_fbgn_ana")]
			}
		}
	}
	
	return(result)
}

dnt_dat_stranded <- strand_flipper(dnt_dat, symb_id_conv, mel_gtf)

dnt_L2_stnd <- gen_burst_table(symb_id_conv_rev[dnt_dat_stranded[,"nei_L2_fbgn_mel"]], summary_all_freq_mean, summary_all_size_mean, "sc_liftoff/bursting/dnt_L2_bursting_comb.csv")
dnt_L2_freq_mean_stnd <- dnt_L2_stnd[[1]]
dnt_L2_size_mean_stnd <- dnt_L2_stnd[[2]]

dnt_L1_stnd <- gen_burst_table(symb_id_conv_rev[dnt_dat_stranded[,"nei_L1_fbgn_mel"]], summary_all_freq_mean, summary_all_size_mean, "sc_liftoff/bursting/dnt_L1_bursting_comb.csv")
dnt_L1_freq_mean_stnd <- dnt_L1_stnd[[1]]
dnt_L1_size_mean_stnd <- dnt_L1_stnd[[2]]

dnt_R1_stnd <- gen_burst_table(symb_id_conv_rev[dnt_dat_stranded[,"nei_R1_fbgn_mel"]], summary_all_freq_mean, summary_all_size_mean, "sc_liftoff/bursting/dnt_R1_bursting_comb.csv")
dnt_R1_freq_mean_stnd <- dnt_R1_stnd[[1]]
dnt_R1_size_mean_stnd <- dnt_R1_stnd[[2]]

dnt_R2_stnd <- gen_burst_table(symb_id_conv_rev[dnt_dat_stranded[,"nei_R2_fbgn_mel"]], summary_all_freq_mean, summary_all_size_mean, "sc_liftoff/bursting/dnt_R2_bursting_comb.csv")
dnt_R2_freq_mean_stnd <- dnt_R2_stnd[[1]]
dnt_R2_size_mean_stnd <- dnt_R2_stnd[[2]]

#plot_onlyone <- function(in_mat, in_mat_lab, whi_spec, y_low, y_high, ylabel, in_title, fix_ana=FALSE){

if(plot_new){
	png("sc_liftoff/bursting/all_cells_burstfreq.png", width=2400*3*0.92*2, height=2400*3*0.8, res=600, pointsize=3)

	p1 <- plot_onlyone(all_burst_size_mean, "All_Genes", "mel", 0, 32, "Burst Size", "Burst Size in D. melanogaster", TRUE)
	p2 <- plot_onlyone(all_burst_size_mean, "All_Genes", "yak", 0, 85, "Burst Size", "Burst Size in D. yakuba", TRUE)
	p3 <- plot_onlyone(all_burst_size_mean, "All_Genes", "ana", 1.2, 4.75, "Burst Size", "Burst Size in D. ananassae", TRUE)

	p4 <- plot_onlyone(all_burst_freq_mean, "All_Genes", "mel", 0, 7.2, "Burst Frequency", "Burst Frequency in D. melanogaster", TRUE)
	p5 <- plot_onlyone(all_burst_freq_mean, "All_Genes", "yak", 0, 2.5, "Burst Frequency", "Burst Frequency in D. yakuba", TRUE)
	p6 <- plot_onlyone(all_burst_freq_mean, "All_Genes", "ana", 0, 6, "Burst Frequency", "Burst Frequency in D. ananassae", TRUE)

	p1 <- p1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	p2 <- p2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	p3 <- p3 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	p4 <- p4 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	p5 <- p5 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	p6 <- p6 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


	grid.arrange(p1, p4, p2, p5, p3, p6, ncol=2)

	dev.off()

	svg("sc_liftoff/bursting/all_cells_burstfreq.svg", width=22.08, height=9.6, pointsize=3)

	grid.arrange(p1, p4, p2, p5, p3, p6, ncol=2)

	dev.off()
}

extract_paired <- function(in1, in2, spec_label, ct_label){
	this_df <- format_df(in1, in2, "a", "b")
	subdf1 <- this_df[this_df[,"spec"] == spec_label,]
	subdf2 <- subdf1[subdf1[,"cell_type"] == ct_label,]
	delta1 <- subdf2[subdf2[,"status"] == "a", "parameter"]
	delta2 <- subdf2[subdf2[,"status"] == "b", "parameter"]

	return(list(delta1, delta2))
}

run_pval_mannwhitney <- function(in1, in2, spec_label, ct_label){
	delta <- extract_paired(in1, in2, spec_label, ct_label)

	return(wilcox.test(delta[[1]], delta[[2]])$p.value)
}

run_pval_mannwhitney(all_burst_size_mean, dnt_L1_size_mean, "mel", "Early_spermatocyte")


#cols will be celltype, rows will be all, L2...R2, the collate mel, yak, ana, make matrix for p.vals, and N
cell_type_vec <- c("Somatic", "GSC_Early_spermatogonia", "Late_spermatogonia", "Early_spermatocyte", "Late_spermatocyte", "Early_spermatid",  "Late_spermatid", "ananassae_spermatid")
spec_vec <- c("mel", "yak", "ana")
dat_size_list <- list(dnt_L2_size_mean, dnt_L1_size_mean, dnt_R1_size_mean, dnt_R2_size_mean)
dat_freq_list <- list(dnt_L2_freq_mean, dnt_L1_freq_mean, dnt_R1_freq_mean, dnt_R2_freq_mean)

row_labs <- c("mel_L2", "mel_L1", "mel_R1", "mel_R2", "yak_L2", "yak_L1", "yak_R1", "yak_R2", "ana_L2", "ana_L1", "ana_R1", "ana_R2")

gen_row_by_ct <- function(spec_lab, in_dat_ref, in_dat, ct_list, cat_ret){
	result <- rep(0, length(ct_list))
	if(cat_ret == "pval") result <- unlist(lapply(ct_list, function(x) run_pval_mannwhitney(in_dat_ref, in_dat, spec_lab, x)))
	if(cat_ret == "N_main") result <- unlist(lapply(ct_list, function(x) length(extract_paired(in_dat_ref, in_dat, spec_lab, x)[[1]])))
	if(cat_ret == "N_exp") result <- unlist(lapply(ct_list, function(x) length(extract_paired(in_dat_ref, in_dat, spec_lab, x)[[2]])))
	
	return(result)
}

gen_rows_specs <- function(spec_lab, in_dat_ref, in_dat_list, ct_list, cat_ret){
	result <- lapply(in_dat_list, function(x) gen_row_by_ct(spec_lab, in_dat_ref, x, ct_list, cat_ret))
	return(do.call("rbind", result))
}

gen_rows_all <- function(spec_lab_list, in_dat_ref, in_dat_list, ct_list, cat_ret){
	result <- lapply(spec_lab_list, function(x) gen_rows_specs(x, in_dat_ref, in_dat_list, ct_list, cat_ret))
	return(do.call("rbind", result))
}

mtx_size_pval <- gen_rows_all(spec_vec, all_burst_size_mean, dat_size_list, cell_type_vec, "pval")
mtx_size_N <- gen_rows_all(spec_vec, all_burst_size_mean, dat_size_list, cell_type_vec, "N_exp")
mtx_freq_pval <- gen_rows_all(spec_vec, all_burst_freq_mean, dat_freq_list, cell_type_vec, "pval")
mtx_freq_N <- gen_rows_all(spec_vec, all_burst_freq_mean, dat_freq_list, cell_type_vec, "N_exp")

mtx_size_N_ref <- gen_rows_all(spec_vec, all_burst_size_mean, dat_size_list, cell_type_vec, "N_main")
mtx_freq_N_ref <- gen_rows_all(spec_vec, all_burst_freq_mean, dat_freq_list, cell_type_vec, "N_main")

rownames(mtx_size_pval) <- row_labs
rownames(mtx_size_N) <- row_labs
rownames(mtx_freq_pval) <- row_labs
rownames(mtx_freq_N) <- row_labs

rownames(mtx_size_N_ref) <- row_labs
rownames(mtx_freq_N_ref) <- row_labs

colnames(mtx_size_pval) <- cell_type_vec
colnames(mtx_size_N) <- cell_type_vec
colnames(mtx_freq_pval) <- cell_type_vec
colnames(mtx_freq_N) <- cell_type_vec

colnames(mtx_size_N_ref) <- cell_type_vec
colnames(mtx_freq_N_ref) <- cell_type_vec

mtx_size_N_ref <- do.call("rbind", list(mtx_size_N_ref[1,], mtx_size_N_ref[5,], mtx_size_N_ref[9,]))
mtx_freq_N_ref <- do.call("rbind", list(mtx_freq_N_ref[1,], mtx_freq_N_ref[5,], mtx_freq_N_ref[9,]))

rownames(mtx_size_N_ref) <- c("mel", "yak", "ana")
rownames(mtx_freq_N_ref) <- c("mel", "yak", "ana")

mtx_size_N_ref_rep <- do.call("rbind", unlist(list(rep(list(mtx_size_N_ref[1,]),4), rep(list(mtx_size_N_ref[2,]), 4), rep(list(mtx_size_N_ref[3,]), 4)), recursive=F))
mtx_freq_N_ref_rep <- do.call("rbind", unlist(list(rep(list(mtx_freq_N_ref[1,]),4), rep(list(mtx_freq_N_ref[2,]), 4), rep(list(mtx_freq_N_ref[3,]), 4)), recursive=F))
rownames(mtx_size_N_ref_rep) <- row_labs
rownames(mtx_freq_N_ref_rep) <- row_labs

write.csv(mtx_size_pval, file="sc_liftoff/bursting/mtx_size_pval.csv")
write.csv(mtx_freq_pval, file="sc_liftoff/bursting/mtx_freq_pval.csv")
write.csv(mtx_size_N, file="sc_liftoff/bursting/mtx_size_N.csv")
write.csv(mtx_freq_N, file="sc_liftoff/bursting/mtx_freq_N.csv")

write.csv(mtx_size_N_ref, file="sc_liftoff/bursting/mtx_size_N_ref.csv")
write.csv(mtx_freq_N_ref, file="sc_liftoff/bursting/mtx_freq_N_ref.csv")

write.csv(all_burst_size_mean, file="sc_liftoff/bursting/all_burst_size_mean.csv")
write.csv(all_burst_freq_mean, file="sc_liftoff/bursting/all_burst_freq_mean.csv")
write.csv(dnt_L2_size_mean, file="sc_liftoff/bursting/dnt_L2_size_mean.csv")
write.csv(dnt_L1_size_mean, file="sc_liftoff/bursting/dnt_L1_size_mean.csv") 
write.csv(dnt_R1_size_mean, file="sc_liftoff/bursting/dnt_R1_size_mean.csv") 
write.csv(dnt_R2_size_mean, file="sc_liftoff/bursting/dnt_R2_size_mean.csv") 
write.csv(dnt_L2_freq_mean, file="sc_liftoff/bursting/dnt_L2_freq_mean.csv")
write.csv(dnt_L1_freq_mean, file="sc_liftoff/bursting/dnt_L1_freq_mean.csv") 
write.csv(dnt_R1_freq_mean, file="sc_liftoff/bursting/dnt_R1_freq_mean.csv") 
write.csv(dnt_R2_freq_mean, file="sc_liftoff/bursting/dnt_R2_freq_mean.csv") 


mtx_comb_list <- unlist(lapply(1:8, function(x) return(list(mtx_size_pval[,x], mtx_freq_pval[,x], mtx_size_N[,x], mtx_size_N_ref_rep[,x]))), recursive=F)
mtx_comb_colnames <- unlist(lapply(cell_type_vec, function(x) return(c(paste(x, "_size_p", sep=""), paste(x, "_freq_p", sep=""), paste(x, "_N_genes", sep=""), paste(x, "_N_allgenes", sep="")))))
 
mtx_comb <- do.call("cbind", mtx_comb_list)
colnames(mtx_comb) <- mtx_comb_colnames
rownames(mtx_comb) <- rownames(mtx_size_pval)



write.csv(mtx_comb, file="sc_liftoff/bursting/mtx_comb.csv")


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

bh_size <- benj_hoch(mtx_size_pval, 0.05)
bh_freq <- benj_hoch(mtx_freq_pval, 0.05)

write.csv(bh_size, file="sc_liftoff/bursting/bh_size_dnt.csv")
write.csv(bh_freq, file="sc_liftoff/bursting/bh_freq_dnt.csv")

add_N_labs <- function(gplot, xs, ys, Ns){
	for(i in 1: length(xs)){
		gplot <- gplot + geom_text(x=xs[i], y=ys[i], label=paste("(N=", Ns[i], ")", sep=""))
	}

	return(gplot)
}

get_pval_stars <- function(pval_vec){
	result <- as.character(pval_vec)
	for(i in 1:length(pval_vec)){
		result[i] <- "NS"
		if(pval_vec[i] < 0.05) result[i] <- "*"
		if(pval_vec[i] < 0.01) result[i] <- "**"
		if(pval_vec[i] < 0.001) result[i] <- "***"
		if(pval_vec[i] < 0.0001) result[i] <- "****"
	}

	return(result)
}


plot_neigh <- function(dat_1, dat_2, lab_1, lab_2, lab_spec, y_low, y_high, y_signif, y_locn, lab_ylab, lab_mainlab, ana_fix, pv_vals, n_vals_1, n_vals_2){
	idxs <- 1:length(n_vals_1)
	xs_idx <- 0:(length(n_vals_1)-1)
	if(ana_fix && lab_spec %in% c("mel", "yak")){
		idxs <- 1:7
		xs_idx <- 0:6
	}
	if(ana_fix && lab_spec == "ana"){
		idxs <- c(1:3, 8, 5:7)
		xs_idx <- 0:6
	}
	
	vec_n_1 <- n_vals_1[idxs]
	vec_n_2 <- n_vals_2[idxs]
	vec_pv <- pv_vals[idxs]
	
	output <- plot_ind(dat_1, dat_2, lab_1, lab_2, lab_spec, y_low, y_high, lab_ylab, lab_mainlab, ana_fix)
	output <- output + geom_signif(y_position=rep(y_signif, length(idxs)), xmin=xs_idx + 0.8, xmax=xs_idx + 1.2, annotation=get_pval_stars(vec_pv), tip_length=0, vjust=2) 
	output <- add_N_labs(output, xs_idx + 0.75, rep(y_locn, length(idxs)), vec_n_1)
	output <- add_N_labs(output, xs_idx + 1.25, rep(y_locn, length(idxs)), vec_n_2)

	return(output)
}


plot_neigh1 <- function(dat_1, dat_2, lab_1, lab_2, lab_spec, y_low, y_high, y_signif, y_locn, lab_ylab, lab_mainlab, ana_fix, pv_vals, n_vals_1, n_vals_2){
	idxs <- 1:length(n_vals_1)
	xs_idx <- 0:(length(n_vals_1)-1)
	if(ana_fix && lab_spec %in% c("mel", "yak")){
		idxs <- 1:7
		xs_idx <- 0:6
	}
	if(ana_fix && lab_spec == "ana"){
		idxs <- c(1:3, 8, 5:7)
		xs_idx <- 0:6
	}
	
	vec_n_1 <- n_vals_1[idxs]
	vec_n_2 <- n_vals_2[idxs]
	vec_pv <- pv_vals[idxs]
	
	output <- plot_ind1(dat_1, dat_2, lab_1, lab_2, lab_spec, y_low, y_high, lab_ylab, lab_mainlab, ana_fix)
	output <- output + geom_signif(y_position=rep(y_signif, length(idxs)), xmin=xs_idx + 0.8, xmax=xs_idx + 1.2, annotation=get_pval_stars(vec_pv), tip_length=0, vjust=2) 
	output <- add_N_labs(output, xs_idx + 0.75, rep(y_locn, length(idxs)), vec_n_1)
	output <- add_N_labs(output, xs_idx + 1.25, rep(y_locn, length(idxs)), vec_n_2)

	return(output)
}


plot_neigh_mult <- function(data_1, data_2_list, label_1, label_2_vec, label_species, y_lbound, y_hbound, y_sig, y_locaN, label_y, label_main, bool_fix_ana, pv_mat, mtx_nvals_1, n_vals_ref){
	neigh_L2 <- plot_neigh(data_1, data_2_list[[1]], label_1, label_2_vec[1], label_species, y_lbound, y_hbound, y_sig, y_locaN, label_y, label_main[1], bool_fix_ana, pv_mat[1,], mtx_nvals_1[1,], n_vals_ref)
	neigh_L2 <- neigh_L2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	
	neigh_L1 <- plot_neigh(data_1, data_2_list[[2]], label_1, label_2_vec[2], label_species, y_lbound, y_hbound, y_sig, y_locaN, label_y, label_main[2], bool_fix_ana, pv_mat[2,], mtx_nvals_1[2,], n_vals_ref)
	neigh_L1 <- neigh_L1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	
	neigh_R1 <- plot_neigh(data_1, data_2_list[[3]], label_1, label_2_vec[3], label_species, y_lbound, y_hbound, y_sig, y_locaN, label_y, label_main[3], bool_fix_ana, pv_mat[3,], mtx_nvals_1[3,], n_vals_ref)
	neigh_R1 <- neigh_R1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	
	neigh_R2 <- plot_neigh(data_1, data_2_list[[4]], label_1, label_2_vec[1], label_species, y_lbound, y_hbound, y_sig, y_locaN, label_y, label_main[4], bool_fix_ana, pv_mat[4,], mtx_nvals_1[4,], n_vals_ref)
	neigh_R2 <- neigh_R2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

	return(list(neigh_L2, neigh_L1, neigh_R1, neigh_R2))
}

if(plot_new){

	png("sc_liftoff/bursting/comb_size_all.png", width=2400*3.5*3.90, height=2400*3, res=600, pointsize=3)

	p_size_mel_list <- plot_neigh_mult(all_burst_size_mean, list(dnt_L2_size_mean, dnt_L1_size_mean, dnt_R1_size_mean, dnt_R2_size_mean), "All_Genes", c("Neighbor_L2", "Neighbor_L1", "Neighbor_R1", "Neighbor_R2"), "mel", -11*45/75, 45, -2.5*45/75, -11*45/75, "Burst Size", unlist(lapply(c("L2", "L1", "R1", "R2"), function(x) paste("Burst Size in D. melanogaster, ", x, " Neighbor", sep=""))), TRUE, mtx_size_pval[1:4,], mtx_size_N[1:4,], mtx_size_N_ref[1,])

	p_size_yak_list <- plot_neigh_mult(all_burst_size_mean, list(dnt_L2_size_mean, dnt_L1_size_mean, dnt_R1_size_mean, dnt_R2_size_mean), "All_Genes", c("Neighbor_L2", "Neighbor_L1", "Neighbor_R1", "Neighbor_R2"), "yak", -11*100/75, 100, -2.5*100/75, -11*100/75, "Burst Size", unlist(lapply(c("L2", "L1", "R1", "R2"), function(x) paste("Burst Size in D. yakuba, ", x, " Neighbor", sep=""))), TRUE, mtx_size_pval[1:4 + 1*4,], mtx_size_N[1:4 + 1*4,], mtx_size_N_ref[2,])

	p_size_ana_list <- plot_neigh_mult(all_burst_size_mean, list(dnt_L2_size_mean, dnt_L1_size_mean, dnt_R1_size_mean, dnt_R2_size_mean), "All_Genes", c("Neighbor_L2", "Neighbor_L1", "Neighbor_R1", "Neighbor_R2"), "ana", -11*12/75, 12, -2.5*15/75, -11*12/75, "Burst Size", unlist(lapply(c("L2", "L1", "R1", "R2"), function(x) paste("Burst Size in D. ananassae, ", x, " Neighbor", sep=""))), TRUE, mtx_size_pval[1:4 + 2*4,], mtx_size_N[1:4 + 2*4,], mtx_size_N_ref[3,])

	grid.arrange(p_size_mel_list[[1]], p_size_mel_list[[2]], p_size_mel_list[[3]], p_size_mel_list[[4]], p_size_yak_list[[1]], p_size_yak_list[[2]], p_size_yak_list[[3]], p_size_yak_list[[4]], p_size_ana_list[[1]], p_size_ana_list[[2]], p_size_ana_list[[3]], p_size_ana_list[[4]], ncol=4)

	dev.off()

	svg("sc_liftoff/bursting/comb_size_all_mel.svg", width=4*3.5*3.90, height=4*1, pointsize=3)
	grid.arrange(p_size_mel_list[[1]], p_size_mel_list[[2]], p_size_mel_list[[3]], p_size_mel_list[[4]], ncol=4)

	dev.off()

	svg("sc_liftoff/bursting/comb_size_all_yak.svg", width=4*3.5*3.90, height=4*1, pointsize=3)
	grid.arrange(p_size_yak_list[[1]], p_size_yak_list[[2]], p_size_yak_list[[3]], p_size_yak_list[[4]], ncol=4)

	dev.off()

	svg("sc_liftoff/bursting/comb_size_all_ana.svg", width=4*3.5*3.90, height=4*1, pointsize=3)
	grid.arrange(p_size_ana_list[[1]], p_size_ana_list[[2]], p_size_ana_list[[3]], p_size_ana_list[[4]], ncol=4)

	dev.off()


	png("sc_liftoff/bursting/comb_freq_all.png", width=2400*3.5*3.90, height=2400*3, res=600, pointsize=3)

	p_freq_mel_list <- plot_neigh_mult(all_burst_freq_mean, list(dnt_L2_freq_mean, dnt_L1_freq_mean, dnt_R1_freq_mean, dnt_R2_freq_mean), "All_Genes", c("Neighbor_L2", "Neighbor_L1", "Neighbor_R1", "Neighbor_R2"), "mel", -11*10/75, 10, -2.5*10/75, -11*10/75, "Burst Frequency", unlist(lapply(c("L2", "L1", "R1", "R2"), function(x) paste("Burst Frequency in D. melanogaster, ", x, " Neighbor", sep=""))), TRUE, mtx_freq_pval[1:4,], mtx_freq_N[1:4,], mtx_freq_N_ref[1,])

	p_freq_yak_list <- plot_neigh_mult(all_burst_freq_mean, list(dnt_L2_freq_mean, dnt_L1_freq_mean, dnt_R1_freq_mean, dnt_R2_freq_mean), "All_Genes", c("Neighbor_L2", "Neighbor_L1", "Neighbor_R1", "Neighbor_R2"), "yak", -11*4/75, 4, -2.5*4/75, -11*4/75, "Burst Frequency", unlist(lapply(c("L2", "L1", "R1", "R2"), function(x) paste("Burst Frequency in D. yakuba, ", x, " Neighbor", sep=""))), TRUE, mtx_freq_pval[1:4 + 1*4,], mtx_freq_N[1:4 + 1*4,], mtx_freq_N_ref[2,])

	p_freq_ana_list <- plot_neigh_mult(all_burst_freq_mean, list(dnt_L2_freq_mean, dnt_L1_freq_mean, dnt_R1_freq_mean, dnt_R2_freq_mean), "All_Genes", c("Neighbor_L2", "Neighbor_L1", "Neighbor_R1", "Neighbor_R2"), "ana", -11*7.5/75, 7.5, -2.5*7.5/75, -11*7.5/75, "Burst Frequency", unlist(lapply(c("L2", "L1", "R1", "R2"), function(x) paste("Burst Frequency in D. ananassae, ", x, " Neighbor", sep=""))), TRUE, mtx_freq_pval[1:4 + 2*4,], mtx_freq_N[1:4 + 2*4,], mtx_freq_N_ref[3,])

	grid.arrange(p_freq_mel_list[[1]], p_freq_mel_list[[2]], p_freq_mel_list[[3]], p_freq_mel_list[[4]], p_freq_yak_list[[1]], p_freq_yak_list[[2]], p_freq_yak_list[[3]], p_freq_yak_list[[4]], p_freq_ana_list[[1]], p_freq_ana_list[[2]], p_freq_ana_list[[3]], p_freq_ana_list[[4]], ncol=4)

	dev.off()


	svg("sc_liftoff/bursting/comb_freq_all_mel.svg", width=4*3.5*3.90, height=4*1, pointsize=3)
	grid.arrange(p_freq_mel_list[[1]], p_freq_mel_list[[2]], p_freq_mel_list[[3]], p_freq_mel_list[[4]], ncol=4)

	dev.off()

	svg("sc_liftoff/bursting/comb_freq_all_yak.svg", width=4*3.5*3.90, height=4*1, pointsize=3)
	grid.arrange(p_freq_yak_list[[1]], p_freq_yak_list[[2]], p_freq_yak_list[[3]], p_freq_yak_list[[4]], ncol=4)

	dev.off()

	svg("sc_liftoff/bursting/comb_freq_all_ana.svg", width=4*3.5*3.90, height=4*1, pointsize=3)
	grid.arrange(p_freq_ana_list[[1]], p_freq_ana_list[[2]], p_freq_ana_list[[3]], p_freq_ana_list[[4]], ncol=4)

	dev.off()
}
#run dedup.py
#run svgo on all svgs

extract_paired <- function(in1, in2, spec_label, ct_label){
	this_df <- format_df(in1, in2, "a", "b")
	subdf1 <- this_df[this_df[,"spec"] == spec_label,]
	subdf2 <- subdf1[subdf1[,"cell_type"] == ct_label,]
	delta1 <- subdf2[subdf2[,"status"] == "a", "parameter"]
	delta2 <- subdf2[subdf2[,"status"] == "b", "parameter"]

	return(list(delta1, delta2))
}

paired_mel_L1_e_spermatocyte <- extract_paired(all_burst_size_mean, dnt_L1_size_mean, "mel", "Early_spermatocyte")
p_m_l1_es <- data.frame(vals=c(paired_mel_L1_e_spermatocyte[[1]], paired_mel_L1_e_spermatocyte[[2]]), sourc=factor(c(rep("All Genes", length(paired_mel_L1_e_spermatocyte[[1]])), rep("L1 Neighbor", length(paired_mel_L1_e_spermatocyte[[2]])) ), levels=c("L1 Neighbor", "All Genes")))

paired_yak_L1_e_spermatocyte <- extract_paired(all_burst_size_mean, dnt_L1_size_mean, "yak", "Early_spermatocyte")
p_y_l1_es <- data.frame(vals=c(paired_yak_L1_e_spermatocyte[[1]], paired_yak_L1_e_spermatocyte[[2]]), sourc=factor(c(rep("All Genes", length(paired_yak_L1_e_spermatocyte[[1]])), rep("L1 Neighbor", length(paired_yak_L1_e_spermatocyte[[2]])) ), levels=c("L1 Neighbor", "All Genes")))

paired_yak_L1_l_spermatocyte <- extract_paired(all_burst_size_mean, dnt_L1_size_mean, "yak", "Late_spermatocyte")
p_y_l1_ls <- data.frame(vals=c(paired_yak_L1_l_spermatocyte[[1]], paired_yak_L1_l_spermatocyte[[2]]), sourc=factor(c(rep("All Genes", length(paired_yak_L1_l_spermatocyte[[1]])), rep("L1 Neighbor", length(paired_yak_L1_l_spermatocyte[[2]]))), levels=c("L1 Neighbor", "All Genes")))

if(plot_new){
	png("sc_liftoff/bursting/comb_sigs_dens.png", width=2400*3, height=2400*3, res=600, pointsize=3)
	p1 <- ggplot(p_m_l1_es, aes(x=vals, fill=sourc)) + geom_density(alpha=0.4, outline.type = "full", adjust=1/2, n=2^13) + scale_x_continuous(limits=c(0,575)) + coord_cartesian(xlim=c(0, 70))
	p2 <- ggplot(p_y_l1_es, aes(x=vals, fill=sourc)) + geom_density(alpha=0.4, outline.type = "full", adjust=1/2, n=2^13) + scale_x_continuous(limits=c(0,575)) + coord_cartesian(xlim=c(0, 35))
	p3 <- ggplot(p_y_l1_ls, aes(x=vals, fill=sourc)) + geom_density(alpha=0.4, outline.type = "full", adjust=1/2, n=2^13) + scale_x_continuous(limits=c(0,575)) + coord_cartesian(xlim=c(0, 15))
	grid.arrange(p1, p2, p3, ncol=1)
	dev.off()


	svg("sc_liftoff/bursting/comb_sigs_dens.svg", width=4*3, height=4*3, pointsize=3)
	p1 <- ggplot(p_m_l1_es, aes(x=vals, fill=sourc)) + geom_density(alpha=0.4, outline.type = "full", adjust=1/2, n=2^13) + scale_x_continuous(limits=c(0,575)) + coord_cartesian(xlim=c(0, 70))
	p1 <- p1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

	p2 <- ggplot(p_y_l1_es, aes(x=vals, fill=sourc)) + geom_density(alpha=0.4, outline.type = "full", adjust=1/2, n=2^13) + scale_x_continuous(limits=c(0,575)) + coord_cartesian(xlim=c(0, 35))
	p2 <- p2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

	p3 <- ggplot(p_y_l1_ls, aes(x=vals, fill=sourc)) + geom_density(alpha=0.4, outline.type = "full", adjust=1/2, n=2^13) + scale_x_continuous(limits=c(0,575)) + coord_cartesian(xlim=c(0, 15))
	p3 <- p3 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

	grid.arrange(p1, p2, p3, ncol=1)
	dev.off()


	png("sc_liftoff/bursting/comb_size_zoom.png", width=2400*3.5*3.90/2, height=2400*3, res=600, pointsize=3)

	p_size_mel_list1 <- plot_neigh_mult(all_burst_size_mean, list(dnt_L2_size_mean, dnt_L1_size_mean, dnt_R1_size_mean, dnt_R2_size_mean), "All_Genes", c("Neighbor_L2", "Neighbor_L1", "Neighbor_R1", "Neighbor_R2"), "mel", -11*40/75, 40, -2.5*40/75, -11*45/75, "Burst Size", unlist(lapply(c("L2", "L1", "R1", "R2"), function(x) paste("Burst Size in D. melanogaster, ", x, " Neighbor", sep=""))), TRUE, mtx_size_pval[1:4,], mtx_size_N[1:4,], mtx_size_N_ref[1,])

	p_size_yak_list1 <- plot_neigh_mult(all_burst_size_mean, list(dnt_L2_size_mean, dnt_L1_size_mean, dnt_R1_size_mean, dnt_R2_size_mean), "All_Genes", c("Neighbor_L2", "Neighbor_L1", "Neighbor_R1", "Neighbor_R2"), "yak", -11*23.5/75, 23.5, -2.5*23.5/75, -11*23.5/75, "Burst Size", unlist(lapply(c("L2", "L1", "R1", "R2"), function(x) paste("Burst Size in D. yakuba, ", x, " Neighbor", sep=""))), TRUE, mtx_size_pval[1:4 + 1*4,], mtx_size_N[1:4 + 1*4,], mtx_size_N_ref[2,])

	p_size_yak_list2 <- plot_neigh_mult(all_burst_size_mean, list(dnt_L2_size_mean, dnt_L1_size_mean, dnt_R1_size_mean, dnt_R2_size_mean), "All_Genes", c("Neighbor_L2", "Neighbor_L1", "Neighbor_R1", "Neighbor_R2"), "yak", -11*10/75, 10, -2.5*10/75, -11*10/75, "Burst Size", unlist(lapply(c("L2", "L1", "R1", "R2"), function(x) paste("Burst Size in D. yakuba, ", x, " Neighbor", sep=""))), TRUE, mtx_size_pval[1:4 + 1*4,], mtx_size_N[1:4 + 1*4,], mtx_size_N_ref[2,])

	grid.arrange(p_size_mel_list1[[2]], p_size_yak_list1[[2]], p_size_yak_list2[[2]], ncol=1)

	dev.off()

	svg("sc_liftoff/bursting/comb_size_zoom.svg", width=27.3, height=12, pointsize=3)

	grid.arrange(p_size_mel_list1[[2]], p_size_yak_list1[[2]], p_size_yak_list2[[2]], ncol=1)

	dev.off()

	pdf("sc_liftoff/bursting/comb_size_zoom.pdf", width=27.3, height=12, pointsize=3)

	grid.arrange(p_size_mel_list1[[2]], p_size_yak_list1[[2]], p_size_yak_list2[[2]], ncol=1)

	dev.off()
}

#X:A
mel_symb_autosome <- unlist(lapply(mel_gtf[which(mel_gtf[,1] %in% c("2L", "2R", "3L", "3R", "4")),9], function(x) extract_genesymbol(x)))
mel_symb_x <- unlist(lapply(mel_gtf[which(mel_gtf[,1] == "X"),9], function(x) extract_genesymbol(x)))

mel_symb_autosome_pcg <- intersect(mel_symb_autosome, mel_symb_genes)
mel_symb_x_pcg <- intersect(mel_symb_x, mel_symb_genes)

autosome_burst <- gen_burst_table(mel_symb_autosome_pcg, summary_all_freq_mean, summary_all_size_mean, "sc_liftoff/bursting/autosome_bursting_comb.csv")
autosome_burst_freq_mean <- autosome_burst[[1]]
autosome_burst_size_mean <- autosome_burst[[2]]

x_burst <- gen_burst_table(mel_symb_x_pcg, summary_all_freq_mean, summary_all_size_mean, "sc_liftoff/bursting/x_bursting_comb.csv")
x_burst_freq_mean <- x_burst[[1]]
x_burst_size_mean <- x_burst[[2]]

png("sc_liftoff/bursting/burst_size_mel_x_to_a.png", width=2400*3.5, height=2400, res=600, pointsize=3)
p_xa_size_mel <- plot_ind(autosome_burst_size_mean, x_burst_size_mean, "Autosomal", "X_linked", "mel", 0, 40, "Burst Size", "Burst Size in D. melanogaster, X:A", TRUE)
p_xa_size_mel
dev.off()

png("sc_liftoff/bursting/burst_freq_mel_x_to_a.png", width=2400*3.5, height=2400, res=600, pointsize=3)
p_xa_freq_mel <- plot_ind(autosome_burst_freq_mean, x_burst_freq_mean, "Autosomal", "X_linked", "mel", 0, 10, "Burst Frequency", "Burst Frequency in D. melanogaster, X:A", TRUE)
p_xa_freq_mel
dev.off()

png("sc_liftoff/bursting/burst_size_yak_x_to_a.png", width=2400*3.5, height=2400, res=600, pointsize=3)
p_xa_size_yak <- plot_ind(autosome_burst_size_mean, x_burst_size_mean, "Autosomal", "X_linked", "yak", 0, 85, "Burst Size", "Burst Size in D. yakuba, X:A", TRUE)
p_xa_size_yak
dev.off()

png("sc_liftoff/bursting/burst_freq_yak_x_to_a.png", width=2400*3.5, height=2400, res=600, pointsize=3)
p_xa_freq_yak <- plot_ind(autosome_burst_freq_mean, x_burst_freq_mean, "Autosomal", "X_linked", "yak", 0, 4, "Burst Frequency", "Burst Frequency in D. yakuba, X:A", TRUE)
p_xa_freq_yak 
dev.off()

png("sc_liftoff/bursting/burst_size_ana_x_to_a.png", width=2400*3.5, height=2400, res=600, pointsize=3)
p_xa_size_ana <- plot_ind(autosome_burst_size_mean, x_burst_size_mean, "Autosomal", "X_linked", "ana", 0, 6, "Burst Size", "Burst Size in D. ananassae, X:A", TRUE)
p_xa_size_ana
dev.off()

png("sc_liftoff/bursting/burst_freq_ana_x_to_a.png", width=2400*3.5, height=2400, res=600, pointsize=3)
p_xa_freq_ana <- plot_ind(autosome_burst_freq_mean, x_burst_freq_mean, "Autosomal", "X_linked", "ana", 0, 7.5, "Burst Frequency", "Burst Frequency in D. ananassae, X:A", TRUE)
p_xa_freq_ana
dev.off()

	if(plot_new){
	png("sc_liftoff/bursting/comb_xa_all.png", width=2400*3.5*3.90/2, height=2400*3, res=600, pointsize=3)

	grid.arrange(p_xa_size_mel, p_xa_freq_mel, p_xa_size_yak, p_xa_freq_yak, p_xa_size_ana, p_xa_freq_ana, ncol=2)

	dev.off()

	svg("sc_liftoff/bursting/comb_xa_all.svg", width=27.3, height=12, pointsize=3)

	grid.arrange(p_xa_size_mel, p_xa_freq_mel, p_xa_size_yak, p_xa_freq_yak, p_xa_size_ana, p_xa_freq_ana, ncol=2)

	dev.off()
}

dat_size_x_list <- list(x_burst_size_mean)
dat_freq_x_list <- list(x_burst_freq_mean)
dat_size_autosome_list <- list(autosome_burst_size_mean)
dat_freq_autosome_list <- list(autosome_burst_freq_mean)

row_labs_xa <- c("mel_XA", "yak_XA", "ana_XA")

mtx_size_pval_xa_list <- lapply(spec_vec, function(x) gen_row_by_ct(x, autosome_burst_size_mean, x_burst_size_mean, cell_type_vec, "pval"))
mtx_size_N_xa_list <- lapply(spec_vec, function(x) gen_row_by_ct(x, autosome_burst_size_mean, x_burst_size_mean, cell_type_vec, "N_exp"))
mtx_size_Nref_xa_list <- lapply(spec_vec, function(x) gen_row_by_ct(x, autosome_burst_size_mean, x_burst_size_mean, cell_type_vec, "N_main"))

mtx_freq_pval_xa_list <- lapply(spec_vec, function(x) gen_row_by_ct(x, autosome_burst_freq_mean, x_burst_freq_mean, cell_type_vec, "pval"))
mtx_freq_N_xa_list <- lapply(spec_vec, function(x) gen_row_by_ct(x, autosome_burst_freq_mean, x_burst_freq_mean, cell_type_vec, "N_exp"))
mtx_freq_Nref_xa_list <- lapply(spec_vec, function(x) gen_row_by_ct(x, autosome_burst_freq_mean, x_burst_freq_mean, cell_type_vec, "N_main"))

mtx_size_pval_xa <- do.call("rbind", mtx_size_pval_xa_list)
mtx_size_N_xa <- do.call("rbind", mtx_size_N_xa_list)
mtx_size_Nref_xa <- do.call("rbind", mtx_size_Nref_xa_list)

mtx_freq_pval_xa <- do.call("rbind", mtx_freq_pval_xa_list)
mtx_freq_N_xa <- do.call("rbind", mtx_freq_N_xa_list)
mtx_freq_Nref_xa <- do.call("rbind", mtx_freq_Nref_xa_list)


rownames(mtx_size_pval_xa) <- row_labs_xa
rownames(mtx_size_N_xa) <- row_labs_xa
rownames(mtx_size_Nref_xa) <- row_labs_xa
rownames(mtx_freq_pval_xa) <- row_labs_xa
rownames(mtx_freq_N_xa) <- row_labs_xa
rownames(mtx_freq_Nref_xa) <- row_labs_xa

colnames(mtx_size_pval_xa) <- cell_type_vec
colnames(mtx_size_N_xa) <- cell_type_vec
colnames(mtx_size_Nref_xa) <- cell_type_vec
colnames(mtx_freq_pval_xa) <- cell_type_vec
colnames(mtx_freq_N_xa) <- cell_type_vec
colnames(mtx_freq_Nref_xa) <- cell_type_vec

if(plot_new){
	write.csv(mtx_size_pval_xa, file="sc_liftoff/bursting/mtx_size_pval_xa.csv")
	write.csv(mtx_freq_pval_xa, file="sc_liftoff/bursting/mtx_freq_pval_xa.csv")
	write.csv(mtx_size_N_xa, file="sc_liftoff/bursting/mtx_size_N_xa.csv")
	write.csv(mtx_freq_N_xa, file="sc_liftoff/bursting/mtx_freq_N_xa.csv")
	write.csv(mtx_size_Nref_xa, file="sc_liftoff/bursting/mtx_size_Nref_xa.csv")
	write.csv(mtx_freq_Nref_xa, file="sc_liftoff/bursting/mtx_freq_Nref_xa.csv")
}

bh_size_xa <- benj_hoch(mtx_size_pval_xa, 0.05)
bh_freq_xa <- benj_hoch(mtx_freq_pval_xa, 0.05)

if(plot_new){
	write.csv(bh_size_xa, file="sc_liftoff/bursting/bh_size_xa.csv")
	write.csv(bh_freq_xa, file="sc_liftoff/bursting/bh_freq_xa.csv")


	png("sc_liftoff/bursting/comb_xa_all_stats.png", width=2400*3.5*3.90/2, height=2400*3, res=600, pointsize=3)

	p_xa_size_mel_stat <- plot_neigh(autosome_burst_size_mean, x_burst_size_mean, "Autosomal", "X_linked", "mel", -11*32/75, 32, -2.5*32/75, -11*32/75, "Burst Size", "Burst Size in D. melanogaster, X:A", TRUE, mtx_size_pval_xa[1,], mtx_size_N_xa[1,], mtx_size_Nref_xa[1,])
	p_xa_freq_mel_stat <- plot_neigh(autosome_burst_freq_mean, x_burst_freq_mean, "Autosomal", "X_linked", "mel", -11*7.2/75, 7.2, -2.5*7.2/75, -11*7.2/75, "Burst Frequency", "Burst Frequency in D. melanogaster, X:A", TRUE, mtx_freq_pval_xa[1,], mtx_freq_N_xa[1,], mtx_freq_Nref_xa[1,])

	p_xa_size_mel_stat <- p_xa_size_mel_stat + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	p_xa_freq_mel_stat <- p_xa_freq_mel_stat + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

	p_xa_size_yak_stat <- plot_neigh(autosome_burst_size_mean, x_burst_size_mean, "Autosomal", "X_linked", "yak", -11*85/75, 85, -2.5*85/75, -11*85/75, "Burst Size", "Burst Size in D. yakuba, X:A", TRUE, mtx_size_pval_xa[2,], mtx_size_N_xa[2,], mtx_size_Nref_xa[2,])
	p_xa_freq_yak_stat <- plot_neigh(autosome_burst_freq_mean, x_burst_freq_mean, "Autosomal", "X_linked", "yak", -11*2.5/75, 2.5, -2.5*2.5/75, -11*2.5/75, "Burst Frequency", "Burst Frequency in D. yakuba, X:A", TRUE, mtx_freq_pval_xa[2,], mtx_freq_N_xa[2,], mtx_freq_Nref_xa[2,])

	p_xa_size_yak_stat <- p_xa_size_yak_stat + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	p_xa_freq_yak_stat <- p_xa_freq_yak_stat + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

	p_xa_size_ana_stat <- plot_neigh(autosome_burst_size_mean, x_burst_size_mean, "Autosomal", "X_linked", "ana", -11*6/75+1.15, 6, -2.5*6/75+1.05, -11*6/75+1.15, "Burst Size", "Burst Size in D. ananassae, X:A", TRUE, mtx_size_pval_xa[3,], mtx_size_N_xa[3,], mtx_size_Nref_xa[3,])
	p_xa_freq_ana_stat <- plot_neigh(autosome_burst_freq_mean, x_burst_freq_mean, "Autosomal", "X_linked", "ana", -11*6/75, 6, -2.5*6/75, -11*6/75, "Burst Frequency", "Burst Frequency in D. ananassae, X:A", TRUE, mtx_freq_pval_xa[3,], mtx_freq_N_xa[3,], mtx_freq_Nref_xa[3,])

	p_xa_size_ana_stat <- p_xa_size_ana_stat + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	p_xa_freq_ana_stat <- p_xa_freq_ana_stat + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


	grid.arrange(p_xa_size_mel_stat, p_xa_freq_mel_stat, p_xa_size_yak_stat, p_xa_freq_yak_stat, p_xa_size_ana_stat, p_xa_freq_ana_stat, ncol=2)

	dev.off()

	svg("sc_liftoff/bursting/comb_xa_all_stats.svg", width=27.3, height=12, pointsize=3)

	grid.arrange(p_xa_size_mel_stat, p_xa_freq_mel_stat, p_xa_size_yak_stat, p_xa_freq_yak_stat, p_xa_size_ana_stat, p_xa_freq_ana_stat, ncol=2)

	dev.off()



	svg("sc_liftoff/bursting/comb_xa_all_stats_ugg.svg", width=27.3, height=12, pointsize=3)

	p_xa_size_mel_stat <- plot_neigh1(autosome_burst_size_mean, x_burst_size_mean, "Autosomal", "X_linked", "mel", -11*32/75, 32, -2.5*32/75, -11*32/75, "Burst Size", "Burst Size in D. melanogaster, X:A", TRUE, mtx_size_pval_xa[1,], mtx_size_N_xa[1,], mtx_size_Nref_xa[1,])
	p_xa_freq_mel_stat <- plot_neigh1(autosome_burst_freq_mean, x_burst_freq_mean, "Autosomal", "X_linked", "mel", -11*7.2/75, 7.2, -2.5*7.2/75, -11*7.2/75, "Burst Frequency", "Burst Frequency in D. melanogaster, X:A", TRUE, mtx_freq_pval_xa[1,], mtx_freq_N_xa[1,], mtx_freq_Nref_xa[1,])

	p_xa_size_mel_stat <- p_xa_size_mel_stat + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	p_xa_freq_mel_stat <- p_xa_freq_mel_stat + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

	p_xa_size_yak_stat <- plot_neigh1(autosome_burst_size_mean, x_burst_size_mean, "Autosomal", "X_linked", "yak", -11*85/75, 85, -2.5*85/75, -11*85/75, "Burst Size", "Burst Size in D. yakuba, X:A", TRUE, mtx_size_pval_xa[2,], mtx_size_N_xa[2,], mtx_size_Nref_xa[2,])
	p_xa_freq_yak_stat <- plot_neigh1(autosome_burst_freq_mean, x_burst_freq_mean, "Autosomal", "X_linked", "yak", -11*2.5/75, 2.5, -2.5*2.5/75, -11*2.5/75, "Burst Frequency", "Burst Frequency in D. yakuba, X:A", TRUE, mtx_freq_pval_xa[2,], mtx_freq_N_xa[2,], mtx_freq_Nref_xa[2,])

	p_xa_size_yak_stat <- p_xa_size_yak_stat + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	p_xa_freq_yak_stat <- p_xa_freq_yak_stat + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

	p_xa_size_ana_stat <- plot_neigh1(autosome_burst_size_mean, x_burst_size_mean, "Autosomal", "X_linked", "ana", -11*6/75+1.15, 6, -2.5*6/75+1.05, -11*6/75+1.15, "Burst Size", "Burst Size in D. ananassae, X:A", TRUE, mtx_size_pval_xa[3,], mtx_size_N_xa[3,], mtx_size_Nref_xa[3,])
	p_xa_freq_ana_stat <- plot_neigh1(autosome_burst_freq_mean, x_burst_freq_mean, "Autosomal", "X_linked", "ana", -11*6/75, 6, -2.5*6/75, -11*6/75, "Burst Frequency", "Burst Frequency in D. ananassae, X:A", TRUE, mtx_freq_pval_xa[3,], mtx_freq_N_xa[3,], mtx_freq_Nref_xa[3,])

	p_xa_size_ana_stat <- p_xa_size_ana_stat + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	p_xa_freq_ana_stat <- p_xa_freq_ana_stat + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


	# png("sc_liftoff/bursting/comb_xa_all_stats_ugg1.png", width=27.3, height=12, units="in", res=300)

	grid.arrange(p_xa_size_mel_stat, p_xa_freq_mel_stat, p_xa_size_yak_stat, p_xa_freq_yak_stat, p_xa_size_ana_stat, p_xa_freq_ana_stat, ncol=2)



	dev.off()
}




gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

col_wheel <- gg_color_hue(8)

ct_cols <- col_wheel[c(7, 3, 6, 2, 5, 1, 4, 8)]

ct_cols1 <- ct_cols[1:7]			
ct_cols2 <- ct_cols[1:7]			
ct_cols3 <- ct_cols[1:7]			




src <- factor(c(rep("N total", 7), rep("N X-linked", 7), rep("N Autosomal", 7)), levels=c("N total", "N X-linked", "N Autosomal"))
src1 <- factor(c(rep("N total", 7), rep("N X-linked", 7), rep("N Autosomal", 7)), levels=c("N total", "N X-linked", "N Autosomal"))

nm <- c(mtx_size_N_ref[1,-8], mtx_size_N_xa[1,-8], mtx_size_Nref_xa[1,-8])
ny <- c(mtx_size_N_ref[2,-8], mtx_size_N_xa[2,-8], mtx_size_Nref_xa[2,-8])
ena <-  c(mtx_size_N_ref[3,-4], mtx_size_N_xa[1,-4], mtx_size_Nref_xa[3,-4])

ceetee <- rep(colnames(mtx_size_N)[-8], 3)
ceetee1 <- rep(colnames(mtx_size_N)[-4], 3)

ceetee[which(ceetee == "GSC_Early_spermatogonia")] <- "GSC/Early Spermatogonia"
ceetee[which(ceetee == "Late_spermatogonia")] <- "Late Spermatogonia"
ceetee[which(ceetee == "Early_spermatocyte")] <- "Early Spermatocyte"
ceetee[which(ceetee == "Late_spermatocyte")] <- "Late Spermatocyte"
ceetee[which(ceetee == "Early_spermatid")] <- "Early Spermatid"
ceetee[which(ceetee == "Late_spermatid")] <- "Late Spermatid"

ceetee1[which(ceetee1 == "GSC_Early_spermatogonia")] <- "GSC/Early Spermatogonia"
ceetee1[which(ceetee1 == "Late_spermatogonia")] <- "Late Spermatogonia"
ceetee1[which(ceetee1 == "Early_spermatocyte")] <- "Early Spermatocyte"
ceetee1[which(ceetee1 == "Late_spermatocyte")] <- "Late Spermatocyte"
ceetee1[which(ceetee1 == "ananassae_spermatid")] <- "ananassae Spermatocyte"
ceetee1[which(ceetee1 == "Early_spermatid")] <- "Early Spermatid"
ceetee1[which(ceetee1 == "Late_spermatid")] <- "Late Spermatid"

ceetee <- factor(ceetee, levels=c("Somatic", "GSC/Early Spermatogonia", "Late Spermatogonia", "Early Spermatocyte", "Late Spermatocyte", "Early Spermatid", "Late Spermatid"))
ceetee1 <- factor(ceetee1, levels=c("Somatic", "GSC/Early Spermatogonia", "Late Spermatogonia", "ananassae Spermatocyte", "Late Spermatocyte", "Early Spermatid", "Late Spermatid"))

df_N_mel <- data.frame(nm, src, ceetee)
df_N_yak <- data.frame(ny, src, ceetee)
df_N_ana <- data.frame(ena, src1, ceetee1)

if(plot_new){
	png("sc_liftoff/bursting/N_boxes.png", width=2400*4, height=3600*1.5, res=600, pointsize=3)

	p1 <- ggplot(df_N_mel, aes(x=ceetee, y=nm, fill=src), outlier.shape = NA) + geom_bar(position="dodge", stat="identity", color="black") + ylab("N") + xlab("Cell Type") + ggtitle("melanogaster")
	p1 <- p1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

	p2 <- ggplot(df_N_yak, aes(x=ceetee, y=ny, fill=src), outlier.shape = NA) + geom_bar(position="dodge", stat="identity", color="black") + ylab("N") + xlab("Cell Type") + ggtitle("yakuba")
	p2 <- p2 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

	p3 <- ggplot(df_N_ana, aes(x=ceetee1, y=ena, fill=src1), outlier.shape = NA) + geom_bar(position="dodge", stat="identity", color="black") + ylab("N") + xlab("Cell Type") + ggtitle("ananassae")
	p3 <- p3 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

	grid.arrange(p1, p2, p3, ncol=1)

	dev.off()

	svg("sc_liftoff/bursting/N_boxes.svg", width=16, height=9, pointsize=3)

	grid.arrange(p1, p2, p3, ncol=1)

	dev.off()
}

dnt_dat_plus <- read.csv("sc_liftoff/dnt_dat_plus.csv")

all_burst_size_mel <- all_burst_size_mean[all_burst_size_mean[,3] == "mel",]
all_burst_size_yak <- all_burst_size_mean[all_burst_size_mean[,3] == "yak",]
all_burst_size_ana <- all_burst_size_mean[all_burst_size_mean[,3] == "ana",]

all_burst_freq_mel <- all_burst_freq_mean[all_burst_freq_mean[,3] == "mel",]
all_burst_freq_yak <- all_burst_freq_mean[all_burst_freq_mean[,3] == "yak",]
all_burst_freq_ana <- all_burst_freq_mean[all_burst_freq_mean[,3] == "ana",]

png("sc_liftoff/test.png")
hist(unlist(all_burst_size_mel[,4:11]))
dev.off()


library(actuar)
library(fitdistr)



pareto_mtx <- function(x, shape, scale){
	return(x)
	out <- matrix(0, nrow=nrow(x), ncol=ncol(x))
	
	for(i in 1:ncol(out)){
		out[,i] <- ppareto(x[,i], shape, scale)
	}
	
	return(out)
}

mtx_uniques <- function(in_mat, genes){
	out <- in_mat
	out <- out[match(genes, in_mat[,1]),]
	out <- out[!is.na(out[,1]),]
		
	return(out)

}

common_genes <- intersect(all_burst_size_mel[,1], intersect(all_burst_size_yak[,1], all_burst_size_ana[,1]))
common_genes <- intersect(common_genes, mel_onlygenes)

burst_size_mel_pareto <- mtx_uniques(all_burst_size_mel, common_genes)
dist_mel_size <- unlist(burst_size_mel_pareto[,4:11])
dist_mel_size <- dist_mel_size[dist_mel_size >= 1]
fit.gamma_mel_size <- fitdist(dist_mel_size, distr="pareto", method="mle")
burst_size_mel_pareto[, 4:11] <- pareto_mtx(burst_size_mel_pareto[,4:11], fit.gamma_mel_size$estimate["shape"], fit.gamma_mel_size$estimate["scale"])

burst_size_yak_pareto <- mtx_uniques(all_burst_size_yak, common_genes)
dist_yak_size <- unlist(burst_size_yak_pareto[,4:11])
dist_yak_size <- dist_yak_size[dist_yak_size >= 1]
fit.gamma_yak_size <- fitdist(dist_yak_size, distr="pareto", method="mle")
burst_size_yak_pareto <- mtx_uniques(all_burst_size_yak, common_genes)
burst_size_yak_pareto[, 4:11] <- pareto_mtx(burst_size_yak_pareto[,4:11], fit.gamma_yak_size$estimate["shape"], fit.gamma_yak_size$estimate["scale"])

burst_size_ana_pareto <- mtx_uniques(all_burst_size_ana, common_genes)
dist_ana_size <- unlist(burst_size_ana_pareto[,4:11])
dist_ana_size <- dist_ana_size[dist_ana_size >= 1]
fit.gamma_ana_size <- fitdist(dist_ana_size, distr="pareto", method="mle")
burst_size_ana_pareto <- mtx_uniques(all_burst_size_ana, common_genes)
burst_size_ana_pareto[, 4:11] <- pareto_mtx(burst_size_ana_pareto[,4:11], fit.gamma_ana_size$estimate["shape"], fit.gamma_ana_size$estimate["scale"])

burst_freq_mel_pareto <- mtx_uniques(all_burst_freq_mel, common_genes)
dist_mel_freq <- unlist(burst_freq_mel_pareto[,4:11])
dist_mel_freq <- dist_mel_freq[dist_mel_freq >= 1]
fit.gamma_mel_freq <- fitdist(dist_mel_freq, distr="pareto", method="mle")
burst_freq_mel_pareto[, 4:11] <- pareto_mtx(burst_freq_mel_pareto[,4:11], fit.gamma_mel_freq$estimate["shape"], fit.gamma_mel_freq$estimate["scale"])

burst_freq_yak_pareto <- mtx_uniques(all_burst_freq_yak, common_genes)
dist_yak_freq <- unlist(burst_freq_yak_pareto[,4:11])
dist_yak_freq <- dist_yak_freq[dist_yak_freq >= 1]
fit.gamma_yak_freq <- fitdist(dist_yak_freq, distr="pareto", method="mle")
burst_freq_yak_pareto <- mtx_uniques(all_burst_freq_yak, common_genes)
burst_freq_yak_pareto[, 4:11] <- pareto_mtx(burst_freq_yak_pareto[,4:11], fit.gamma_yak_freq$estimate["shape"], fit.gamma_yak_freq$estimate["scale"])

burst_freq_ana_pareto <- mtx_uniques(all_burst_freq_ana, common_genes)
dist_ana_freq <- unlist(burst_freq_ana_pareto[,4:11])
dist_ana_freq <- dist_ana_freq[dist_ana_freq >= 1]
fit.gamma_ana_freq <- fitdist(dist_ana_freq, distr="pareto", method="mle")
burst_freq_ana_pareto <- mtx_uniques(all_burst_freq_ana, common_genes)
burst_freq_ana_pareto[, 4:11] <- pareto_mtx(burst_freq_ana_pareto[,4:11], fit.gamma_ana_freq$estimate["shape"], fit.gamma_ana_freq$estimate["scale"])


png("sc_liftoff/pareto_mel_size.png", width=2400, height=2400, res=600, pointsize=3)
plot(fit.gamma_mel_size)
dev.off()

png("sc_liftoff/pareto_yak_size.png", width=2400, height=2400, res=600, pointsize=3)
plot(fit.gamma_yak_size)
dev.off()

png("sc_liftoff/pareto_ana_size.png", width=2400, height=2400, res=600, pointsize=3)
plot(fit.gamma_ana_size)
dev.off()

png("sc_liftoff/pareto_mel_freq.png", width=2400, height=2400, res=600, pointsize=3)
plot(fit.gamma_mel_freq)
dev.off()

png("sc_liftoff/pareto_yak_freq.png", width=2400, height=2400, res=600, pointsize=3)
plot(fit.gamma_yak_freq)
dev.off()

png("sc_liftoff/pareto_ana_freq.png", width=2400, height=2400, res=600, pointsize=3)
plot(fit.gamma_ana_freq)
dev.off()



delta_dist <- function(x, y){
	out <- x

	for(i in 1:ncol(x)){
		out[,i] <- y[,i] - x[,i]
		
		idx <- unlist(apply(cbind(y[,i], x[,i]), 1, function(m) return(sum(m == 0) > 0)))
		out[idx,i] <- 0
	}
	
	return(out)
}

burst_size_diff_meya <- burst_size_mel_pareto
burst_size_diff_mean <- burst_size_mel_pareto
burst_freq_diff_meya <- burst_freq_mel_pareto
burst_freq_diff_mean <- burst_freq_mel_pareto

burst_size_diff_meya[,4:11] <- delta_dist(burst_size_yak_pareto[,4:11], burst_size_mel_pareto[,4:11]) 
burst_size_diff_mean[,4:11] <- delta_dist(burst_size_ana_pareto[,4:11], burst_size_mel_pareto[,4:11])
 
burst_freq_diff_meya[,4:11] <- delta_dist(burst_freq_yak_pareto[,4:11], burst_freq_mel_pareto[,4:11]) 
burst_freq_diff_mean[,4:11] <- delta_dist(burst_freq_ana_pareto[,4:11], burst_freq_mel_pareto[,4:11])


obs_genes <- intersect(which(rowSums(burst_size_diff_meya[,4:11]) > 0), which(rowSums(burst_size_diff_mean[,4:11]) > 0))
all_dnt <- unique(unlist(dnt_dat_plus[,c("nei_L2_sym_mel","nei_L1_sym_mel", "nei_R1_sym_mel", "nei_R2_sym_mel")]))
obs_dnt <- intersect(burst_size_diff_mean[obs_genes,1], all_dnt)

#621 all neighbors, down to 36 that have observations in all 3 spec.

dnt_size_diff_mean <- burst_size_diff_mean[match(obs_dnt, burst_size_diff_mean[,1]),]

all_size_diff_mean_noz_list <- list(burst_size_diff_mean[,c(1,4)], burst_size_diff_mean[,c(1,5)], burst_size_diff_mean[,c(1,6)], burst_size_diff_mean[,c(1,7)], burst_size_diff_mean[,c(1,8)], burst_size_diff_mean[,c(1,9)], burst_size_diff_mean[,c(1,10)], burst_size_diff_mean[,c(1,11)])

all_size_diff_mean_noz_list <- lapply(all_size_diff_mean_noz_list, function(x) x[which(!x[,2] == 0),])
names(all_size_diff_mean_noz_list) <- colnames(burst_freq_diff_mean)[4:11]

all_size_diff_mean_noz_all <- unlist(lapply(all_size_diff_mean_noz_list, function(x) x[,2]))

dnt_size_diff_mean_noz_list <- list(dnt_size_diff_mean[,c(1,4)], dnt_size_diff_mean[,c(1,5)], dnt_size_diff_mean[,c(1,6)], dnt_size_diff_mean[,c(1,7)], dnt_size_diff_mean[,c(1,8)], dnt_size_diff_mean[,c(1,9)], dnt_size_diff_mean[,c(1,10)], dnt_size_diff_mean[,c(1,11)])
 
dnt_size_diff_mean_noz_list <- lapply(dnt_size_diff_mean_noz_list, function(x) x[which(!x[,2] == 0),])
names(dnt_size_diff_mean_noz_list) <- colnames(burst_freq_diff_mean)[4:11]

dnt_size_diff_mean_noz_all <- unlist(lapply(all_size_diff_mean_noz_list, function(x) x[,2]))

n_ana_all <- length(unique(unlist(lapply(all_size_diff_mean_noz_list, function(x) x[,1]))))
n_ana_dnt <- length(unique(unlist(lapply(dnt_size_diff_mean_noz_list, function(x) x[,1]))))

for(i in (1:8)[-c(4, 7)]) print(c(names(dnt_size_diff_mean_noz_list)[i], t.test(dnt_size_diff_mean_noz_list[[i]][,2], all_size_diff_mean_noz_list[[i]][, 2])$p.value))

png("sc_liftoff/test.png", width=2400, height=2400, res=600, pointsize=3)
# hist(dnt_size_diff_mean_noz_list[[2]][,2], freq=F, col=rgb(1, 0, 0, 1/4), xlim=c(-1,1), breaks=seq(-1, 1, length.out=40))
# hist(all_size_diff_mean_noz_list[[2]][,2], freq=F, col=rgb(0, 0, 1, 1/4), xlim=c(-1,1), breaks=seq(-1, 1, length.out=40), add=T)
hist(all_size_diff_mean_noz_all, freq=F, breaks=100)
dev.off()

png("sc_liftoff/all_size_diff_mean.png", width=2400, height=2400, res=800, pointsize=3)
# hist(dnt_size_diff_mean_noz_list[[2]][,2], freq=F, col=rgb(1, 0, 0, 1/4), xlim=c(-1,1), breaks=seq(-1, 1, length.out=40))
# hist(all_size_diff_mean_noz_list[[2]][,2], freq=F, col=rgb(0, 0, 1, 1/4), xlim=c(-1,1), breaks=seq(-1, 1, length.out=40), add=T)
hist(all_size_diff_mean_noz_all, freq=F, breaks=seq(min(all_size_diff_mean_noz_all), max(all_size_diff_mean_noz_all), length.out=400), main="Difference in Burst Size (mel-ana)", xlab="Parameter Difference (mel-ana)", xlim=c(-25, 100), cex.lab=1.5, cex.main=1.5)
hist(all_size_diff_mean_noz_all, freq=F, breaks=seq(min(all_size_diff_mean_noz_all), max(all_size_diff_mean_noz_all), length.out=400), main="Difference in Burst Size (mel-ana)", xlab="Parameter Difference (mel-ana)", xlim=c(-25, 100), cex.lab=1.5, cex.main=1.5)
abline(v=sort(all_size_diff_mean_noz_all)[floor(length(all_size_diff_mean_noz_all)*0.05)], lwd=0.5, lty=2)
abline(v=sort(all_size_diff_mean_noz_all)[floor(length(all_size_diff_mean_noz_all)*0.95)], lwd=0.5, lty=2)
dev.off()




dnt_size_diff_meya <- burst_size_diff_meya[match(obs_dnt, burst_size_diff_meya[,1]),]

all_size_diff_meya_noz_list <- list(burst_size_diff_meya[,c(1,4)], burst_size_diff_meya[,c(1,5)], burst_size_diff_meya[,c(1,6)], burst_size_diff_meya[,c(1,7)], burst_size_diff_meya[,c(1,8)], burst_size_diff_meya[,c(1,9)], burst_size_diff_meya[,c(1,10)], burst_size_diff_meya[,c(1,11)])

all_size_diff_meya_noz_list <- lapply(all_size_diff_meya_noz_list, function(x) x[which(!x[,2] == 0),])
names(all_size_diff_meya_noz_list) <- colnames(burst_freq_diff_meya)[4:11]

all_size_diff_meya_noz_all <- unlist(lapply(all_size_diff_meya_noz_list, function(x) x[,2]))

dnt_size_diff_meya_noz_list <- list(dnt_size_diff_meya[,c(1,4)], dnt_size_diff_meya[,c(1,5)], dnt_size_diff_meya[,c(1,6)], dnt_size_diff_meya[,c(1,7)], dnt_size_diff_meya[,c(1,8)], dnt_size_diff_meya[,c(1,9)], dnt_size_diff_meya[,c(1,10)], dnt_size_diff_meya[,c(1,11)])
 
dnt_size_diff_meya_noz_list <- lapply(dnt_size_diff_meya_noz_list, function(x) x[which(!x[,2] == 0),])
names(dnt_size_diff_meya_noz_list) <- colnames(burst_freq_diff_meya)[4:11]

n_yak_all <- length(unique(unlist(lapply(all_size_diff_meya_noz_list, function(x) x[,1]))))
n_yak_dnt <- length(unique(unlist(lapply(dnt_size_diff_meya_noz_list, function(x) x[,1]))))

for(i in (1:8)) print(c(names(dnt_size_diff_meya_noz_list)[i], t.test(dnt_size_diff_meya_noz_list[[i]][,2], all_size_diff_meya_noz_list[[i]][, 2])$p.value))

png("sc_liftoff/test1.png", width=2400, height=2400, res=600, pointsize=3)
hist(dnt_size_diff_meya_noz_list[[3]][,2], freq=F, col=rgb(1, 0, 0, 1/4), xlim=c(-1,1), breaks=seq(-1, 1, length.out=40))
hist(all_size_diff_meya_noz_list[[3]][,2], freq=F, col=rgb(0, 0, 1, 1/4), xlim=c(-1,1), breaks=seq(-1, 1, length.out=40), add=T)
dev.off()

png("sc_liftoff/all_size_diff_meya.png", width=2400, height=2400, res=800, pointsize=3)
# hist(dnt_size_diff_meya_noz_list[[2]][,2], freq=F, col=rgb(1, 0, 0, 1/4), xlim=c(-1,1), breaks=seq(-1, 1, length.out=40))
# hist(all_size_diff_meya_noz_list[[2]][,2], freq=F, col=rgb(0, 0, 1, 1/4), xlim=c(-1,1), breaks=seq(-1, 1, length.out=40), add=T)
hist(all_size_diff_meya_noz_all, freq=F, seq(min(all_size_diff_meya_noz_all), max(all_size_diff_meya_noz_all), length.out=400), main="Difference in Burst Size (mel-yak)", xlab="Parameter Difference (mel-yak)", xlim=c(-175, 50), cex.lab=1.5, cex.main=1.5)
hist(all_size_diff_meya_noz_all, freq=F, seq(min(all_size_diff_meya_noz_all), max(all_size_diff_meya_noz_all), length.out=400), main="Difference in Burst Size (mel-yak)", xlab="Parameter Difference (mel-yak)", xlim=c(-175, 50), cex.lab=1.5, cex.main=1.5)
abline(v=sort(all_size_diff_meya_noz_all)[floor(length(all_size_diff_meya_noz_all)*0.05)], lwd=0.5, lty=2)
abline(v=sort(all_size_diff_meya_noz_all)[floor(length(all_size_diff_meya_noz_all)*0.95)], lwd=0.5, lty=2)
dev.off()



dist_list <- lapply(list(unname(dnt_dat_plus[,c(6, 10, 11,1,2)]), unname(dnt_dat_plus[,c(12, 16, 17,1,2)]), unname(dnt_dat_plus[,c(18, 22, 23, 1,2)]), unname(dnt_dat_plus[,c(24, 28, 29, 1,2)])), function(x) return(x[x[,2] > 0,]))
 
dist_list_mtx <- data.frame(genes=c(dist_list[[1]][,1], dist_list[[2]][,1], dist_list[[3]][,1],dist_list[[4]][,1]), dist=as.numeric(c(dist_list[[1]][,2], dist_list[[2]][,2], dist_list[[3]][,2], dist_list[[4]][,2])), syn=c(dist_list[[1]][,3], dist_list[[2]][,3], dist_list[[3]][,3],dist_list[[4]][,3]),origin=c(dist_list[[1]][,4], dist_list[[2]][,4], dist_list[[3]][,4],dist_list[[4]][,4]), dnt=c(dist_list[[1]][,5], dist_list[[2]][,5], dist_list[[3]][,5],dist_list[[4]][,5]))

dnt_size_diff_mean_noz_list_indist <- lapply(dnt_size_diff_mean_noz_list, function(x){
	if(length(x) <= 1) return(x)
	return(x[which(x[,1] %in% dist_list_mtx[,1]),])})

dnt_size_diff_mean_noz_list_alldist <- lapply(dnt_size_diff_mean_noz_list_indist, function(x) unlist(lapply(x[,1], function(y) min(dist_list_mtx[which(dist_list_mtx[,1] == y),2]))))

dnt_size_diff_mean_noz_list_syndist <- lapply(dnt_size_diff_mean_noz_list_indist, function(x) unlist(lapply(x[,1], function(y) dist_list_mtx[which(dist_list_mtx[,1] == y),3])))

dnt_size_diff_mean_noz_list_oridist <- lapply(dnt_size_diff_mean_noz_list_indist, function(x) unlist(lapply(x[,1], function(y) dist_list_mtx[which(dist_list_mtx[,1] == y),4])))

dnt_size_diff_mean_noz_list_dntdist <- lapply(dnt_size_diff_mean_noz_list_indist, function(x) unlist(lapply(x[,1], function(y) dist_list_mtx[which(dist_list_mtx[,1] == y),5])))

dnt_size_diff_mean_combined <- data.frame(deltap=unlist(lapply(dnt_size_diff_mean_noz_list_indist, function(x) return(x[,2]))), dist=unlist(dnt_size_diff_mean_noz_list_alldist), syn=unlist(dnt_size_diff_mean_noz_list_syndist), origin=unlist(dnt_size_diff_mean_noz_list_oridist), dnt=unlist(dnt_size_diff_mean_noz_list_dntdist))

idx_inc <- dnt_size_diff_mean_combined[,3] %in% c("tdn", "div")

temp2 <- dnt_size_diff_mean_combined[,2] < 40000
temp2 <- intersect(which(temp2), which(dnt_size_diff_mean_combined[,1] < 200))

png("sc_liftoff/dnt_size_diff_dist_mean.png", width=2400, height=2400, res=800, pointsize=4.5)
plot(dnt_size_diff_mean_combined[temp2,2], dnt_size_diff_mean_combined[temp2,1], xlab="Distance to de novo transcript (in mel)", ylab="Difference in Burst Size (mel-ana)", main="Genes neighboring de novo transcripts (ananassae)", cex.lab=1.5, cex.main=1.5, ylim=c(-25, 100))
abline(h=sort(all_size_diff_mean_noz_all)[floor(length(all_size_diff_mean_noz_all)*0.05)], lwd=0.5, lty=2)
abline(h=sort(all_size_diff_mean_noz_all)[floor(length(all_size_diff_mean_noz_all)*0.95)], lwd=0.5, lty=2)
dev.off()


dnt_size_diff_meya_noz_list_indist <- lapply(dnt_size_diff_meya_noz_list, function(x){
	if(length(x) <= 1) return(x)
	return(x[which(x[,1] %in% dist_list_mtx[,1]),])})

dnt_size_diff_meya_noz_list_alldist <- lapply(dnt_size_diff_meya_noz_list_indist, function(x) unlist(lapply(x[,1], function(y) min(dist_list_mtx[which(dist_list_mtx[,1] == y),2]))))

dnt_size_diff_meya_noz_list_syndist <- lapply(dnt_size_diff_meya_noz_list_indist, function(x) unlist(lapply(x[,1], function(y) dist_list_mtx[which(dist_list_mtx[,1] == y),3])))

dnt_size_diff_meya_noz_list_oridist <- lapply(dnt_size_diff_meya_noz_list_indist, function(x) unlist(lapply(x[,1], function(y) dist_list_mtx[which(dist_list_mtx[,1] == y),4])))

dnt_size_diff_meya_noz_list_dntdist <- lapply(dnt_size_diff_meya_noz_list_indist, function(x) unlist(lapply(x[,1], function(y) dist_list_mtx[which(dist_list_mtx[,1] == y),5])))

dnt_size_diff_meya_combined <- data.frame(deltap=unlist(lapply(dnt_size_diff_meya_noz_list_indist, function(x) return(x[,2]))), dist=unlist(dnt_size_diff_meya_noz_list_alldist), syn=unlist(dnt_size_diff_meya_noz_list_syndist), origin=unlist(dnt_size_diff_meya_noz_list_oridist), dnt=unlist(dnt_size_diff_meya_noz_list_dntdist))

idx_inc <- dnt_size_diff_meya_combined[,3] %in% c("tdn", "div")

temp3 <- dnt_size_diff_meya_combined[,2] < 40000
temp3 <- intersect(which(temp3), which(dnt_size_diff_meya_combined[,1] < 200))
temp3 <- intersect(temp3, which(dnt_size_diff_meya_combined[,4] %in% c("000001", "000111")))

png("sc_liftoff/dnt_size_diff_dist_meya.png", width=2400, height=2400, res=800, pointsize=4.5)
plot(dnt_size_diff_meya_combined[temp3,2], dnt_size_diff_meya_combined[temp3,1], xlab="Distance to de novo transcript (in mel)", ylab="Difference in Burst Size (mel-yak)", main="Genes neighboring de novo transcripts (yakuba)", cex.lab=1.5, cex.main=1.5, ylim=c(-175, 50))
abline(h=sort(all_size_diff_meya_noz_all)[floor(length(all_size_diff_meya_noz_all)*0.05)], lwd=0.5, lty=2)
abline(h=sort(all_size_diff_meya_noz_all)[floor(length(all_size_diff_meya_noz_all)*0.95)], lwd=0.5, lty=2)
dev.off()



temp <- dnt_size_diff_meya_combined[,4] %in% c("000001", "000111")
temp1 <- intersect(temp, temp3)






temp1 <- unlist(lapply(temp5[,1], function(x) min(temp[which(temp[,1] == x),2])))
 
temp2 <- which(!temp5[,"Early_spermatocyte"] == 0 )
  
png("sc_liftoff/test1.png", width=2400, height=2400, res=600, pointsize=3)
plot(temp1[temp2], temp5[temp2,"Early_spermatocyte"])
dev.off()

# save(list=ls(), file="sc_liftoff/bursting/wrkspce.rdata"))

# svgo -f raw_svg/ -o raw_svg_svgo
