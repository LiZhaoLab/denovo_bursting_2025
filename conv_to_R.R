#srun --pty -p zhao apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis/sc_liftoff/bursting/:/home --home /home --env "LD_LIBRARY_PATH=/root/miniconda/lib" /lustre/fs4/zhao_lab/scratch/ulee/docka/samap_fig.sif R


library(reticulate)
use_python("/root/miniconda/bin/python")
use_condaenv("/root/miniconda")

source_python("pickle_reader.py")
pickle_data <- read_pickle_file("ana_ananassae_spermatid_ML.pkl")

flatten_pickle <- function(pick_dat){
	df_pick <- as.data.frame(do.call("rbind", pick_dat[,1]))
	colnames(df_pick) <- c("k_on", "k_off", "k_syn")
	df_pick$burst_freq <- df_pick$k_on
	df_pick$burst_size <- df_pick$k_syn / df_pick$k_off
	df_pick$keep <- unlist(pick_dat[,2])
	
	rownames(df_pick) <- rownames(pick_dat)
	
	return(df_pick)
}

import_pickle <- function(pick_string){
	pickle_name <- paste(pick_string, ".pkl", sep="")

	this_pickle <- read_pickle_file(pickle_name)
	this_pickle <- flatten_pickle(this_pickle)

	return(this_pickle)
}

all_picks <- list.files()
all_picks <- all_picks[grepl(".pkl", all_picks)]

all_picks_trim <- unlist(lapply(all_picks, function(x) substr(x, 1, nchar(x)-4)))

all_picks_list <- lapply(all_picks_trim, function(x) import_pickle(x))
names(all_picks_list) <- all_picks_trim

save(list=c("all_picks_list"), file="allpickslist.rdata")