#! /bin/bash

cd /ru-auth/local/home/ulee/scratch_store/analysis/sc_liftoff/bursting

for i in ana_*.csv
do
	srun -p zhao apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis/sc_liftoff/bursting/:/home --home /home --env "LD_LIBRARY_PATH=/root/miniconda/lib" /lustre/fs4/zhao_lab/scratch/ulee/docka/samap_fig.sif python txburstML.py $i
done

for i in mel_*.csv
do
	srun -p zhao apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis/sc_liftoff/bursting/:/home --home /home --env "LD_LIBRARY_PATH=/root/miniconda/lib" /lustre/fs4/zhao_lab/scratch/ulee/docka/samap_fig.sif python txburstML.py $i
done

for i in yak_*.csv
do
	srun -p zhao apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis/sc_liftoff/bursting/:/home --home /home --env "LD_LIBRARY_PATH=/root/miniconda/lib" /lustre/fs4/zhao_lab/scratch/ulee/docka/samap_fig.sif python txburstML.py $i
done
