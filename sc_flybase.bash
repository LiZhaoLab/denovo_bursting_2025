#!/usr/bin/bash

cd /ru-auth/local/home/ulee/scratch_store/analysis/sc_liftoff/

#make 10x ref genomes
srun -p zhao apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis/sc_liftoff/:/home --home /home/ /lustre/fs4/zhao_lab/scratch/ulee/docka/samap.sif cellranger mkref \
--genome=dmel_fb --fasta=/home/genomes_fb/dmel-all-chromosome-r6.44.fasta --genes=/home/genomes_fb/dmel-all-r6.44.gtf

srun -p zhao apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis/sc_liftoff/:/home --home /home/ /lustre/fs4/zhao_lab/scratch/ulee/docka/samap.sif cellranger mkref \
--genome=dyak_fb --fasta=/home/genomes_fb/dyak-all-chromosome-r1.05.fasta --genes=/home/genomes_fb/dyak-all-r1.05.gtf

srun -p zhao apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis/sc_liftoff/:/home --home /home/ /lustre/fs4/zhao_lab/scratch/ulee/docka/samap.sif cellranger mkref \
--genome=dana_fb --fasta=/home/genomes_fb/dana-all-chromosome-r1.06.fasta --genes=/home/genomes_fb/dana-all-r1.06.gtf

#run 10x count
srun -p zhao apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis/:/home --home /home/ /lustre/fs4/zhao_lab/scratch/ulee/docka/samap.sif cellranger count \
--id=dmel_testis_fb --transcriptome=/home/sc_liftoff/dmel_fb --fastqs=/home/samap/dat/raw_mel --sample=517_novogene --force-cells=5000

mv /ru-auth/local/home/ulee/scratch_store/analysis/dmel_testis_fb /ru-auth/local/home/ulee/scratch_store/analysis/sc_liftoff/dmel_testis_fb

srun -p zhao apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis/:/home --home /home/ /lustre/fs4/zhao_lab/scratch/ulee/docka/samap.sif cellranger count \
--id=dyak_testis_fb --transcriptome=/home/sc_liftoff/dyak_fb --fastqs=/home/samap/dat/raw_yak --sample=A01_Cy01 --force-cells=5000

mv /ru-auth/local/home/ulee/scratch_store/analysis/dyak_testis_fb /ru-auth/local/home/ulee/scratch_store/analysis/sc_liftoff/dyak_testis_fb

srun -p zhao apptainer exec --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/analysis/:/home --home /home/ /lustre/fs4/zhao_lab/scratch/ulee/docka/samap.sif cellranger count \
--id=dana_testis_fb --transcriptome=/home/sc_liftoff/dana_fb --fastqs=/home/samap/dat/raw_ana --sample=Dana_CKDL230014380-1A_H75FVDSX7 --force-cells=5000

mv /ru-auth/local/home/ulee/scratch_store/analysis/dana_testis_fb /ru-auth/local/home/ulee/scratch_store/analysis/sc_liftoff/dana_testis_fb
