#!/usr/bin/bash

echo "transcript" > /ru-auth/local/home/ulee/scratch_store/analysis/liftofftools_features.txt
echo "CDS" >> /ru-auth/local/home/ulee/scratch_store/analysis/liftofftools_features.txt

#sim
apptainer exec --no-home --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/:/home /ru-auth/local/home/ulee/scratch_store/docka/liftofftools.sif \
liftofftools -r /home/genomes/dmel-all-chromosome-r6.15.fasta -t /home/genomes/dsim-all-chromosome-r2.02.fasta \
-rg /home/analysis/comb_gr_nochr.gtf -tg /home/analysis/comb_gr_nochr_liftoff_sim.gff_polished \
-f /home/analysis/liftofftools_features.txt -infer-genes -dir /home/analysis/lotools/sim_lo all \

apptainer exec --no-home --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/:/home /ru-auth/local/home/ulee/scratch_store/docka/liftofftools.sif \
liftofftools -r /home/genomes/dmel-all-chromosome-r6.15.fasta -t /home/genomes/dsim-all-chromosome-r2.02.fasta \
-rg /home/analysis/comb_gr_nochr.gtf -tg /home/analysis/comb_gr_nochr_halp_sim.gtf \
-f /home/analysis/liftofftools_features.txt -infer-genes -dir /home/analysis/lotools/sim_halp all \

apptainer exec --no-home --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/:/home /ru-auth/local/home/ulee/scratch_store/docka/liftofftools.sif \
liftofftools -r /home/genomes/dmel-all-chromosome-r6.15.fasta -t /home/genomes/dsim-all-chromosome-r2.02.fasta \
-rg /home/analysis/comb_gr_nochr.gtf -tg /home/analysis/comb_gr_nochr_sort_sim.gtf \
-f /home/analysis/liftofftools_features.txt -infer-genes -dir /home/analysis/lotools/sim_hal_sort all \

#sec
apptainer exec --no-home --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/:/home /ru-auth/local/home/ulee/scratch_store/docka/liftofftools.sif \
liftofftools -r /home/genomes/dmel-all-chromosome-r6.15.fasta -t /home/genomes/dsec-all-chromosome-r1.3.fasta \
-rg /home/analysis/comb_gr_nochr.gtf -tg /home/analysis/comb_gr_nochr_liftoff_sec.gff_polished \
-f /home/analysis/liftofftools_features.txt -infer-genes -dir /home/analysis/lotools/sec_lo all \

apptainer exec --no-home --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/:/home /ru-auth/local/home/ulee/scratch_store/docka/liftofftools.sif \
liftofftools -r /home/genomes/dmel-all-chromosome-r6.15.fasta -t /home/genomes/dsec-all-chromosome-r1.3.fasta \
-rg /home/analysis/comb_gr_nochr.gtf -tg /home/analysis/comb_gr_nochr_halp_sec.gtf \
-f /home/analysis/liftofftools_features.txt -infer-genes -dir /home/analysis/lotools/sec_halp all \

apptainer exec --no-home --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/:/home /ru-auth/local/home/ulee/scratch_store/docka/liftofftools.sif \
liftofftools -r /home/genomes/dmel-all-chromosome-r6.15.fasta -t /home/genomes/dsec-all-chromosome-r1.3.fasta \
-rg /home/analysis/comb_gr_nochr.gtf -tg /home/analysis/comb_gr_nochr_sort_sec.gtf \
-f /home/analysis/liftofftools_features.txt -infer-genes -dir /home/analysis/lotools/sec_hal_sort all \

#yak
apptainer exec --no-home --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/:/home /ru-auth/local/home/ulee/scratch_store/docka/liftofftools.sif \
liftofftools -r /home/genomes/dmel-all-chromosome-r6.15.fasta -t /home/genomes/dyak-all-chromosome-r1.3.fasta \
-rg /home/analysis/comb_gr_nochr.gtf -tg /home/analysis/comb_gr_nochr_liftoff_yak.gff_polished \
-f /home/analysis/liftofftools_features.txt -infer-genes -dir /home/analysis/lotools/yak_lo all \

apptainer exec --no-home --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/:/home /ru-auth/local/home/ulee/scratch_store/docka/liftofftools.sif \
liftofftools -r /home/genomes/dmel-all-chromosome-r6.15.fasta -t /home/genomes/dyak-all-chromosome-r1.3.fasta \
-rg /home/analysis/comb_gr_nochr.gtf -tg /home/analysis/comb_gr_nochr_halp_yak.gtf \
-f /home/analysis/liftofftools_features.txt -infer-genes -dir /home/analysis/lotools/yak_halp all \

apptainer exec --no-home --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/:/home /ru-auth/local/home/ulee/scratch_store/docka/liftofftools.sif \
liftofftools -r /home/genomes/dmel-all-chromosome-r6.15.fasta -t /home/genomes/dyak-all-chromosome-r1.3.fasta \
-rg /home/analysis/comb_gr_nochr.gtf -tg /home/analysis/comb_gr_nochr_sort_yak.gtf \
-f /home/analysis/liftofftools_features.txt -infer-genes -dir /home/analysis/lotools/yak_hal_sort all \

#ere
apptainer exec --no-home --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/:/home /ru-auth/local/home/ulee/scratch_store/docka/liftofftools.sif \
liftofftools -r /home/genomes/dmel-all-chromosome-r6.15.fasta -t /home/genomes/dere-all-chromosome-r1.3.fasta \
-rg /home/analysis/comb_gr_nochr.gtf -tg /home/analysis/comb_gr_nochr_liftoff_ere.gff_polished \
-f /home/analysis/liftofftools_features.txt -infer-genes -dir /home/analysis/lotools/ere_lo all \

apptainer exec --no-home --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/:/home /ru-auth/local/home/ulee/scratch_store/docka/liftofftools.sif \
liftofftools -r /home/genomes/dmel-all-chromosome-r6.15.fasta -t /home/genomes/dere-all-chromosome-r1.3.fasta \
-rg /home/analysis/comb_gr_nochr.gtf -tg /home/analysis/comb_gr_nochr_halp_ere.gtf \
-f /home/analysis/liftofftools_features.txt -infer-genes -dir /home/analysis/lotools/ere_halp all \

apptainer exec --no-home --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/:/home /ru-auth/local/home/ulee/scratch_store/docka/liftofftools.sif \
liftofftools -r /home/genomes/dmel-all-chromosome-r6.15.fasta -t /home/genomes/dere-all-chromosome-r1.3.fasta \
-rg /home/analysis/comb_gr_nochr.gtf -tg /home/analysis/comb_gr_nochr_sort_ere.gtf \
-f /home/analysis/liftofftools_features.txt -infer-genes -dir /home/analysis/lotools/ere_hal_sort all \

#ana
apptainer exec --no-home --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/:/home /ru-auth/local/home/ulee/scratch_store/docka/liftofftools.sif \
liftofftools -r /home/genomes/dmel-all-chromosome-r6.15.fasta -t /home/genomes/dana-all-chromosome-r1.3.fasta \
-rg /home/analysis/comb_gr_nochr.gtf -tg /home/analysis/comb_gr_nochr_liftoff_ana.gff_polished \
-f /home/analysis/liftofftools_features.txt -infer-genes -dir /home/analysis/lotools/ana_lo all \

apptainer exec --no-home --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/:/home /ru-auth/local/home/ulee/scratch_store/docka/liftofftools.sif \
liftofftools -r /home/genomes/dmel-all-chromosome-r6.15.fasta -t /home/genomes/dana-all-chromosome-r1.3.fasta \
-rg /home/analysis/comb_gr_nochr.gtf -tg /home/analysis/comb_gr_nochr_halp_ana.gtf \
-f /home/analysis/liftofftools_features.txt -infer-genes -dir /home/analysis/lotools/ana_halp all \

apptainer exec --no-home --cleanenv --bind /lustre/fs4/zhao_lab/scratch/ulee/:/home /ru-auth/local/home/ulee/scratch_store/docka/liftofftools.sif \
liftofftools -r /home/genomes/dmel-all-chromosome-r6.15.fasta -t /home/genomes/dana-all-chromosome-r1.3.fasta \
-rg /home/analysis/comb_gr_nochr.gtf -tg /home/analysis/comb_gr_nochr_sort_ana.gtf \
-f /home/analysis/liftofftools_features.txt -infer-genes -dir /home/analysis/lotools/ana_hal_sort all \
