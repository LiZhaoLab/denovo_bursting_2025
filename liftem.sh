#! /bin/bash

# srun -p zhao halLiftover dro1.hal mel ../comb_gr_nochr.bed sim ../comb_gr_nochr_sim.bed
# srun -p zhao halLiftover dro1.hal mel ../comb_gr_nochr.bed sec ../comb_gr_nochr_sec.bed
# srun -p zhao halLiftover dro1.hal mel ../comb_gr_nochr.bed yak ../comb_gr_nochr_yak.bed
# srun -p zhao halLiftover dro1.hal mel ../comb_gr_nochr.bed ere ../comb_gr_nochr_ere.bed
# srun -p zhao halLiftover dro1.hal mel ../comb_gr_nochr.bed ana ../comb_gr_nochr_ana.bed

declare -a arr=("sim"
"sec"
"yak"
"ere"
"ana")

# for i in "${arr[@]}"
# do
	# bedtools sort -i ../comb_gr_nochr_$i.bed|bedtools groupby -g 4,1,5,6 -c 2,3 -o min,max|awk -v OFS="\t" '{print $2,$5,$6,$1,$3,$4}' > ../comb_gr_nochr_sort_$i.bed
# done


#halper
#peaksummits

for i in "${arr[@]}"
do
	# halAlignmentDepth --outWiggle ../mel_"$i"_aln_dep.wig dro1.hal $i
	python3 -m getMaxScorePositionFromWig --bedFileName ../comb_gr_nochr.bed --wigFileName ../mel_"$i"_aln_dep.wig --chromSizesFileName /ru-auth/local/home/ulee/scratch_store/genomes/d"$i"/d"$i".chrom.sizes --highestScoreLocationFileName ../comb_gr_nochr_"$i"_sfile.bed --tempDir ~/tmp/
	halLiftover dro1.hal mel ../comb_gr_nochr_"$i"_sfile.bed $i ../comb_gr_nochr_"$i"_sfile_lo.bed
	python3 -m orthologFind -max_len 100000 -min_len 50 -protect_dist 5 -qFile ../comb_gr_nochr.bed -tFile ../comb_gr_nochr_"$i".bed -sFile ../comb_gr_nochr_"$i"_sfile_lo.bed -oFile ../comb_gr_nochr_halp_"$i".bed -narrowPeak -mult_keepone
done


