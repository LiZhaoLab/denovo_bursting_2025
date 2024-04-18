#!/usr/bin/bash


declare -a arr_ref=("dana-all-chromosome-r1.3.fasta"
"dere-all-chromosome-r1.3.fasta"
"dyak-all-chromosome-r1.3.fasta"
"dsec-all-chromosome-r1.3.fasta"
"dsim-all-chromosome-r2.02.fasta"
"dmel-all-chromosome-r6.15.fasta")

declare -a arr=("SRR1024003_trimmed.fastq" 
"SRR1024002_trimmed.fastq"
"SRR1617567_trimmed.fastq"
"SRR486109_trimmed.fastq" 
"SRR5839446_1_trimmed.fastq\ SRR5839446_2_trimmed.fastq"
"SRR16213299_1_trimmed.fastq\ SRR16213299_2_trimmed.fastq" 
"SRR16213300_1_trimmed.fastq\ SRR16213300_2_trimmed.fastq"
"SRR6667440_1_trimmed.fastq\ SRR6667440_2_trimmed.fastq" 
"SRR6667439_1_trimmed.fastq\ SRR6667439_2_trimmed.fastq"
"SRR931549_1_trimmed.fastq\ SRR931549_2_trimmed.fastq")

declare -a arr_clean=("SRR1024003" 
"SRR1024002"
"SRR1617567"
"SRR486109" 
"SRR5839446"
"SRR16213299" 
"SRR16213300"
"SRR6667440" 
"SRR6667439"
"SRR931549")

declare -a arr_spec=("/ru-auth/local/home/ulee/scratch_store/genomes/dana-all-chromosome-r1.3"
"/ru-auth/local/home/ulee/scratch_store/genomes/dana-all-chromosome-r1.3"
"/ru-auth/local/home/ulee/scratch_store/genomes/dere-all-chromosome-r1.3"
"/ru-auth/local/home/ulee/scratch_store/genomes/dyak-all-chromosome-r1.3"
"/ru-auth/local/home/ulee/scratch_store/genomes/dyak-all-chromosome-r1.3"
"/ru-auth/local/home/ulee/scratch_store/genomes/dsec-all-chromosome-r1.3"
"/ru-auth/local/home/ulee/scratch_store/genomes/dsec-all-chromosome-r1.3"
"/ru-auth/local/home/ulee/scratch_store/genomes/dsim-all-chromosome-r2.02"
"/ru-auth/local/home/ulee/scratch_store/genomes/dsim-all-chromosome-r2.02"
"/ru-auth/local/home/ulee/scratch_store/genomes/dmel-all-chromosome-r6.15")

declare -a arr_gtf=("../comb_gr_nochr_halp_ana.bed"
"../comb_gr_nochr_halp_ana.bed"
"../comb_gr_nochr_halp_ere.bed"
"../comb_gr_nochr_halp_yak.bed"
"../comb_gr_nochr_halp_yak.bed"
"../comb_gr_nochr_halp_sec.bed"
"../comb_gr_nochr_halp_sec.bed"
"../comb_gr_nochr_halp_sim.bed"
"../comb_gr_nochr_halp_sim.bed"
"../comb_gr_nochr.bed")

declare -a arr_fb=("../annots/dana-all-r1.3-fb.gff"
"../annots/dana-all-r1.3-fb.gff"
"../annots/dere-all-r1.3-fb.gff"
"../annots/dyak-all-r1.3-fb.gff"
"../annots/dyak-all-r1.3-fb.gff"
"../annots/dsec-all-r1.3-fb.gff"
"../annots/dsec-all-r1.3-fb.gff"
"../annots/dsim-all-r2.02-fb.gff"
"../annots/dsim-all-r2.02-fb.gff"
"../comb_gr_nochr.gtf")

declare -a arr_liftoff=("../comb_gr_nochr_liftoff_ana.gff"
"../comb_gr_nochr_liftoff_ana.gff"
"../comb_gr_nochr_liftoff_ere.gff"
"../comb_gr_nochr_liftoff_yak.gff"
"../comb_gr_nochr_liftoff_yak.gff"
"../comb_gr_nochr_liftoff_sec.gff"
"../comb_gr_nochr_liftoff_sec.gff"
"../comb_gr_nochr_liftoff_sim.gff"
"../comb_gr_nochr_liftoff_sim.gff"
"../comb_gr_nochr.gtf")

declare -a arr_hal=("../comb_gr_nochr_ana_rmerge.gtf"
"../comb_gr_nochr_ana_rmerge.gtf"
"../comb_gr_nochr_ere_rmerge.gtf"
"../comb_gr_nochr_yak_rmerge.gtf"
"../comb_gr_nochr_yak_rmerge.gtf"
"../comb_gr_nochr_sec_rmerge.gtf"
"../comb_gr_nochr_sec_rmerge.gtf"
"../comb_gr_nochr_sim_rmerge.gtf"
"../comb_gr_nochr_sim_rmerge.gtf"
"../comb_gr_nochr.gtf")

cd /ru-auth/local/home/ulee/scratch_store/genomes

# for i in "${arr_ref[@]}"
# do
	# cp $i /ru-auth/local/home/ulee/bowtie2-2.5.1-linux-x86_64
	# bowtie2-build $i "$(basename $i .fasta)" &
# done

# wait;

cd /ru-auth/local/home/ulee/scratch_store/analysis/txn_testis_outgroup/trimmed/

# mkdir fastqc
# mkdir trimmed

# for i in *.fastq
# do
	# fastqc -o fastqc $i
	# outfile="$(basename $i .fastq)"_trimmed.fastq
	# java -jar /ru-auth/local/home/ulee/Trimmomatic-main/dist/jar/trimmomatic-0.40-rc1.jar SE -phred33 -threads 2 $i trimmed/$outfile ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
# done 

# mkdir map_tophat

# for i in ${!arr[@]}
# do
	# outfile="${arr_clean[$i]}"_map
	# tophat -o ../$outfile "${arr_spec[$i]}" "${arr[$i]}"  
# done

# tophat -o ../SRR5839446_map /ru-auth/local/home/ulee/scratch_store/genomes/dyak-all-chromosome-r1.3 SRR5839446_1_trimmed.fastq SRR5839446_2_trimmed.fastq
# tophat -o ../SRR16213299_map /ru-auth/local/home/ulee/scratch_store/genomes/dsec-all-chromosome-r1.3 SRR16213299_1_trimmed.fastq SRR16213299_2_trimmed.fastq
# tophat -o ../SRR16213300_map /ru-auth/local/home/ulee/scratch_store/genomes/dsec-all-chromosome-r1.3 SRR16213300_1_trimmed.fastq SRR16213300_2_trimmed.fastq
# tophat -o ../SRR6667440_map /ru-auth/local/home/ulee/scratch_store/genomes/dsim-all-chromosome-r2.02 SRR6667440_1_trimmed.fastq SRR6667440_2_trimmed.fastq
# tophat -o ../SRR6667439_map /ru-auth/local/home/ulee/scratch_store/genomes/dsim-all-chromosome-r2.02 SRR6667439_1_trimmed.fastq SRR6667439_2_trimmed.fastq
# tophat -o ../SRR931549_map /ru-auth/local/home/ulee/scratch_store/genomes/dmel-all-chromosome-r6.15 SRR931549_1_trimmed.fastq SRR931549_2_trimmed.fastq

cd /ru-auth/local/home/ulee/scratch_store/analysis/txn_testis_outgroup/

# for i in ${!arr_clean[@]}
# do
	# arr_dir="${arr_clean[$i]}_map"
	# arr_file="$arr_dir"/accepted_hits.bam
	
	# stringtie -o $arr_dir/output.gtf $arr_file 
# done

# bedToGenePred ../comb_gr_nochr.bed /dev/stdout | genePredToGtf file /dev/stdin ../comb_gr_nochr.gtf

# stringtie -o SRR931549_map/output_withref.gtf -G ../comb_gr_nochr.gtf SRR931549_map/accepted_hits.bam

for i in ${!arr_clean[@]}
do
	arr_dir="${arr_clean[$i]}_map"
	arr_file="$arr_dir"/accepted_hits.bam
	gtf_clean="$(basename ${arr_gtf[$i]} .bed)"
	
	# Rscript ../narrowPeak_for_ucsc.R ${arr_gtf[$i]} "${gtf_clean}_fix.bed"
	
	# bedToGenePred ${gtf_clean}_fix.bed /dev/stdout | genePredToGtf file /dev/stdin ../"${gtf_clean}.gtf"

	# stringtie -o $arr_dir/output_withref_halp.gtf -G ../"${gtf_clean}.gtf" -e -A $arr_dir/abund_withref_halp.csv $arr_file 
	# stringtie -o $arr_dir/output_withref_halp_deno.gtf -G ../"${gtf_clean}.gtf" -A $arr_dir/abund_withref_halp_denovo.csv $arr_file 
	
	stringtie -o $arr_dir/output_fb_halp.gtf -G "${arr_fb[$i]}" -e -A $arr_dir/abund_fb_halp.csv $arr_file
	stringtie -o $arr_dir/output_fb_halp_deno.gtf -G "${arr_fb[$i]}" -A $arr_dir/abund_fb_halp_denovo.csv $arr_file 

	# stringtie -o $arr_dir/output_liftoff.gtf -G "${arr_liftoff[$i]}" -e -A $arr_dir/abund_liftoff.csv $arr_file
	# stringtie -o $arr_dir/output_liftoff_deno.gtf -G "${arr_liftoff[$i]}" -A $arr_dir/abund_liftoff_denovo.csv $arr_file 
	
	# stringtie -o $arr_dir/output_hal.gtf -G "${arr_hal[$i]}" -e -A $arr_dir/abund_hal.csv $arr_file
	# stringtie -o $arr_dir/output_hal_deno.gtf -G "${arr_hal[$i]}" -A $arr_dir/abund_hal_denovo.csv $arr_file 
done

