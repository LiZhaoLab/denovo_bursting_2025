#! /bin/bash

declare -a arr=("SRR1024003"
"SRR1024002"
"SRR1617567"
"SRR486109"
"SRR5839446_1"
"SRR5839446_2"
"SRR16213299_1"
"SRR16213299_2"
"SRR16213300_1"
"SRR16213300_2"
"SRR6667440_1"
"SRR6667440_2"
"SRR6667439_1"
"SRR6667439_2"
"SRR931549_1"
"SRR931549_2")

declare -a arr_spec=("/ru-auth/local/home/ulee/scratch_store/genomes/dana-all-chromosome-r1.3"
"/ru-auth/local/home/ulee/scratch_store/genomes/dana-all-chromosome-r1.3"
"/ru-auth/local/home/ulee/scratch_store/genomes/dere-all-chromosome-r1.3"
"/ru-auth/local/home/ulee/scratch_store/genomes/dyak-all-chromosome-r1.3"
"/ru-auth/local/home/ulee/scratch_store/genomes/dyak-all-chromosome-r1.3"
"/ru-auth/local/home/ulee/scratch_store/genomes/dyak-all-chromosome-r1.3"
"/ru-auth/local/home/ulee/scratch_store/genomes/dsec-all-chromosome-r1.3"
"/ru-auth/local/home/ulee/scratch_store/genomes/dsec-all-chromosome-r1.3"
"/ru-auth/local/home/ulee/scratch_store/genomes/dsec-all-chromosome-r1.3"
"/ru-auth/local/home/ulee/scratch_store/genomes/dsec-all-chromosome-r1.3"
"/ru-auth/local/home/ulee/scratch_store/genomes/dsim-all-chromosome-r2.02"
"/ru-auth/local/home/ulee/scratch_store/genomes/dsim-all-chromosome-r2.02"
"/ru-auth/local/home/ulee/scratch_store/genomes/dsim-all-chromosome-r2.02"
"/ru-auth/local/home/ulee/scratch_store/genomes/dsim-all-chromosome-r2.02"
"/ru-auth/local/home/ulee/scratch_store/genomes/dmel-all-chromosome-r6.15"
"/ru-auth/local/home/ulee/scratch_store/genomes/dmel-all-chromosome-r6.15"
)

cd /ru-auth/local/home/ulee/scratch_store/genomes/

for i in *.fasta
do
	outfile="$(basename $i .fasta)"
	STAR --runMode genomeGenerate --genomeDir $outfile --genomeFastaFiles $i --runThreadN 2 --genomeSAindexNbases 12
done

cd /ru-auth/local/home/ulee/scratch_store/analysis/txn_testis_outgroup
mkdir fastqc
mkdir trimmed

for i in *.fastq
do
	fastqc -o fastqc $i
	outfile="$(basename $i .fastq)"_trimmed.fastq
	java -jar /ru-auth/local/home/ulee/Trimmomatic-main/dist/jar/trimmomatic-0.40-rc1.jar SE -phred33 -threads 2 $i trimmed/$outfile ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
done 

mkdir mapped

for i in ${!arr[@]}
do
	infile=trimmed/"${arr[$i]}_trimmed.fastq"
	outfile="${arr[$i]}"_map_
	STAR --genomeDir "${arr_spec[$i]}" --runThreadN 2 --readFilesIn $infile --outFileNamePrefix mapped/$outfile --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --outSAMattributes All
done

cd /ru-auth/local/home/ulee/scratch_store/analysis/txn_testis_outgroup/mapped/

for i in *.bam
do
	outfile="$(basename $i .bam)".bed
	bedtools bamtobed -i $i > $outfile
done