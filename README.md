# umi_tools

######################################
#---------Mutation rates Ppr---------#
######################################



#Raw data processing:
#1. Trim Illumina adapters and UMI adapters from raw reads
#2. Mapping trimmed reads

#(base) linda@bast-work-1:/RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6$ cat uni_tools.sh
for i in A006200178_153621_S1 \
A006200178_153622_S2 \
A006200178_153623_S3 \
A006200178_153624_S4 \
A006200178_153625_S5 \
A006200178_153626_S6
do
umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin=${i}_L002_R2_001.fastq.gz --read2-in=${i}_L002_R1_001.fastq.gz --stdout=${i}_add_barcode_R1.fastq.gz --read2-stdout
umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin=${i}_L002_R2_001.fastq.gz --read2-in=${i}_L002_R3_001.fastq.gz --stdout=${i}_add_barcode_R3.fastq.gz --read2-stdout
done


#umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin=A006200178_153626_S6_L002_R2_001.fastq.gz -I A006200178_153626_S6_L002_R1_001.fastq.gz --read2-in=A006200178_153626_S6_L002_R3_001.fastq.gz --stdout=S6_L002_R1_001.fastq.gz --read2-out=S6_L002_R3_001.fastq.gz
(base) linda@bast-work-1:/RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6$ cat clean.sh
for i in A006200178_153621_S1 \
	A006200178_153622_S2 \
	A006200178_153623_S3 \
	A006200178_153624_S4 \
	A006200178_153625_S5
	#A006200178_153626_S6
do
	#~/Software/QC/TrimGalore-0.6.5/trim_galore -j 30 -q 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --fastqc --paired --output_dir clean_data ${i}_add_barcode_R1.fastq.gz ${i}_add_barcode_R3.fastq.gz 
  #bwa mem meaning: software package for mapping low divergent sequences against a large reverence genome.Fast and accurate. 
	#bwa mem -t 40 /home/shangao/Data/hifiasm_tell-seq/Ppr/Ppr.fa clean_data/${i}_add_barcode_R1_val_1.fq.gz clean_data/${i}_add_barcode_R3_val_2.fq.gz > clean_data/${i}.sam
	#python /home/shangao/script/python/umi_tools_change_reads_title.py -s clean_data/${i}.sam -o clean_data/${i}.replace.sam
	#LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH
	#export LD_LIBRARY_PATH
	#bowtie --threads 40 -v 2 -m 10 -a /home/shangao/Data/hifiasm_tell-seq/Ppr/Ppr -1 clean_data/${i}_add_barcode_R1_val_1.fq.gz -2 clean_data/${i}_add_barcode_R3_val_2.fq.gz  --sam > clean_data/${i}.sam
	#samtools view -bS clean_data/${i}.sam > clean_data/${i}.bam
	#samtools sort clean_data/${i}.bam -o clean_data/${i}.sort.bam
	#samtools index clean_data/${i}.sort.bam
	umi_tools dedup -I clean_data/${i}.sort.bam --output-stats=deduplicated --paired -S clean_data/${i}.sort.de.bam --temp-dir=tmp
done


#download the cleaned files: scp  linda@bast-work-1.zoologie.uni-koeln.de:/RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/clean_data/*.html ./

#Looking at the FastQC report.
