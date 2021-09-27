# umi_tools

######################################
#---------Mutation rates Ppr---------#
######################################


#Raw data processing:
#1. Trim Illumina adapters and UMI adapters from raw reads

#pwd
/home/linda/Scratch

#new file in Scratch: mkdir trimmap
#pwd
/home/linda/Scratch/trimmap

#cp reads in my directory for S1-S6
cp /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/A006200178_153621_S1_L* ./


#to look into script:
cat /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/uni_tools.sh
#for i in A006200178_153621_S1 \
#A006200178_153622_S2 \
#A006200178_153623_S3 \
#A006200178_153624_S4 \
#A006200178_153625_S5 \
#A006200178_153626_S6
#do
#umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin=${i}_L002_R2_001.fastq.gz --read2-#in=${i}_L002_R1_001.fastq.gz --stdout=${i}_add_barcode_R1.fastq.gz --read2-stdout
#umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin=${i}_L002_R2_001.fastq.gz --read2-#in=${i}_L002_R3_001.fastq.gz --stdout=${i}_add_barcode_R3.fastq.gz --read2-stdout
#done

#command
/RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/uni_tools.sh
 
 #Result: extracted barcodes, reads are trimmed

-----
#2. Mapping trimmed reads:

#to look into script: 
cat /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6
/clean_linda.sh
#for i in A006200178_153621_S1 \
#	A006200178_153622_S2 \
#	A006200178_153623_S3 \
#	A006200178_153624_S4 \
#	A006200178_153625_S5
#	#A006200178_153626_S6
#do
#/NVME/Software/QC/TrimGalore-0.6.5/trim_galore -j 30 -q 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --fastqc --paired --output_dir clean_data ${i}_add_barcode_R1.fastq.gz ${i}_add_barcode_R3.fastq.gz
        #cd clean_data
        #bwa mem -t 40 /home/shangao/Data/hifiasm_tell-seq/Ppr/Ppr.fa ./${i}_add_barcode_R1_val_1.fq.gz ./${i}_add_barcode_R3_val_2.fq.gz > ./${i}.sam
        #python /home/shangao/script/python/umi_tools_change_reads_title.py -s ./${i}.sam -o ./${i}.replace.sam
        #LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH
        #export LD_LIBRARY_PATH
        #bowtie --threads 40 -v 2 -m 10 -a /home/shangao/Data/hifiasm_tell-seq/Ppr/Ppr -1 ./${i}_add_barcode_R1_val_1.fq.gz -2 ./${i}_add_barcode_R3_val_2.fq.gz  --sam > ./${i}.sam
        #samtools view -bS ./${i}.sam > ./${i}.bam
        #samtools sort ./${i}.bam -o ./${i}.sort.bam
        #samtools index ./${i}.sort.bam
        #umi_tools dedup -I ./${i}.sort.bam --output-stats=deduplicated --paired -S ./${i}.sort.de.bam --temp-dir=tmp
#done


#command: /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6
/clean_linda.sh


#3. Download the cleaned files: scp  linda@bast-work-1.zoologie.uni-koeln.de:/RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/clean_data/*.html ./

#4. Looking at the FastQC report.
