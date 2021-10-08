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
#2. Mapping trimmed reads to the reference genome:
#- the mapped reads were transformed into bam format and sorted with samtools

#to look into script: 
cat /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/clean_linda.sh
#for i in A006200178_153621_S1 \
#A006200178_153622_S2 \
#A006200178_153623_S3 \
#A006200178_153624_S4 \
#A006200178_153625_S5 \
#A006200178_153626_S6
#do
        #/NVME/Software/QC/TrimGalore-0.6.5/trim_galore -j 30 -q 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a2 AGATCGGAA$
        bwa mem -t 40 /home/shangao/Data/juicer/Ppr/Ppr_withoutchange/review/Ppr.FINAL.sort.fasta /home/linda/Data/map$
        #python /home/shangao/script/python/umi_tools_change_reads_title.py -s ./${i}.sam -o ./${i}.replace.sam
        #LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH
        #export LD_LIBRARY_PATH
        #bowtie --threads 40 -v 2 -m 10 -a /home/shangao/Data/juicer/Ppr/Ppr_withoutchange/review/Ppr.FINAL.sort -1 ./$
        #samtools view -bS ./${i}.sam > ./${i}.bam
        #samtools sort ./${i}.bam -o ./${i}.sort.bam
        #samtools index ./${i}.sort.bam
        #umi_tools dedup -I ./${i}.sort.bam --output-stats=deduplicated --paired -S ./${i}.sort.de.bam --temp-dir=tmp
#done

#command: sh clean_linda.sh


#3. Download the cleaned files: scp  linda@bast-work-1.zoologie.uni-koeln.de:/home/linda/Scratch/trimmap/clean_data/*.html ./
#and looking in the FastQC report.


#4. SNP calling with GATK4:
#genotype each line seperately with HaplotypeCaller

#to look into script:
cat/RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/clean_data/snp/gatk.sh 
###software
#gatk=/NVME/Software/popgen/gatk-4.1.9.0/gatk
###data
#ref=/home/shangao/Data/juicer/Ppr/test/review/blobtools/genome/Ppr.FINAL.fa
#bam=$1
#out1=$2
#out2=$3
###build dict for genome
#/NVME/Software/popgen/gatk-4.1.9.0/gatk CreateSequenceDictionary -R Ppr.FINAL.fa -O Ppr.FINAL.dict
###call gvcf
#$gatk HaplotypeCaller \
       # -R $ref \
       # --emit-ref-confidence GVCF \
       #-I $bam \
       # -O $out1

###detect SNPs
#$gatk GenotypeGVCFs \
       # -R $ref \
       # -V $out1 \
       # -O $out2

###compress
#bgzip -f $out2
#tabix -p vcf $out2.gz

#command: sh /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/clean_data/snp/gatk.sh


#5. Setting hard filters with parameters as follows:
#QD (variant confidence standardized by depth) <2
#MQ (mapping quality of a SNP) <30
#FS (strand bias in support for REF versus ALT allele calls) >60
#SOR (sequencing bias in which one DNA strand is favored over the other) >5
#ReadPosRankSumTest < -8

#also Indel:
#QD <2
#FS >100
#SOR >5
#ReadPosRankSumTest < -8

#to look into script:
cat /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/clean_data/snp/filter.sh
#vcf=$1
#snpvcf=$2
#indelvcf=$3
#filterSNP=$4
#filterINDEL=$5
#finalvcf=$6
###SelectVariants SNP
#$gatk SelectVariants \
       # -select-type SNP \
       # -V $vcf \
       # -O $snpvcf
###SelectVariants INDEL
#$gatk SelectVariants \
       # -select-type INDEL \
       # -V $vcf \
       # -O $indelvcf
###filter SNP
#$gatk VariantFiltration \
       # -V $snpvcf \
       # --filter-expression "QD <2.0 || MQ <30.0 || FS >60.0 || SOR >5.0 || ReadPosRankSum < -8.0" \
       # --filter-name "PASS" \
       # -O $filterSNP
###filter INDEL
#$gatk VariantFiltration \
       # -V $indelvcf \
       # --filter-expression "QD <2.0 || FS >100.0 || SOR >5.0 || ReadPosRankSum < -8.0" \
       # --filter-name "PASS" \
       # -O $filterINDEL
###merge SNP INDEL
#$gatk MergeVcfs \
       # -I $filterSNP \
       # -I $filterINDEL \
       # -O $finalvcf
###delete temp
#rm -f $snpvcf $indelvcf $filterSNP $filterINDEL

#command: sh filter.sh
#output: VCF file

#6. VCFtools for mutation rates: comparing mothers and daughters to get the mutation rate
#to look into script:
cat /RAID/Data/mites/reads/linda_umi/bastet.ccg.uni-koeln.de/downloads/jbast_JB03_September6/clean_data/snp/merged_gvcf/filter.BIALLELIC.sh 
#/NVME/Software/popgen/gatk-4.1.9.0/gatk SelectVariants  -V merge.snp.f.vcf.gz --restrict-alleles-to BIALLELIC -O merge.snp.f.bi.vcf.gz

#command: sh filter.BIALLELIC.sh

