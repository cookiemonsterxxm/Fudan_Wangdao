#BMTagger_Filter
```
#STEP1: ALIGN THE RAW SEQUENCING DATA TO HUMAN GENOME SEQUENCES (BMTAGGER) AND REMOVE THE HOST CONTAMINATIONS

in_dir=/home/xmixu/infant_resistome_2020/2019_Moran_cellhostmicrobes/RAW_DATA
out_dir=/home/xmixu/infant_resistome_2020/2019_Moran_cellhostmicrobes/STEP1_BMTagger
###set the input/output raw sequences directory; subjected to change everytime while you analyze a different study 

module load blast srprism bmtools

Read1=$(ls $in_dir/*fastq | head -n $SLURM_ARRAY_TASK_ID | tail -1)
Read2=$(echo $Read1 | sed 's/_1/_2/')
outname=$(echo $Read1 | cut -f 7 -d '/' | cut -f 1 -d '_')

##the human bitmask file; The Dec. 2013 human reference sequence (GRCh38) was produced by the Genome Reference Consortium: http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/

hm_db=/group/kmkalangrp/databases/genome_references/human_12_16

bmtagger.sh -b $hm_db/allchr_allalt.bitmask  -x $hm_db/allchr_allalt.srprism -q 1 -1 $Read1 -2 $Read2 -o $out_dir/$outname.human.txt

###Remove human DNA 
module load java bbmap

filterbyname.sh in=$Read1 in2=$Read2 out=$out_dir/$outname.R1_nohuman.fastq out2=$out_dir/$outname.R2_nohuman.fastq names=$out_dir/$outname.human.txt include=f

```


#Trimmomatic
```
# USE TRIMMOMATIC TO REMOVE LOW-QUALITY READS
# Take input from BMTagger (removal of the human DNA)
# Trimmomatic parameters followed this paper (Taft et al., 2018, mSphere) and using the paired end mode 
# from link: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

echo "NOW STARTING READ CLEANING WITH TRIMMOMATIC AT: "; date 

#input and output file directory
in_dir=/home/xmixu/infant_resistome_2020/2019_Moran_cellhostmicrobes/STEP1_BMTagger
out_dir=/home/xmixu/infant_resistome_2020/2019_Moran_cellhostmicrobes/STEP2_Trim

#identify the sequences and file names 
file1=$(ls $in_dir/*fastq | head -n $SLURM_ARRAY_TASK_ID | tail -1)
file2=$(echo $file1 | sed 's/R1_/R2_/')
outname=$(echo $file1 | cut -f 7 -d '/' | cut -f 1 -d '.')

module load java trimmomatic

#path to the sequencing adapter sequences 
trim_adapter=/group/kmkalangrp/databases/trimmomatic-0.36_adapters

java -jar /usr/share/java/trimmomatic-0.36.jar PE -phred33 -threads 16 \
-trimlog $out_dir/trimmomatic_log.txt $file1 $file2 \
$out_dir/$outname.R1_paired.fastq $out_dir/$outname.R1_unpaired.fastq.gz $out_dir/$outname.R2_paired.fastq $out_dir/$outname.R2_unpaired.fastq.gz \
ILLUMINACLIP:$trim_adapter/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40

echo "STEP 2 DONE AT: "; date
```

#Fastuniq
```
# STEP 3: Remove duplicated reads
# For FastUniq options: https://wiki.gacrc.uga.edu/wiki/FastUniq

echo "NOW STARTING REMOVING DUPLICATE READS AT: "; date 

#input and output file directory
in_dir=/home/xmixu/infant_resistome_2020/2019_Moran_cellhostmicrobes/STEP2_Trim
out_dir=/home/xmixu/infant_resistome_2020/2019_Moran_cellhostmicrobes/STEP3_FastUniq

#identify the sequences and file names 
file1=$(ls $in_dir/*paired.fastq | head -n $SLURM_ARRAY_TASK_ID | tail -1)
file2=$(echo $file1 | sed 's/R1_/R2_/')
basename=$(echo $file1 | cut -f 7 -d '/' | cut -f 1 -d '.')

touch $out_dir/$basename.input_list.txt
echo $file1 >> $out_dir/$basename.input_list.txt
echo $file2 >> $out_dir/$basename.input_list.txt

#call the program, fastuniq is installed in the conda environment under bio3
module load bio3

fastuniq -i $out_dir/$basename.input_list.txt -t q -o $out_dir/$basename.R1_dedup.fastq -p $out_dir/$basename.R2_dedup.fastq

echo "STEP 3 DONE AT: "; date
```

#FLASH
```
# STEP 5: Merge paired-end reads with FLASH
# The FLASH manual link: http://ccb.jhu.edu/software/FLASH/MANUAL
echo "NOW STARTING PAIRED-END MERGING WITH FLASH AT: "; date

#input and output file directory
in_dir=/home/xmixu/infant_resistome_2020/2019_Moran_cellhostmicrobes/STEP3_FastUniq
out_dir=/home/xmixu/infant_resistome_2020/2019_Moran_cellhostmicrobes/STEP5_FLASH

#identify the sequences and file names 
file1=$(ls $in_dir/*_dedup.fastq | head -n $SLURM_ARRAY_TASK_ID | tail -1)
file2=$(echo $file1 | sed 's/R1_/R2_/')
basename=$(echo $file1 | cut -f 7 -d '/' | cut -f 1 -d '.')

#call the program directly a module 
module load flash

        # -m: minium overlap length 10bp to be similar to pear 
        # -M: max overlap length 
        # -x: mismatch ratio, default is 0.25, which is quite high (e.g: 50bp overlap --> 12.5 mismatch by default)
        
flash $file1 $file2 -m 10 -M 65 -x 0.1 -O -r 100 -f 180 -o $basename -d $out_dir -t 16


echo "STEP 5 DONE AT: "; date


```

#MEGAHIT
```
# STEP 6. Assemble reads into contigs with megaHIT

echo "NOW STARTING ASEEMBLY WITH MEGAHIT AT: "; date

# assemble with paired end reads input (after fastuniq deduplicate)

#input and output file directory
in_dir=/home/xmixu/infant_resistome_2020/2019_Moran_cellhostmicrobes/STEP3_FastUniq
out_dir=/home/xmixu/infant_resistome_2020/2019_Moran_cellhostmicrobes/STEP6_MEGAHIT

#identify the sequences and file names 
file1=$(ls $in_dir/*_dedup.fastq | head -n $SLURM_ARRAY_TASK_ID | tail -1)
file2=$(echo $file1 | sed 's/R1_/R2_/')
basename=$(echo $file1 | cut -f 7 -d '/' | cut -f 1 -d '.')

###call the program, megahit was installed under bio 

module load bio

#-m use up to 95% of the assigned memory; 
#-t run the job using 16 threads 

megahit -1 $file1 -2 $file2 -m 0.95 -t 16 -o $out_dir/$basename.megahit

echo "STEP 6 DONE AT: "; date

```

#ANTIMICROBIAL RESISTANCE ANALYSI WITH MEAGRES
```
# STEP 8: START ANALYZING THE ANTIBIOTIC RESISTOME 
# 
echo "NOW STARTING RESISTOME AT: "; date

#input and output file directory
in_dir=/home/xmixu/infant_resistome_2020/2019_Moran_cellhostmicrobes/STEP3_FastUniq
out_dir=/home/xmixu/infant_resistome_2020/2019_Moran_cellhostmicrobes/STEP8_RESISTOME2
megares_dir=/group/kmkalangrp/databases/antibiotic_references/MEGAResv2.0.0

#identify the sequences and file names 
file1=$(ls $in_dir/*_dedup.fastq | head -n $SLURM_ARRAY_TASK_ID | tail -1)
file2=$(echo $file1 | sed 's/R1_/R2_/')
basename=$(echo $file1 | cut -f 7 -d '/' | cut -f 1 -d '.')

# use the bwa mem method for fastest speed and accuracy
        # BWA aligner manual: http://bio-bwa.sourceforge.net/bwa.shtml
        # -t: thread
module load bwa

# bwa index $megares_dir/megares_modified_database_v2.00_wtsnp.fasta

bwa mem -t 16 $megares_dir/megares_modified_database_v2.00_wtsnp.fasta $file1 $file2 > $out_dir/$basename.align.sam


echo "STEP 8 DONE AT: "; date
```

#METAPHLAN
```
# STEP 10: MetaPhlAn2 TAXONOMIC PROFILING
# The MetaPhlAn2 help page: https://bitbucket.org/biobakery/metaphlan2/src/default/

echo "NOW STARTING METAGENOMIC TAXONOMY PROFILING AT: "; date

#input and output file directory
in_dir=/home/xmixu/infant_resistome_2020/2019_Moran_cellhostmicrobes/STEP5_FLASH
out_dir=/home/xmixu/infant_resistome_2020/2019_Moran_cellhostmicrobes/STEP10_TAXONOMY
# metaphlan2_DB=/home/jxnliu/metaphlan2020/metaphlan2/db_v20

#identify the sequences and file names 
file=$(ls $in_dir/*.extendedFrags.fastq | head -n $SLURM_ARRAY_TASK_ID | tail -1)

basename=$(echo $file | cut -f 7 -d '/' | cut -f 1 -d '.')

module load bio
#"module load bio" includes metaphlan v2.7.7, bowtie2 v2.3.4.3, and numpy v1.15.4 all built in.
# --mpa_pkl MPA_PKL     The metadata pickled MetaPhlAn file [deprecated]
# --bowtie2db METAPHLAN_BOWTIE2_DB; The BowTie2 database file of the MetaPhlAn database. 
# Used if --input_type is fastq, fasta, multifasta, or multifastq [default /share/apps/bio/bio/bin/metaphlan_databases]

metaphlan2.py $file \
--input_type fastq --nproc 20 > $out_dir/$basename.metaphlan.txt

echo "STEP 10 DONE AT: "; date
```

