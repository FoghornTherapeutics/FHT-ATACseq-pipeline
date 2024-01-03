#!/bin/bash

# Run inside the directory with the name of the experiment_id containing one folder: raw_data


set -e  # exit script when hitting an error (non-zero exit code)
#set -x # debug mode when needed


if [[ -e alignment ]]
then
	rm -r alignment
fi
mkdir alignment

cd alignment


##### Copy code from ATACseqCode #####
ln -s ~/code/bioinformatics-NGS-and-Analysis/NGS/ATAC_seq/alignment src_code
cp src_code/atac_seq_snakemake_conf.json .


##### Copy fastqs from raw_data #####
mkdir fastqs code  #star_output fastqc code log_files
cd fastqs

ln -s ../../raw_data/*.fq.gz .
ln -s ../../raw_data/*.fastq.gz .

#ln -s /mnt/wkdir/raw_data/raw_data/*/*.fq.gz *.fastq.gz .

#for a in *.fq.gz; do b=$(echo $a | awk -F _ '{print $1""$2"_"$NF}'); echo $a $b; done #print new names next to old names
#for a in *.fq.gz; do b=$(echo $a | awk -F _ '{print $1"" $2"_"$NF}'); echo $a $b; mv $a $b; done # rename files


##### Create sample_group_IDs #####
#echo check for pairs
#ls -1 | sed 's/_[1,2].fq.gz//' | grep -v check_for_paired_fastq | sort | uniq -c &> check_for_paired_fastq
#cat check_for_paired_fastq 
#echo check if any of above are not 2:
#grep -v " 2 " check_for_paired_fastq &> num_check_for_paired_fastq || true  # added || true at the end because grep has exit code of 1 when nothing found which is the normal expected behavior for this search
#cat num_check_for_paired_fastq
#wc -l num_check_for_paired_fastq
#mv check_for_paired_fastq num_check_for_paired_fastq ../code

#ls -1 | sed 's/_[1,2].fq.gz//' | sort -u > ../sample_IDs.grp

#cd ..
