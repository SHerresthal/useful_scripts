#!/bin/bash

# This script goes to a directory, looks for .gz files in all subfolders and merges all files that share the first identifier. 
# For example, 3624_S6_L001_R1_001.fastq.gz and 3624_S6_L001_R1_002.fastq.gz is merged to 3624_S6_L001_R1.fastq. 

# go to the directory that contains subfolders for each run
cd /dir/fastq
# make a new directory for the merged runs
mkdir ./merged_runs

for filename in $( find -name *.gz -type f | sort)
do
  echo "filename = $filename"
  filebase=$(basename $filename)
  IFS='_' read -r -a parts <<< "$filebase"
  echo "parts[0] = ${parts[0]}"
  gunzip -c $filename >> "./merged_runs/${parts[0]}_${parts[1]}_${parts[2]}_${parts[3]}.fastq"
done


exit 0
