#!/bin/bash
###
### Sam2Bam
### Author:Congjia.chen
### Usage: Sam2Bam
###   sh sam2bam.sh <SAM FILE> 
### Options:
###   <SAM FILE>   SAM file
###   -h        Show this message.

help() {
	awk -F'### ' '/^###/ { print $2 }' "$0"
}

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
	help
	exit 1
fi
SAM=$1
/home/fanyucai/software/samtools/samtools-1.8/bin/samtools view -uS $SAM |/home/fanyucai/software/samtools/samtools-1.8/bin/samtools sort  - -o $SAM.bam && /home/fanyucai/software/samtools/samtools-1.8/bin/samtools index $SAM.bam
