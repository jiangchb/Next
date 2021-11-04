SAM=$1
/home/fanyucai/software/samtools/samtools-1.8/bin/samtools view -uS $SAM |/home/fanyucai/software/samtools/samtools-1.8/bin/samtools sort  - -o $SAM.bam && /home/fanyucai/software/samtools/samtools-1.8/bin/samtools index $SAM.bam
