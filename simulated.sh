dir=
pbsim3 --strategy wgs --method qshmm --qshmm $dir/QSHMM-RSII.model --depth 20 --genome ./ragtag.sca.fasta --pass-num 10

#n=number of chromosome
samtools view -@32 -b sd.sam > sd.bam

ccs sd.bam sd.ccs.bam

pbmerge -o sd.ccs.bam sd*.ccs.bam
pbindex sd.ccs.bam

samtools fastq sd.ccs.bam > output.fastq
