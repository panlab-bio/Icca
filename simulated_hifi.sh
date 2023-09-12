###
Generating simulated hifi data using pbsim3
###

#1 Run PBSIM with the sample data
   #To run multi-pass sequencing:
 pbsim --strategy wgs
      --method qshmm
      --qshmm data/QSHMM-RSII.model
      --depth 20
      --genome sample/sample.fasta
      --pass-num 10
#PBSIM does not directly simulate HiFi reads, but only simulate generation of CLR by multi-pass sequencing;
#the output of PBSIM simulation (SAM format data) is converted to BAM format data, and input to ccs software(https://github.com/PacificBiosciences/ccs),which generates HiFi reads.

#2 Convert the compared sam file into a bam file
samtools view [options] <in.bam>|<in.sam>
#eg: samtools view -@32 -b sd_0001.sam > sd_0001.bam

#3 ccs combines multiple subreads of the same SMRTbell molecule using a statistical model to produce one highly accurate consensus sequence, also called a HiFi read, along with base quality values.
ccs movie.subreads.bam movie.ccs.bam
# You can parallelize by chunking. which can be installed with conda install pbbam
pbindex movie.subreads.bam
# An example workflow, all ccs invocations can run simultaneously:
ccs movie.subreads.bam sd.ccs.1.bam --chunk 1/10 -j <THREADS>
ccs movie.subreads.bam sd.ccs.2.bam --chunk 2/10 -j <THREADS>
...
ccs movie.subreads.bam sd.ccs.10.bam --chunk 10/10 -j <THREADS>

#4 Merge chunks with pbmerge and index with pbindex
pbmerge -o sd.ccs.bam sd*.ccs.bam
pbindex sd.ccs.bam

#5 bam2fq
samtools fastq sd.ccs.bam > output.fastq
