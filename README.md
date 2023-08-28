# Icca
We propose a novel method that is able to identify and fill the gaps in the chromosome-level assembly by recalling the sequences in the lost small contigs.

The algorithm comprises four phases.To organize and prevent data disorder, we're setting up four folders to execute scripts for their respective stages.
In the first phase, the completeness of the initial assembly results is assessed using BUSCO, and a list of single-copy genes missing in the chromosome-level assembly but appearing in contigs is identified. 
The second phase locates these missing single-copy genes in contigs and identified the entire missing regions. 
In the third phase, the HiFi reads related to the missing genomic regions are recalled, and these reads are then input into an assembly tool to regenerate new sequences. 
In the fourth phase, the newly obtained sequences are aligned one by one with the chromosome-level assembly, gradually in-creasing the level of completeness. 
