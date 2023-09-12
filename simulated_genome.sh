'''
Based on the simulation haploid genome, 
the script SimSID.py is used to artificially construct a diploid genome by introducing different numbers of SNPs, inser-tions, and deletions using various parameters
'''
#hap1
SimSID.py --random_length -v -r ./name.sca.fasta -o name_hap1
#hap2
#SimSID.py [-h] [-s SNP] [-i INSERTION] [--insert_length INSERT_LENGTH] [-d DELETION] [--delete_length DELETE_LENGTH] [--random_length] [-v] -r REF -o OUT
SimSID.py -s 0.02 -i 0.02 --insert_length 12 -d 0.02 --delete_length 12 --random_length -v -r ./name.sca.fasta -o name_hap2
