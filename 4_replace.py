#!/usr/bin/python
import os
import glob
import pandas as pd
import numpy as np
from collections import defaultdict
import sys, getopt
import os.path
import subprocess
from Bio import SeqIO
import multiprocessing


os.system("ln -s ../01_pre/*_sca.fa ./sca.fa")
os.system("ln -s ../01_pre/*_odb10 ./")
os.system("ln -s ../03_recall_final/*_test/*_single_sequence.fasta ./")
os.system("ln -s ../03_recall_final/all_success.bed ./all_locate_gene.info")

bed_file = "all_locate_gene.info"
bed_dict = {}
with open("all_locate_gene.info", "r") as f:
    for line in f:
        line = line.strip().split("\t")
        print(line,type(line))
        key=line[0]
        #value = [line[1], line[2], line[3]]
        value = [line[2], line[3], line[4]]
                #seq_name  start     end
        bed_dict[key] = value

def read_genome(genome_file):
    genome = {}
    with open(genome_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                seq_name = line.strip().lstrip(">")
                genome[seq_name] = ""
            else:
                genome[seq_name] += line.strip()
    return genome

def replace_contig(Genome_file, contig_file, sort_filter_map_paf,step,key):
    
    file1 = Genome_file
    file_name1 = os.path.basename(file1)
    
    file2 = contig_file
    file_name2 = os.path.basename(file2)
    
    file3 = sort_filter_map_paf
    file_name3 = os.path.basename(file3)
    
    print("file_name1:",file_name1,"\n","file_name2:",file_name2,"\n","file_name3:",file_name3) 
    
    genome = read_genome(Genome_file)
    
    with open(contig_file, "r") as f:
        contig_seq = f.read().strip().split("\n")[1]

    with open(sort_filter_map_paf, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            key = fields[8]

            chr_name = fields[0]
            seq_name = fields[4]

            start1 = int(fields[2])
            end1 = int(fields[3])
            
            seq_len = int(fields[5])
            
            start2 = int(fields[6])
            end2 = int(fields[7])
            
            print("chr_name:",chr_name,"start1:",start1,"end1:",end1,"end1-start1",end1-start1)
            print("seq_name:",seq_name,"start2:",start2,"end2:",end2,"end2-start2",end2-start2)
            
            seq1 = genome[chr_name][:start1-1]
            print("seq1:",len(seq1))
            
            seq2 = contig_seq[start2:end2]
            print("seq2:",len(seq2))
             
            seq3 = genome[chr_name][end1:]
            print("seq3:",len(seq3))
            
            v_con, g_s, g_e = bed_dict[key]
            if(int(g_s)>=start2 and int(g_e)<=end2):
                genome[chr_name] = genome[chr_name][:start1] + contig_seq[start2:end2] + genome[chr_name][end1:]
                print("genome:",len(genome[chr_name]))
    
    new_genome_file = "replaced_"+str(step)+"_sca.fa"  # raw:sca
    
    with open(new_genome_file, "w") as f:
        for seq_name, seq in genome.items():
            print("seq_name:",seq_name)
            print("len:",len(seq))
            f.write(f">{seq_name}\n{seq}\n")
            
    print("Inside the function new_genome_file:",os.path.basename(new_genome_file))        
    return new_genome_file

step=0
with open(bed_file, "r") as f:
    for line in f:
        print("line:",line)
        fields = line.strip().split("\t")
        
        key = fields[0]
        seq_name = fields[2]
        
        start = int(fields[3])
        end = int(fields[4])
        
        #minimap
        contig_file = f"{key}_{seq_name}_single_sequence.fasta"
        map_paf = f"{key}_{seq_name}_sca_tqcon_aln.paf"
        
        if(step==0):
          genome_file="./sca.fa"
        else:
          genome_file=new_genome_file
            
        command = f"minimap2 -x asm5 -c {genome_file} {contig_file} > {map_paf}"
        subprocess.run(command, shell=True)
        print("minmap2 finish")
        
        sort_filter_map_paf = f"{key}_{seq_name}_sort_sca_tqcon_aln.paf"
        
        find_flag=0
        f_n = map_paf 
        with open(f_n, 'r') as f_in, open(sort_filter_map_paf, 'w') as f_out: 
            max_value = None
            max_row = None
            
            for line in f_in:
                line_list = line.strip().split('\t')
            
                map_con_s = int(line_list[2])
                map_con_e = int(line_list[3])
                    
                if (start >= map_con_s and end <= map_con_e):
                    max_row = line
                    find_flag=1
                    break
            if(find_flag):     
                max_row = max_row.strip().split("\t")
                
                contig_name = max_row[0]
                contig_length = max_row[1]
                contig_start = max_row[2]
                contig_end = max_row[3]
                
                chr_name   = max_row[5]
                chr_length = max_row[6]
                chr_start  = max_row[7]
                chr_end    = max_row[8]
                
                f_out.write(f'{chr_name}\t{chr_length}\t{chr_start}\t{chr_end}\t{contig_name}\t{contig_length}\t{contig_start}\t{contig_end}\t{key}\n')       
                
        if(find_flag):     
            #replace
            new_genome_file=replace_contig(genome_file, contig_file, sort_filter_map_paf,step,key)
            print("Outside the function new_genome_file:",os.path.basename(new_genome_file))
            step+=1

#command = f"busco -i ./sca.fa -c 60 -o sca_busco -m geno -l ./*_odb10 --offline"
#os.system(command)

def busco_check(filename):
    command = f"busco -i {filename} -c 60 -o {filename}_busco -m geno -l ./*_odb10 --offline"
    os.system(command)

filenames = ['./replaced_{}_sca.fa'.format(i) for i in range(step)]
with multiprocessing.Pool(processes=4) as pool:
    pool.map(busco_check, filenames)
        
print("Number of replacements:",step)


