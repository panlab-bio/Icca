#!/usr/bin/python
import os
import re
import glob
import pandas as pd
import numpy as np
from collections import defaultdict
import sys, getopt
import os.path
import subprocess
from Bio import SeqIO
import multiprocessing

os.system("ln -s ../02_locate_gene/all_con_filter.bed ./")
os.system("ln -s ../01_pre/sca_reads_aln.paf ./sca_reads_aln.paf")
os.system("ln -s ../01_pre/all_con_reads_aln.paf ./all_con_reads_aln.paf")
os.system("ln -s ../01_pre/*_reads.fa ./")
os.system("ln -s ../01_pre/*_odb10 ./")

####### get_recall_reads ######################
with open("all_con_filter.bed","r")as f:
    pwd=os.getcwd()
    for line in f:
        
        fields = line.strip().split("\t")
        key    = fields[8]
        
        chr_name = fields[0]  #sca
        strat1 = fields[2]
        end1   = fields[3]
        
        con_name = fields[4]  #con
        strat2 = fields[6]
        end2   = fields[7]
        
        #print("step:",step,"con_name:",con_name,"chr_name:",chr_name)
        
        file_name=str(key)+'_test'
        os.makedirs(file_name, exist_ok=True)
        subdir_path = os.path.join(pwd, file_name)
        os.chdir(subdir_path)
        
        sca_file_path = os.path.join(pwd, "sca_reads_aln.paf")
        with open(sca_file_path,"r") as f_sca:
            
            for line in f_sca:
                fields = line.strip().split("\t")
                query_name = fields[5]  #sca
                sca_start = fields[7]
                sca_end   = fields[8]
                #print("query_name:",query_name)
                if query_name == chr_name:
                    if ((sca_end >= strat1 and sca_end <= end1) or (sca_start >= strat1 and sca_end <= end1) or (sca_start >= strat1 and sca_start <= end1)):
                        
                        with open(str(key) + "_sca_reads.paf", "a") as out:
                             out.write(line)
        
        con_file_path = os.path.join(pwd, "all_con_reads_aln.paf")
        with open(con_file_path,"r") as f_con:
            
            for line in f_con:
                fields = line.strip().split("\t")
                query_name = fields[5]  #con
                con_start = fields[7]
                con_end   = fields[8]
                #print("query_name:",query_name)
                
                if query_name == con_name:
                    #if ((con_end >= strat1 and con_end <= end1) or (con_start >= strat1 and con_end <= end1) or (con_start >= strat1 and con_start <= end1)):
                       
                        with open(str(key) + "_con_reads.paf", "a") as out:
                             out.write(line)   
        
        file1 = str(key) + "_sca_reads.paf"
        file2 = str(key) + "_con_reads.paf"
        
       
        subprocess.run(f"awk '{{print $1}}' {file1} | sort | uniq > su_file1", shell=True)

        
        subprocess.run(f"awk '{{print $1}}' {file2} | sort | uniq > su_file2", shell=True)
        
        
        file4 = str(key) + "_sca_con_bing_reads.txt"
        subprocess.run(f"sort su_file1 su_file2 | uniq > {file4}", shell=True)
        os.chdir(pwd)    

####### recall_hifi_reads ######################
for folder in os.listdir(pwd):
    if os.path.isdir(folder) and folder.endswith("_test"):
        os.chdir(os.path.join(pwd, folder))
        
        sca_con_reads_file = None
        for file_name in os.listdir():
            if file_name.endswith("_sca_con_bing_reads.txt"):
                sca_con_reads_file = file_name
                break
        
        if sca_con_reads_file:

            cmd = f"seqtk subseq ../*_reads.fa {sca_con_reads_file} > tq_reads.fa"
            os.system(cmd)
        
        os.chdir(pwd)
        
def run_hicanu(folder):
    
    print("run_hicanu,folder:",folder)
    os.chdir(folder)
    cmd = f"canu -p asm genomeSize=20m -pacbio-hifi ./tq_reads.fa"
    os.system(cmd)
    os.chdir("..")
    print(f"{folder} hicanu done")
    
def run_hifiasm(folder):
   
    print("run_hifiasm,folder:",folder)
    try:
        os.chdir(folder)
        cmd = f"hifiasm -t 32 -o tq ./tq_reads.fa"
        os.system(cmd)
        os.chdir("..")
        print(f"{folder} hifiasm done")
    except Exception as e:
        print(f"run_hifiasm error: {e}")

def run_busco_hifiasm(folder):
    
    os.chdir(folder)
    cmd = f"busco -i ./tq.bp.p_ctg.fa -c 60 -o hifiasm_tq_busco -m geno -l ../*_odb10 --offline"
    os.system(cmd)
    os.chdir("..")
    print(f"{folder} run_busco_hifiasm done")

def run_busco_hicanu(folder):
    
    os.chdir(folder)
    cmd = f"busco -i asm.contigs.fasta -c 60 -o hicanu_tq_busco -m geno -l ../*_odb10 --offline"
    os.system(cmd)
    os.chdir("..")
    print(f"{folder} run_busco_hicanu done")

def busco_check(dif):
    
    pwd = os.getcwd()
    
    for folder in os.listdir(pwd):
        if os.path.isdir(folder) and folder.endswith("_test"):
            folder_path = os.path.join(pwd, folder)
            os.chdir(folder_path)
            folder_name = folder
            prefix = folder_name.split("_")[0]
            
            if dif==1:
                #hifiasm_tq_busco
                tq_busco_path = glob.glob(os.path.join(folder_path, "hifiasm_tq_busco"))[0] 
            if dif==2:
                # #hicanu_tq_busco
                tq_busco_path = glob.glob(os.path.join(folder_path, "hicanu_tq_busco"))[0]  
            
            print("tq_busco_path:",tq_busco_path)
            
            
            if os.path.exists(tq_busco_path):
                
                for run_folder in os.listdir(tq_busco_path):
                    if os.path.isdir(os.path.join(tq_busco_path, run_folder)) and run_folder.startswith("run_"):
                        file_path = os.path.join(tq_busco_path, run_folder, "full_table.tsv")
                      
                        if os.path.exists(file_path):
                            print("file_path:",file_path)
                            print("The file exists.")
                                            
                        with open(file_path, "r") as f:
                            for line in f:
                                if prefix in line and "Complete" in line:
                                    columns = line.split("\t")
                                    out_file = f"{prefix}_success.bed"
                                    with open(out_file, "w") as f:
                                        f.write(line)
                                    print(f"The matching lines have been written to the file {out_file}.")
                                  
                                    seq_name = columns[2]
        
                                    print("folder_path:",folder_path)
                                    
                                    if dif==1:
                                        fasta_file = os.path.join(folder_path,"tq.bp.p_ctg.fa")
                                    #tq.bp.p_ctg.fa
                                    if dif==2:
                                        fasta_file = os.path.join(folder_path,"asm.contigs.fasta")
                                        
                                    
                                    with open(fasta_file, "r") as f:
                                        records = SeqIO.parse(f, "fasta")
                                        for record in records:
                                            if record.id == seq_name:
                                                seq_record = record
                                                break
                                    out_file = f"{prefix}_{seq_name}_sequence.fasta"
                                    with open(out_file, "w") as f:
                                        SeqIO.write(seq_record, f, "fasta")
                                    n_out_file = f"{prefix}_{seq_name}_single_sequence.fasta"
                                    convert_fasta_to_single_line(out_file, n_out_file)
                                    print(f"The sequence {seq_name} has been written to the file {n_out_file}.")
                                    print("\n")                     

            os.chdir(pwd)

def convert_fasta_to_single_line(infile, outfile):
    with open(infile, "r") as f, open(outfile, "w") as out_f:
        for record in SeqIO.parse(f, "fasta"):
            sequence = str(record.seq).replace("\n", "")
            out_f.write(f">{record.id}\n{sequence}\n")

def is_file_empty(file_path):
    return os.path.isfile(file_path) and os.path.getsize(file_path) == 0

pwd = os.getcwd()
folders = [os.path.join(pwd, folder) for folder in os.listdir(pwd) if os.path.isdir(folder) and folder.endswith("_test")]

hifiasm_flag=1
with multiprocessing.Pool(processes=4) as pool:
    pool.map(run_hifiasm, folders)
    
for folder in folders:
    folder_path = os.path.join(pwd, folder)
    os.chdir(folder_path)

    if not os.path.exists("tq.bp.p_ctg.gfa"):
        hifiasm_flag = 0
        break

    os.chdir(pwd)

if hifiasm_flag==1:
    pwd = os.getcwd()
    for i, folder in enumerate(os.listdir(pwd)):
        if os.path.isdir(os.path.join(pwd, folder)) and folder.endswith("_test"):
            os.chdir(os.path.join(pwd, folder))
            print(os.getcwd()) 
            cmd=f"awk \'/^S/{{print \">\"$2;print $3}}\' tq.bp.p_ctg.gfa > tq.bp.p_ctg.fa"
            os.system(cmd)
            os.system("chmod +x tq.bp.p_ctg.fa")
            os.chdir(pwd)
    print("hifiasm_2 done")

    with multiprocessing.Pool(processes=4) as pool:
        pool.map(run_busco_hifiasm, folders)

    dif=1
    busco_check(dif)
    
os.chdir(pwd) 
print("pwd:",os.getcwd()) 

os.system("cat *_test/*_success.bed > all_success.bed")

file_path = 'all_success.bed'

if is_file_empty(file_path):
    
    with multiprocessing.Pool(processes=4) as pool:
        pool.map(run_hicanu, folders)

    with multiprocessing.Pool(processes=4) as pool:
        pool.map(run_busco_hicanu, folders)
    dif=2
    busco_check(dif)

os.system("cat *_test/*_success.bed > all_success.bed")

    
