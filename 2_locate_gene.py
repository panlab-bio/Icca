#!/usr/bin/python
import os
import glob
import pandas as pd
import numpy as np
from collections import defaultdict
import sys, getopt
import os.path

os.system("ln -s ../01_pre/sca_con_mis.txt ./")
mis_gene_set=set()
with open("sca_con_mis.txt","r") as f:
    for line in f.readlines():                          
           line = line.strip()
           mis_gene_set.add(line)
print("sca_con_misgene_save_finish\n")
print("mis_gene_set:",mis_gene_set)

os.system("ln -s ../01_pre/all_con_full_table.tsv ./")

filename = "all_con_full_table.tsv"
new_filename = "locate_gene.info"

for i in mis_gene_set:
    keyword = i
    command = f"grep '\\b{keyword}\\b' {filename} > {keyword}_{filename}"
    os.system(command)
    if os.path.isfile(f"{keyword}_{filename}"):
        os.chdir(os.path.dirname(os.path.abspath(f"{keyword}_{filename}")))
        os.system(f"awk '{{print $1,$3,$4,$5}}' {keyword}_{filename} > {keyword}_{new_filename}")
        os.remove(f"{keyword}_{filename}")
    else:
        print(f"File {keyword}_{filename} not found.")
   
os.system("ln -s ../01_pre/all_con_sca_aln.paf ./")
data = {}
folder_path = "./"

for file_name in os.listdir(folder_path):
    if file_name.endswith("_locate_gene.info"):
        node_name = file_name.split("_")[0]
        data[node_name] = {}
        with open(os.path.join(folder_path, file_name), "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(node_name):
                    contig_name, start_pos, end_pos = line.split()[1:]
                    if contig_name not in data[node_name]:
                        data[node_name][contig_name] = []
                    data[node_name][contig_name].append((int(start_pos), int(end_pos)))

print(data)

matches = []
for key in data.keys():
    pwd = os.getcwd()
    folder_name = key + '_test'
    os.makedirs(folder_name, exist_ok=True)
    os.chdir(os.path.join(pwd, folder_name))
    con_id = list(data[key].keys())
    
    for con in con_id:
        flag=0
        value_list=list(data[key][con][0])
        gene_s = int(value_list[0])
        gene_e = int(value_list[1])
        
        with open("../all_con_sca_aln.paf","r") as f:
            for line in f:
                fields = line.strip().split("\t")
                query_name = fields[5]  #con
                target_name = fields[0] #chr
                map_s = int(fields[7])
                map_e = int(fields[8])
                
                if query_name == con and (gene_s >= map_s and gene_e <= map_e):
                    flag=1
                    with open(con + "_con.paf", "a") as out:
                        out.write(line)
                    with open("choose_gene.info", "a") as out2:
                        out2.write(f'{key}\t{con}\t{gene_s}\t{gene_e}\n')
                if(flag==1):
                    break
        
        for filename in glob.glob('*_con.paf'):
            with open(filename, 'r') as f_in, open(filename.replace('_con.paf', '_filter.paf'), 'a') as f_out:
                first_line = f_in.readline().strip().split()
                
                chr_name, chr_length = first_line[0], int(first_line[1])
                chr_start, chr_end = int(first_line[2]), int(first_line[3])
                
                contig_name, contig_length = first_line[5], int(first_line[6])
                contig_start, contig_end = int(first_line[7]), int(first_line[8])
                
                f_out.write(f'{chr_name}\t{chr_length}\t{chr_start}\t{chr_end}\t{contig_name}\t{contig_length}\t{contig_start}\t{contig_end}\t{key}\n')
        
        if(flag==1):
            break
        
    os.chdir(pwd)
    os.system("cat *_test/*_filter.paf > all_con_filter.bed")
    os.system("cat *_test/choose_gene.info > all_locate_gene.info")
