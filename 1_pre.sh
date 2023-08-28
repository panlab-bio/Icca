#You need to self-label your data for easier subsequent operations.
name=tair
#sca_dir
sca_dir=
ln -s $sca_dir/sca* ./${name}_sca.fa
#con_dir
con_dir=
cat $con_dir/*.fasta > ${name}_all_con.fa
#reads_dir
reads_dir=
ln -s $reads_dir/sca* ./${name}_reads.fa

##########pre_busco###########
#busco
odb_name=embryophyta_odb10
bdir=/data/liujunyang/31_merge_scaffold
ln -s $bdir/${odb_name} ./
#ln -s /data/liujunyang/31_merge_scaffold/embryophyta_odb10/ ./

busco -i ${name}_all_con.fa -c 60 -o all_con_busco -m geno -l embryophyta_odb10 --offline
busco -i ${name}_sca.fa -c 60 -o sca_busco -m geno -l embryophyta_odb10 --offline

ln -s all_con_busco/run_${odb_name}/missing_busco_list.tsv ./all_con_missing_busco_list.tsv
ln -s all_con_busco/run_${odb_name}/full_table.tsv ./all_con_full_table.tsv
ln -s sca_busco/run_${odb_name}/missing_busco_list.tsv ./sca_missing_busco_list.tsv

#busco_cha_mis
sed 1,3d sca_missing_busco_list.tsv > sca_missing_busco_list.txt
sed 1,3d all_con_missing_busco_list.tsv > all_con_missing_busco_list.txt
#Difference set
sort sca_missing_busco_list.txt all_con_missing_busco_list.txt all_con_missing_busco_list.txt | uniq -u > sca_con_mis.txt
chmod +x sca_con_mis.txt

##########pre_minimap2###########
#minimap2 -c ref.fa query.fq > alignment.paf
#con_sac_map
minimap2 -x asm5 -c ${name}_all_con.fa ${name}_sca.fa > all_con_sca_aln.paf
minimap2 -x asm5 -c ${name}_all_con.fa ${name}_reads.fa > all_con_reads_aln.paf
minimap2 -x asm5 -c ${name}_sca.fa ${name}_reads.fa > sca_reads_aln.paf
