#!/bin/bash
# by Yuying Chen on 2022/01/08
# Arctic metagenomic analysis


cd /mnt/nfs/software

#wget -c http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
#unzip fastqc_v0.11.5.zip
#cd FastQC
#chmod +x fastqc

cd /mnt/nfs/chenyy/TP_Arctic

# Trimmomatic对raw reads的修剪和过滤
mkdir clean_data
#source activate /home/zhangzh/anaconda3/envs/trimmomatic

#for samplename in $(cat samplename.txt)
#do 
#trimmomatic PE -phred33 raw_data/${samplename}_1.fq.gz raw_data/${samplename}_2.fq.gz clean_data/${samplename}_trimmed_R1.fq.gz clean_data/${samplename}_trimmed_U_R1.fq.gz clean_data/${samplename}_trimmed_R2.fq.gz clean_data/${samplename}_trimmed_U_R2.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15  MINLEN:36 -threads 20
#done &
#conda deactivate

#fastp对raw reads的修剪和过滤
#conda activate megahit
#conda install -c bioconda fastp

#for samplename in $(cat samplename.txt)
#do 
#fastp -i raw_data/${samplename}_1.fq.gz -I raw_data/${samplename}_2.fq.gz -o clean_data/${samplename}_trimmed_R1.fq.gz -O clean_data/${samplename}_trimmed_R2.fq.gz
#done &


cd clean_data

#
nohup /mnt/nfs/software/FastQC/fastqc *.fq.gz -t 80 &


for samplename in `ls *_trimmed_R1.fq.gz`
do
    samplename=${samplename/_trimmed_R1.fq.gz/}
    echo $samplename
done > ../samplename.txt
cd ..

mkdir clean_data_gunzip
cp clean_data/* clean_data_gunzip/
gunzip *



#MEGAHIT one sample assembly
mkdir magehit_assembly
conda activate megahit
for samplename in $(cat samplename.txt)
do
    R1=clean_data/${samplename}_trimmed_R1.fq.gz
    R2=clean_data/${samplename}_trimmed_R2.fq.gz
    echo "megahit -1 $R1 -2 $R2 --k-min 21 --k-max 141 --k-step 12 --out-dir magehit_assembly/${samplename} --out-prefix ${samplename} -t 40"
done  >generate_megahit.sh

ln -s /mnt/nfs/chenyy/TP_Arctic/magehit_assembly/*.fa /mnt/nfs/chenyy/TP_Arctic/contig/
#change file name 

#MEGAHIT co sample assembly

ls clean_data/A_cryoconite*_R1.fq.gz > Acryo_R1.txt
ls clean_data/A_cryoconite*_R2.fq.gz > Acryo_R2.txt
ls clean_data/A_Ice*_R1.fq.gz > Aice_R1.txt
ls clean_data/A_Ice*_R2.fq.gz > Aice_R2.txt
ls clean_data/T_*C*_R1.fq.gz > Tcryo_R1.txt
ls clean_data/T_*C*_R2.fq.gz > Tcryo_R2.txt
ls clean_data/T_*I*_R1.fq.gz > Tice_R1.txt
ls clean_data/T_*I*_R2.fq.gz > Tice_R2.txt


for i in *_R1.txt
do 
  prefix=$(basename $i _R1.txt)
  j=${prefix}_R2.txt
  R1=$(cat "$i"| cut -f3 -d,)
  R1_All=$(echo $R1|sed 's/ /,/g')
  R2=$(cat "$j"| cut -f3 -d,)
  R2_All=$(echo $R2|sed 's/ /,/g')
  megahit -1 "${R1_All}" -2 "${R2_All}" --k-min 21 --k-max 141 --k-step 10 --mem-flag 2 --out-dir ${prefix}_Coassembly --out-prefix ${prefix}_Coassembly -t 60
done

conda deactivate

ln /mnt/nfs/chenyy/TP_Arctic/*Coassembly/*.fa /mnt/nfs/chenyy/TP_Arctic/contig/

#Quast report
conda activate metawrap-env
cd contig
mkdir contig_quality
ls *Coassembly*.fa > ../Coassembly_contigs.txt

cd ..
python /mnt/nfs/software/quast-5.0.2/quast.py -t 30 -o contig/contig_quality/ contig/*.fa &

#remove short reads (#quality of contigs with 1000bp are great, default for metabat2, we can use 1000bp)
mkdir Contig_removed_short

for samplename in $(cat samplename.txt)
do
pullseq -i contig/${samplename}.contigs.fa -m 1000 > Contig_removed_short/${samplename}_contigs_1kb.fa
done &


for samplename in $(cat Coassembly_contigs.txt)
do
pullseq -i contig/${samplename}.contigs.fa -m 1000 > Contig_removed_short/${samplename}_contigs_1kb.fa
done &



conda deactivate


#Non-reduency gene catloge
mkdir NR

bash Multigenome_deredundancy.sh Contig_removed_short fa NR
cd NR
grep ">" Merged_GENE_deredundancy_rep_protein_seq.fasta > rep_protein.id
bash Extract_fasta_sequence_by_ID.sh Merged_genome_prodigal.fna rep_protein.id Merged_GENE_deredundancy_rep_seq


####Salmon定量####
conda activate metawrap-env
mkdir gene_salmon
cd gene_salmon
#1.使用salmon 对基因集构建index
nohup salmon index -p 80 -t ../NR/Merged_GENE_deredundancy_rep_seq.fasta -i Merged_GENE_deredundancy_rep_seq.fasta_index &
#2.对样本比对定量：
cd ..
for samp in $(cat samplename.txt)
do 
salmon quant -i gene_salmon/Merged_GENE_deredundancy_rep_seq.fasta_index -l A  -1 clean_data/${samp}_trimmed_R1.fq.gz   -2 clean_data/${samp}_trimmed_R2.fq.gz  -p 30 --validateMappings -o gene_salmon/${samp}_quant >gene_salmon/${samp}.salmon.log 2>&1; 
done &


cd gene_salmon
find . -name quant.sf
chmod +x gather-counts.py
./gather-counts.py


for file in *counts
do
  # 提取样品名
  name=${file%%.*}
  echo $name
  # 将每个文件中的count列改为样品列
  sed -e "s/count/$name/g" $file > tmp
  mv tmp $file
done
# 合并所有样品
paste *counts |cut -f 1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48 > Combined-counts-TPM.txt



awk ' {print $1,$2+$3+$4,$5+$6+$7,$8+$9+$10+$12+$13+$14+$18+$19+$20,$11+$15+$16+$17+$21+$22+$23+$24}' Combined-counts-TPM.txt > Combined-Sum.txt
awk ' {print $1,$2}' Combined-Sum.txt > Combined-Sum-Arctic_Cryoconite.txt
awk ' {print $1,$4}' Combined-Sum.txt > Combined-Sum-TP_Cryoconite.txt
awk ' {print $1,$3}' Combined-Sum.txt > Combined-Sum-Arctic-Ice.txt
awk ' {print $1,$5}' Combined-Sum.txt > Combined-Sum-TP-Ice.txt
#Calculate the shred between two column
grep -F -f a.txt b.txt | sort | uniq
comm a.file b.file  
a(find(sum(a,2)==0),:)=[]%

awk '$2!=0' Combined-Sum-Arctic_Cryoconite.txt > Combined-Sum-Arctic_Cryoconite_remove.txt
awk '$2!=0' Combined-Sum-TP_Cryoconite.txt > Combined-Sum-TP_Cryoconite_remove.txt
awk '$2!=0' Combined-Sum-Arctic-Ice.txt > Combined-Sum-Arctic-Ice_remove.txt
awk '$2!=0' Combined-Sum-TP-Ice.txt > Combined-Sum-TP-Ice_remove.txt

awk ' {print $1}' Combined-Sum-Arctic_Cryoconite_remove.txt > Combined-Sum-Arctic_Cryoconite_remove_name.txt
awk ' {print $1}' Combined-Sum-TP_Cryoconite_remove.txt > Combined-Sum-TP_Cryoconite_remove_name.txt
awk ' {print $1}' Combined-Sum-Arctic-Ice_remove.txt > Combined-Sum-Arctic-Ice_remove_name.txt
awk ' {print $1}' Combined-Sum-TP-Ice_remove.txt > Combined-Sum-TP-Ice_remove_name.txt

comm -12 <(sort Combined-Sum-Arctic_Cryoconite_remove_name.txt|uniq) <(sort Combined-Sum-TP_Cryoconite_remove_name.txt|uniq) > Cryoconite_comm_test.txt
wc -l Cryoconite_comm_test.txt
comm -12 <(sort Combined-Sum-Arctic-Ice_remove_name.txt|uniq) <(sort Combined-Sum-TP-Ice_remove_name.txt|uniq) > Ice_comm_test.txt
wc -l Ice_comm_test.txt

#grep -Ff Combined-Sum-Arctic_Cryoconite_remove_name.txt Combined-Sum-TP_Cryoconite_remove_name.txt > Cryoconite_comm.txt
#grep -Ff Combined-Sum-Arctic-Ice_remove_name.txt Combined-Sum-TP-Ice_remove_name.txt > Ice_comm.txt

wc -l Combined-Sum-Arctic_Cryoconite_remove_name.txt Combined-Sum-TP_Cryoconite_remove_name.txt
wc -l Combined-Sum-Arctic-Ice_remove_name.txt Combined-Sum-TP-Ice_remove_name.txt



#R

data<-read.table('Combined-counts-TPM.tab',header=T,row.names=1)
data<-data[which(rowSums(data) > 0),]
write.csv("Combined-counts-TPM.csv")
data_R<-data/1000000

library(vegan)
library(ggplot2)

otu<-read.table("Combined-counts-TPM.tab",header = T,row.names = 1,sep = "\t")
colSums(otu != 0)
otu_t<-data.frame(t(otu))
group<-read.csv("Group.csv",header = T,row.names = 1)

vare.mds <- metaMDS(otu_t)
vare.mds 

data.scores <- as.data.frame(scores(vare.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$name <- rownames(data.scores)  # create a column of species, from the rownames of species.scores
data.scores$Region <- group$Region  #  add the grp variable created earlier
data.scores$Habitat <- group$Habitat
head(data.scores)  #look at the data



ggplot() + 
  # geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=name),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=Habitat,colour=Region),size=5) + # add the point markers
  # geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=name),size=5,vjust=0) +  # add the site labels
  annotate("text",x=0.3,y=0.3,label=paste("Stress=",format(vare.mds$stress, digits=4), sep=""))+
  xlim(-0.5,0.5)+ylim(-0.5,0.5)+
  labs(title = "16S")+
  coord_equal() +
  theme_bw()+
  theme(axis.text.x = element_text(size=16),  # remove x-axis text
        axis.text.y = element_text(size=16), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        plot.title = element_text(size=20,hjust = 0.5),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())


cd ..
mkdir gene_function_annotation


#KEGG
cd gene_function_annotation
mkdir KEGG
cd ..
/home/chenyy/Soft/kofamscan/bin/kofam_scan-1.3.0/exec_annotation --cpu 80 -f detail-tsv -o gene_function_annotation/KEGG/NR_kegg_raw.txt NR/Merged_GENE_deredundancy_rep_seq.fasta &> kegg.log
cd ../KEGG

grep '*' NR_kegg_raw.txt > kegg_sort.txt
awk '{print $2,$3}' kegg_sort.txt > kegg_sort_gene_ko_name.txt
awk -F, 'NR==1 {print "transcript","KO"} 
                 {gsub(/"/,""); 
                  print $1,$2}' kegg_sort_gene_ko_name.txt | column -t > kegg_sort_gene_ko_name_new.txt

cat ../../gene_salmon/Combined-counts-TPM.tab | grep -f kegg_sort_gene_name.txt > Marker_gene_all_gene_salmon.txt



R
library(vegan)
library(ggplot2)


kegg_test<-read.csv("kegg_sort_gene_ko_salmon.csv",header=T,row.names=1)
colSums(kegg_test != 0)
otu<-kegg_test[,-c(1:2)]

otu_t<-data.frame(t(otu))
group<-read.csv("Group.csv",header = T,row.names = 1)

vare.mds <- metaMDS(otu_t)
vare.mds 

data.scores <- as.data.frame(scores(vare.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$name <- rownames(data.scores)  # create a column of species, from the rownames of species.scores
data.scores$Region <- group$Region  #  add the grp variable created earlier
data.scores$Habitat <- group$Habitat
head(data.scores)  #look at the data



ggplot() + 
  # geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=name),alpha=0.5) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=Habitat,colour=Region),size=5) + # add the point markers
  # geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=name),size=5,vjust=0) +  # add the site labels
  annotate("text",x=0.3,y=0.3,label=paste("Stress=",format(vare.mds$stress, digits=4), sep=""))+
  xlim(-0.5,0.5)+ylim(-0.5,0.5)+
  labs(title = "16S")+
  coord_equal() +
  theme_bw()+
  theme(axis.text.x = element_text(size=16),  # remove x-axis text
        axis.text.y = element_text(size=16), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        plot.title = element_text(size=20,hjust = 0.5),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
		
data<-read.table("/mnt/nfs/chenyy/TP_Arctic/gene_salmon/Combined-counts-TPM.tab",header=T)
rs<-read.table("kegg_sort_gene_ko_name_new.txt",header=T)


df <- merge(rs,data,by="transcript")


kegg<-read.table("KEGG_pathway_ko.txt",header=T, na.strings = "",sep="\t",quote = "", fill=TRUE)
kegg_path<-merge(kegg,df,by="KO")

kegg_aggP1<-aggregate(kegg_path[,12:34],by=list(kegg_path$P1),sum)
kegg_aggP2<-aggregate(kegg_path[,12:34],by=list(kegg_path$P2),sum)
kegg_aggP3<-aggregate(kegg_path[,12:34],by=list(kegg_path$P3),sum)
kegg_aggKO<-aggregate(kegg_path[,12:34],by=list(kegg_path$KO),sum)

write.csv(df,"kegg_sort_gene_ko_salmon.csv")
write.csv(kegg_path,"kegg_pathway_salmon.csv")
write.csv(kegg_aggP1,"kegg_pathway_P1_salmon_sum.csv")
write.csv(kegg_aggP2,"kegg_pathway_P2_salmon_sum.csv")
write.csv(kegg_aggP3,"kegg_pathway_P3_salmon_sum.csv")
write.csv(kegg_aggKO,"kegg_pathway_KO_salmon_sum.csv")




cat kegg_sort.txt | grep -f Marker_gene_83.txt > Marker_gene_from_dataset.txt
awk '{print $2}' Marker_gene_from_dataset.txt > Marker_gene_gene_name.txt
cat ../../gene_salmon/Combined-counts-TPM.tab | grep -f Marker_gene_gene_name.txt > Marker_gene_from_dataset_salmon.txt






ls Contig_removed_short/* >contig_removed_short.txt



#######################################################################
########################################################################
######################################################################
###################################NR_FUNCTION_BLAST######################
conda activate metawrap-env
ln -s /mnt/nfs/chenyy/cryoconite/Reference_database/Greening_lab/* /mnt/nfs/chenyy/TP_Arctic/Reference/
cd /mnt/nfs/chenyy/TP_Arctic/Reference

for file in *.dmnd
do name=$(basename $file .dmnd)
   echo $name
   diamond blastp -d ${name} -q /mnt/nfs/chenyy/TP_Arctic/NR/Merged_GENE_deredundancy_rep_protein_seq.fasta -o /mnt/nfs/chenyy/TP_Arctic/Diamond_NR_blast_out/${name}_out -p 60 -f 6 --query-cover 80
done

cd /mnt/nfs/chenyy/TP_Arctic/Diamond_NR_blast_out
#identity>50,只保留第一次出现的字段
for file in *_out
do
  name=$(basename $file _out)
  cat $file | awk '$3>50{print $0}' >${name}_identity
done


for file in *_identity
do
  name=$(basename $file _identity)
  time awk 'BEGIN{FS="|"}{print $0}' $file|awk '!x[$1]++' >${name}_uniq
done

cat FeFe_hydrogenase_uniq | awk '$3>60{print $0}' >FeFe_hydrogenase_uniq1
cat Carbon_monoxide_dehydrogenase_CoxL_uniq | awk '$3>60{print $0}' >Carbon_monoxide_dehydrogenase_CoxL_uniq1

cat Nitrite_oxidoreductase_NxrA_uniq | awk '$3>60{print $0}' >Nitrite_oxidoreductase_NxrA_uniq1
cat NADH-ubiquinone_oxidoreductase_NuoF_uniq | awk '$3>60{print $0}' >NADH-ubiquinone_oxidoreductase_NuoF_uniq1
cat F-type_ATP_synthase_AtpA_uniq | awk '$3>70{print $0}' >F-type_ATP_synthase_AtpA_uniq1
cat Selenate_reductase_YgfK_uniq | awk '$3>70{print $0}' >Selenate_reductase_YgfK_uniq1
cat Thaumarchaeota_4-hydroxybutyryl-CoA_synthetase_HbsT_sequences_uniq | awk '$3>70{print $0}' >Thaumarchaeota_4-hydroxybutyryl-CoA_synthetase_HbsT_sequences_uniq1
cat Arsenite_oxidase_ARO_uniq | awk '$3>70{print $0}' >Arsenite_oxidase_ARO_uniq1
cat Photosystem_II_reaction_centre_PsbA_sequences_uniq | awk '$3>70{print $0}' >Photosystem_II_reaction_centre_PsbA_sequences_uniq1
cat Photosystem_I_reaction_centre_PsaA_sequences_uniq | awk '$3>80{print $0}' >Photosystem_I_reaction_centre_PsaA_sequences_uniq1


cat *_uniq > All_marker_blast

awk -F, 'NR==1 {print "transcript","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"} 
                 {gsub(/"/,""); 
                  print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' All_marker_blast  | column -t > All_marker_blast.txt

R
data1<-read.table("All_marker_blast.txt",header=T,sep="")
data2<-read.table("markergene_name.txt",sep="",header=T)
df1 <- merge(data1,data2,by="sseqid")
write.csv(df1,"All_marker_blast_name.csv",row.names=FALSE)


data3<-read.table("/mnt/nfs/chenyy/TP_Arctic/gene_salmon/Combined-counts-TPM.txt",header=T)
colnames(data3)[1] <- 'transcript'
df2 <- merge(df1,data3,by="transcript")
write.csv(df2,"All_marker_blast_salmon.csv",row.names=FALSE)

df3<-aggregate(df2[,16:38],list(df2[,14]),sum)
write.csv(df3,"All_marker_blast_Name_salmon_aggre.csv",row.names=FALSE)

q()






#################################################################
#################################################################
#################################################################


####Phyloflash####
mkdir phyloflash
conda activate phyloflash
cp samplename.txt phyloflash/
cd phyloflash
for samplename in $(cat samplename.txt)
do 
    echo  "time phyloFlash.pl -dbhome /mnt/nfs/software/Silva_138.1_PhyloFlash_database -CPUs 10 -lib $samplename -almosteverything -log -read1 ../clean_data/${samplename}_trimmed_R1.fq.gz -read2 ../clean_data/${samplename}_trimmed_R2.fq.gz"  
done >Generate_phyloflash.sh

nohup cat Generate_phyloflash.sh |parallel -j 10 &>phyloflash/Generate_phyloflash.log &

#rename 's/_trimmed_R1.fq.gz.SSU.1.fq/_1.fastq/' *
#rename 's/_trimmed_R1.fq.gz.SSU.2.fq/_2.fastq/' *


for phyloname in `ls *_trimmed_R1.fq.gz.SSU.1.fq`
do
    phyloname=${phyloname/_trimmed_R1.fq.gz.SSU.1.fq/}
    echo $phyloname
done > phyloname.txt
cd ..
conda deactivate



#PM#
mkdir PM
cd PM
mkdir Seqs
ln -s /mnt/nfs/chenyy/TP_Arctic/phyloflash/*.fq /mnt/nfs/chenyy/TP_Arctic/PM/Seqs

source deactivate
PM-pipeline -D S -m meta_pair.txt -i seqs_pair.list -o Silva_16S_1 -d 0.97 -L 123456 -F 1234 &>log&
PM-pipeline -D E -m meta_pair.txt -i seqs_pair.list -o Silva_18s1 -d 0.97 -L 123456 -F 1234 &>log&



#########################################################
####################Binning##############################
#########################################################
conda activate metawrap-env
mkdir binning

for samplename in $(cat samplename.txt)
do 
    echo  "time metawrap binning -o binning/${samplename}_bin -t 40 -a Contig_removed_short/${samplename}_contigs_1kb.fa \
          --metabat2 --maxbin2 --universal clean_data_gunzip/${samplename}_trimmed_1.fastq clean_data_gunzip/${samplename}_trimmed_2.fastq &>binning/${samplename}_binning.log"  
done >metawrap_binning.sh


nohup cat metawrap_binning.sh |parallel -j 5 &> binning/metawrap_binning.log &

nohup time metawrap binning -o binning/Acryo_Coassembly_bin -t 50 -a Contig_removed_short/Acryo_Coassembly_contigs_1kb.fa \
          --metabat2 --maxbin2 --universal clean_data_gunzip/A_cryo*_trimmed_1.fastq clean_data_gunzip/A_cryo*_trimmed_2.fastq &>binning/Acryo_Coassembly_binning.log &

nohup time metawrap binning -o binning/Aice_Coassembly_bin -t 50 -a Contig_removed_short/Aice_Coassembly_contigs_1kb.fa \
          --metabat2 --maxbin2 --universal clean_data_gunzip/A_Ice*_trimmed_1.fastq clean_data_gunzip/A_Ice*_trimmed_2.fastq &>binning/Aice_Coassembly_binning.log &

nohup  time metawrap binning -o binning/Tcryo_Coassembly_bin -t 50 -a Contig_removed_short/Tcryo_Coassembly_contigs_1kb.fa \
          --metabat2 --maxbin2 --universal clean_data_gunzip/T_*C*_trimmed_1.fastq clean_data_gunzip/T_*C*_trimmed_2.fastq &>binning/Tcryo_Coassembly_binning.log &
  
nohup  time metawrap binning -o binning/Tice_Coassembly_bin -t 50 -a Contig_removed_short/Tice_Coassembly_contigs_1kb.fa \
          --metabat2 --maxbin2 --universal clean_data_gunzip/T_*I*_trimmed_1.fastq clean_data_gunzip/T_*I*_trimmed_2.fastq &>binning/Tice_Coassembly_binning.log &





#把文件名中的AA替换成aa
rename "s/bin/A_Ice1_meta/" *
rename "s/bin/A_Ice1_max/" *





#dRep去重复
mkdir All_MAGs
ln -s /mnt/nfs/chenyy/TP_Arctic/binning/*/*/*.fa /mnt/nfs/chenyy/TP_Arctic/All_MAGs/
mkdir dRep_MAGs
conda activate drep
dRep bonus output_directory --check_dependencies

cd dRep_MAGs
dRep dereplicate dereplicate_wf -p 40 -comp 50 -con 10 -g ../Greenlan/* &
#聚类OTU
dRep dereplicate dereplicate_otu -p 40 -comp 50 -con 10 -nc 0.3 -sa 0.95 -g All_178samples_MAGs/* &


conda deactivate

####GTDB注释####
mkdir GTDB
conda activate gtdbtk


gtdbtk classify_wf --genome_dir dereplicate_otu4/dereplicated_genomes/ \
    --out_dir GTDB/classify_wf \
    --extension fa \
    --prefix bin \
    --cpus 40 &


####Tree####
conda deactivate
cd /mnt/nfs/chenyy/TP_Arctic/Tree
muscle -in bin.user_msa.fasta -out bin.user_msa_algin.fasta
FastTree -gtr -gamma < bin.user_msa_algin.fasta > bin.user_msa_algin_tree.nwk



#CheckM
conda activate drep
checkm lineage_wf -t 100 -x fa --nt --tab_table -f CheckM/bins_qa.txt MAG CheckM
#checkm ssu_finder -t 20 -x fa ALL_DKMD_maxbin2_bin.0.fa bin ssu_finder







#构建进化树Tree
cd /mnt/nfs/chenyy/TP_Arctic/Tree
#从fasta序列提取特定序列
seqkit grep -f MAG_grep.txt bin.bac120.user_msa.fasta -o MAG_grep.fasta

muscle -in MAG_grep.fasta -out MAG_grep_algin.fasta
FastTree -gtr -gamma < MAG_grep_algin.fasta > MAG_grep_algin_tree.nwk


grep ">" MAG_grep_algin.fasta >  MAG_grep_algin_name.txt


 # 生成注释文件中每列为单独一个文件

## 方案2. 生成丰度柱形图注释文件
Rscript table2itol.R -a -d -c none -D plan2 -b Phylum -i OTUID -l Genus -t %s -w 0.5 annotation.txt

## 方案3. 生成热图注释文件
Rscript table2itol.R -c keep -D plan3 -i OTUID -t %s otutab.txt

## 方案4. 将整数转化成因子生成注释文件
Rscript table2itol.R -a -c factor -D plan4 -i OTUID -l Genus -t %s -w 0 annotation.txt

## 方案5. 自定义颜色
Rscript table2itol.R -a -C table2itol/tests/INPUT/colours_1.yml -c double -D plan5 \
  -i OTUID -l Genus -t %s -w 0.5 annotation.txt



####CoverM####
conda activate coverm
mkdir MAGs
mkdir MAGs_coverm
ln -s /mnt/nfs/chenyy/TP_Arctic/dRep_MAGs/dereplicate_wf1/dereplicated_genomes/* /mnt/nfs/chenyy/TP_Arctic/MAGs/

cd /mnt/nfs/chenyy/TP_Arctic/clean_data_gunzip
####
#for i in *_1.fastq ;do j=${i%_*}_2.fastq ; m=${i%_*}; n=${m#*/};echo "coverm genome --coupled $i $j -t 10 --methods mean --genome-fasta-files ../MAGs/*.fa -o ../MAGs_coverm/${n}_counts.tsv " ;done >coverm2.sh
#nohup cat coverm2.sh |parallel -j 10 &>all.log & 
#for i in *_1.fastq ;do j=${i%_*}_2.fastq ; m=${i%_*}; n=${m#*/};echo "coverm genome --coupled $i $j -t 10 --methods covered_fraction --genome-fasta-files ../MAGs/*.fa -o ../MAGs_coverm/${n}_counts.tsv " ;done >coverm2.sh
#for i in *_1.fastq ;do j=${i%_*}_2.fastq ; m=${i%_*}; n=${m#*/};echo "coverm genome --coupled $i $j -t 10 --methods rpkm --genome-fasta-files ../MAGs/*.fa -o ../MAGs_coverm/${n}_counts.tsv " ;done >coverm2.sh




####
for i in *_1.fastq ;do j=${i%_*}_2.fastq ; m=${i%_*}; n=${m#*/};echo "coverm genome --coupled $i $j -t 10 --genome-fasta-files ../MAGs/*.fa -o ../MAGs_coverm/${n}_coverm.tsv " ;done >coverm3.sh
nohup cat coverm3.sh |parallel -j 10 &>all.log & 

cd /mnt/nfs/chenyy/TP_Arctic/MAGs_coverm
paste *tsv |cut -f 1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44 > Combined-CoverM_relative.txt








####EnrichM####
conda activate enrichm_0.5.0

cd /mnt/nfs/chenyy/TP_Arctic/MAGs
rename "s/.fa/.fna/" *


cd ..
enrichm annotate --log annotate_LOG --output EnrichM/Output --force --genome_directory MAGs --ko_hmm --threads 10 --parallel 4 
enrichm classify --log classify_LOG --output EnrichM/classify_Output --force --genome_and_annotation_matrix Output/ko_hmm_frequency_table.tsv 



#prokka注释
for i in *.fa
do 
  prefix=$(basename $i .fa)
  nohup time prokka $i --prefix ${prefix} --metagenome -cpus 20 --kingdom Bacteria
done



#KEGG注释
/home/chenyy/Soft/kofamscan/bin/kofam_scan-1.3.0/exec_annotation --cpu 20 -f detail-tsv -o test_kegg_raw.txt test.fna &> kegg.log
cd /mnt/nfs/chenyy/cryoconite/dRep/dereplicate_wf/dereplicated_genomes
for i in *.fa;
do 
  prefix=$(basename $i .fa)
  /home/chenyy/Soft/kofamscan/bin/kofam_scan-1.3.0/exec_annotation --cpu 10 -f detail-tsv -o /mnt/nfs/chenyy/cryoconite/dRep/kegg/${prefix}_kegg_raw.txt /mnt/nfs/chenyy/cryoconite/dRep/dereplicate_wf/data/prodigal/${prefix}.fa.faa 
done &

cd /mnt/nfs/chenyy/cryoconite/dRep/kegg
for i in *.txt;
do 
  prefix=$(basename $i .txt)
  echo "grep '*' $i > result/${prefix}_result.txt & "
done >kegg_result.log
nohup cat kegg_result.log |parallel -j 20 &>All.log &
cd result/
for i in *.txt;
do 
  prefix=$(basename $i .txt)
  echo "cat $i | grep -f KO_number.txt > ../kegg_extrac/${prefix}_number.txt &"
done >kegg_extract.log
nohup cat kegg_extract.log |parallel -j 20 &>All.log &
cd ../kegg_extrac
for i in *.txt;
do 
  prefix=$(basename $i .txt)
  cat $i | awk -F '     ' '{print $3}'| sort | uniq -c | sort -nr > ../uniq/${prefix}_uni.txt
done &
cd ../uniq
for i in *.txt
do 
  prefix=$(basename $i .txt)
  awk '{print $3,$2,$1}' $i > ../vlookup/${prefix}_sort.txt
done
cd ../vlookup
for i in *_sort.txt
do 
  prefix=$(basename $i _sort.txt)
  awk 'NR==FNR{a[$1]=$0;next}NR>FNR{if($1 in a){print $0 "\t" a[$1]} else {print $0 "\t" "0"}}' $i KO_number.txt >final/${prefix}_vlookup.txt
  awk '{print $3}' final/${prefix}_vlookup.txt > result/${prefix}_resilt_final.txt
done
cd result/
ls | grep ".txt" >list
rm list
paste *.txt > merge.txt













####TGG_GGG_NR_blast####
cd /mnt/nfs/chenyy/cryoconite/Reference_database/Greening_lab
for file in *.dmnd
do name=$(basename $file .dmnd)
   echo $name
   diamond blastp -d ${name} -q /mnt/nfs/zhangzh/TGG_and_GGG/GGG/Gene/GGG_ORF/GGG_135Meta_and_1005ISO_prodiagl_protein_over500_dreped_rep_seq.fasta -o /mnt/nfs/chenyy/cryoconite/Diamond_GGG_blast_out/${name}_out -p 160 -f 6 --query-cover 80
done

cd /mnt/nfs/chenyy/cryoconite/Diamond_GGG_blast_out
#identity>50,只保留第一次出现的字段
for file in *_out
do
  name=$(basename $file _out)
  cat $file | awk '$3>50{print $0}' >${name}_identity
done

for file in *_identity
do
  name=$(basename $file _identity)
  time awk 'BEGIN{FS="|"}{print $0}' $file|awk '!x[$1]++' >${name}_uniq
done

cat FeFe_hydrogenase_uniq | awk '$3>60{print $0}' >FeFe_hydrogenase_uniq1
cat Carbon_monoxide_dehydrogenase_CoxL_uniq | awk '$3>60{print $0}' >Carbon_monoxide_dehydrogenase_CoxL_uniq1

cat Nitrite_oxidoreductase_NxrA_uniq | awk '$3>60{print $0}' >Nitrite_oxidoreductase_NxrA_uniq1
cat NADH-ubiquinone_oxidoreductase_NuoF_uniq | awk '$3>60{print $0}' >NADH-ubiquinone_oxidoreductase_NuoF_uniq1
cat F-type_ATP_synthase_AtpA_uniq | awk '$3>70{print $0}' >F-type_ATP_synthase_AtpA_uniq1
cat Selenate_reductase_YgfK_uniq | awk '$3>70{print $0}' >Selenate_reductase_YgfK_uniq1
cat Thaumarchaeota_4-hydroxybutyryl-CoA_synthetase_HbsT_sequences_uniq | awk '$3>70{print $0}' >Thaumarchaeota_4-hydroxybutyryl-CoA_synthetase_HbsT_sequences_uniq1
cat Arsenite_oxidase_ARO_uniq | awk '$3>70{print $0}' >Arsenite_oxidase_ARO_uniq1
cat Photosystem_II_reaction_centre_PsbA_sequences_uniq | awk '$3>70{print $0}' >Photosystem_II_reaction_centre_PsbA_sequences_uniq1
cat Photosystem_I_reaction_centre_PsaA_sequences_uniq | awk '$3>80{print $0}' >Photosystem_I_reaction_centre_PsaA_sequences_uniq1


cat *_uniq > All_marker_blast

awk -F, 'NR==1 {print "transcript","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"} 
                 {gsub(/"/,""); 
                  print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' All_marker_blast  | column -t > All_marker_blast.txt

R
data1<-read.table("All_marker_blast.txt",header=T,sep="")
data2<-read.table("markergene_name.txt",sep="",header=T)
df1 <- merge(data1,data2,by="sseqid")
write.csv(df1,"All_marker_blast_name.csv",row.names=FALSE)


data3<-read.table("/mnt/nfs/zhangzh/TGG_and_GGG/GGG/Gene/Mapping_salmon/result/133sample_tpm.tsv",header=T)
colnames(data3)[1] <- 'transcript'
df2 <- merge(df1,data3,by="transcript")
write.csv(df2,"All_marker_blast_salmon.csv",row.names=FALSE)

df3<-aggregate(df2[,16:148],list(df2[,14]),sum)
write.csv(df3,"All_marker_blast_Name_salmon_aggre.csv",row.names=FALSE)

q()







cd /mnt/nfs/chenyy/cryoconite/Diamond_blast_OTU_out
grep "Rubisco_RbcL_95_Filtered" *_diamond_out.txt.csv > Rubisco_RbcL_95_Filtered.txt

seqkit grep -f Rubisco_RbcL_95.txt *.faa -o Rubisco_RbcL_95.fasta





#dRep去重复
conda activate drep
dRep bonus output_directory --check_dependencies
#把文件名中的AA替换成aa
rename "s/bin/TGL1.bin/" *  
cd /mnt/nfs/chenyy/cryoconite/dRep
dRep dereplicate dereplicate_wf -p 30 -comp 50 -con 10 -g A_cryoconite1/* &
nohup dRep dereplicate dereplicate_wf -p 30 -comp 50 -con 10 -g A_cryoconite2/* &
nohup dRep dereplicate dereplicate_wf -p 30 -comp 50 -con 10 -g A_cryoconite3/* &

nohup dRep dereplicate dereplicate_wf -p 30 -comp 50 -con 10 -g A_Ice1/* &
nohup dRep dereplicate dereplicate_wf -p 30 -comp 50 -con 10 -g A_Ice2/* &
nohup dRep dereplicate dereplicate_wf -p 30 -comp 50 -con 10 -g A_Ice3/* &

nohup dRep dereplicate dereplicate_wf -p 30 -comp 50 -con 10 -g Acryo_Coassembly/* &
nohup dRep dereplicate dereplicate_wf -p 30 -comp 50 -con 10 -g Aice_Coassembly/* &





####OTU 功能注释
conda activate metawrap-env
mkdir MAGs_protein
cd /mnt/nfs/chenyy/TP_Arctic/All_data_MAG/dereplicate_otu4/dereplicated_genomes
for file in *.fa
do
  # 提取样品名
  name=$(basename $file .fa)
  ln -s /mnt/nfs/chenyy/TP_Arctic/All_data_MAG/dereplicate_otu4/data/prodigal/${name}.fa.faa /mnt/nfs/chenyy/TP_Arctic/All_data_MAG/MAGs_protein/
done



cd /mnt/nfs/chenyy/TP_Arctic/All_data_MAG/MAGs_protein
for file in *.faa
do
  # 提取样品名
  name=$(basename $file .faa)
  echo $name
  diamond blastp -d /mnt/nfs/chenyy/TP_Arctic/Reference/All/Greening_db -q $file -o /mnt/nfs/chenyy/TP_Arctic/Diamond_OTU_blast_out/${name}_out -p 30 -f 6 --query-cover 80
done



cd /mnt/nfs/chenyy/TP_Arctic/Diamond_OTU_blast_out
#identity>50,只保留第一次出现的字段
for file in *.fa_out
do
  name=$(basename $file .fa_out)
  cat $file | awk '$3>50{print $0}' >${name}_identity
done

for file in *_identity
do
  name=$(basename $file _identity)
  time awk 'BEGIN{FS="|"}{print $0}' $file|awk '!x[$1]++' >${name}_uniq
done


for file in *_uniq
do
  name=$(basename $file _uniq)
  awk -F, 'NR==1 {print "qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"} 
                 {gsub(/"/,""); 
                  print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' $file  | column -t > ${name}_diamond_out.txt
done


R
a = list.files(pattern = "_diamond_out.txt")
n = length(a)
 
for (i in 1:n){
   rs = read.table(a[i], header=T, sep="")
   data<-read.table("markergene_name.txt",sep="",header=T)
   df <-merge(rs,data,by="sseqid")
   filter(df, df$Name == "FeFe_hydrogenase" & df$pident >= 60)
   filter(df, df$Name == "Ammonia_monooxygenase_amoA" & df$pident >= 60)
   filter(df, df$Name == "Carbon_monoxide_dehydrogenase_CoxL" & df$pident >= 60)
   filter(df, df$Name == "Ammonia_monooxygenase_amoA" & df$pident >= 60)
   filter(df, df$Name == "Nitrite_oxidoreductase_NxrA" & df$pident >= 60)
   filter(df, df$Name == "NADH-ubiquinone_oxidoreductase_NuoF" & df$pident >= 60)
   filter(df, df$Name == "F-type_ATP_synthase_AtpA" & df$pident >= 70)
   filter(df, df$Name == "Selenate_reductase_YgfK" & df$pident >= 70)
   filter(df, df$Name == "Thaumarchaeota_4-hydroxybutyryl-CoA_synthetase_HbsT_sequences" & df$pident >= 70)
   filter(df, df$Name == "Arsenite_oxidase_ARO" & df$pident >= 70)
   filter(df, df$Name == "Photosystem_II_reaction_centre_PsbA_sequences" & df$pident >= 70)
   filter(df, df$Name == "Photosystem_I_reaction_centre_PsaA_sequences" & df$pident >= 80)
   write.csv(df,paste(a[i],".csv",sep=""),row.names=FALSE)
}


a = list.files(pattern = "_diamond_out.txt.csv")
n = length(a)  
for (i in 1:n){
   data = read.csv(a[i], header=T)
   aggreat<-table(data$Name)
   aggreat2<-table(data$Full)
   write.csv(aggreat,paste(a[i],"Name.csv",sep=""),row.names=FALSE)
   write.csv(aggreat2,paste(a[i],"Full.csv",sep=""),row.names=FALSE)
}


a = list.files(pattern = "_diamond_out.txt.csvName.csv")
n = length(a)  
for (i in 1:n){
   name=basename(a[i])
   data = read.csv(a[i], header=T)
   colnames(data)<-c("Name",name)
   write.csv(data,paste(a[i],"_rename.csv",sep=""),row.names=FALSE)
}

a = list.files(pattern = "_diamond_out.txt.csvFull.csv")
n = length(a)  
for (i in 1:n){
   name=basename(a[i])
   data = read.csv(a[i], header=T)
   colnames(data)<-c("Full",name)
   write.csv(data,paste(a[i],"_rename.csv",sep=""),row.names=FALSE)
}


a = list.files(pattern = "_diamond_out.txt.csvName.csv_rename.csv") 
n = length(a)                                    
merge.data = read.csv("Greening_lab_maker_name_merge.csv",header=T)  

for (i in 1:n){
   new.data = read.csv(a[i], header=T, sep=",")
   merge.data = merge(merge.data,new.data,by="Name",all=T)
}

write.csv(merge.data,file = "Maker_gene_OTU_Name_merge.csv",row.names=FALSE)



a = list.files(pattern = "_diamond_out.txt.csvFull.csv_rename.csv") 
n = length(a)                                    
merge.data = read.csv("Greening_lab_maker_Full_merge.csv",header=T)  

for (i in 1:n){
   new.data = read.csv(a[i], header=T, sep=",")
   merge.data = merge(merge.data,new.data,by="Full",all=T)
   
}

write.csv(merge.data,file = "Maker_gene_OTU_Full_merge.csv",row.names=FALSE)

q()


####基因富集分析
#install.packages('BiocManager')  #需要首先安装 BiocManager，如果尚未安装请先执行该步
BiocManager::install('DESeq2')

#（2）在 GitHub 中获取最新版本的方式
#install.packages('devtools')  #需要首先安装 devtools，如果尚未安装请先执行该步
devtools::install_github('mikelove/DESeq2@ae7c6bd')

#其它备注项：若中间提示有其它依赖 R 包的旧版包冲突的话，先删除旧包再安装新的
remove.packages('xxx')
BiocManager::install('xxx')


dat <- read.delim('control_treat.count.txt', row.names = 1, sep = '\t', check.names = FALSE)

#指定分组因子顺序
#注意要保证表达矩阵中的样本顺序和这里的分组顺序是一一对应的
coldata <- data.frame(condition = factor(rep(c('control', 'treat'), each = 3), levels = c('control', 'treat')))

library(DESeq2)

#第一步，构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata, design= ~condition)

#第二步，计算差异倍数并获得 p 值
#备注：parallel = TRUE 可以多线程运行，在数据量较大时建议开启
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)

#注意，需将 treat 在前，control 在后，意为 treat 相较于 control 中哪些基因上调/下调
res <- results(dds1, contrast = c('condition', 'treat', 'control'))

res

#输出表格至本地
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
write.table(res1, 'control_treat.DESeq2.txt', col.names = NA, sep = '\t', quote = FALSE)


res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

#log2FC≥1 & padj<0.01 标识 up，代表显著上调的基因
#log2FC≤-1 & padj<0.01 标识 down，代表显著下调的基因
#其余标识 none，代表非差异的基因
res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.01),'sig'] <- 'up'
res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.01),'sig'] <- 'down'
res1[which(abs(res1$log2FoldChange) <= 1 | res1$padj >= 0.01),'sig'] <- 'none'

#输出选择的差异基因总表
res1_select <- subset(res1, sig %in% c('up', 'down'))
write.table(res1_select, file = 'control_treat.DESeq2.select.txt', sep = '\t', col.names = NA, quote = FALSE)

#根据 up 和 down 分开输出
res1_up <- subset(res1, sig == 'up')
res1_down <- subset(res1, sig == 'down')

write.table(res1_up, file = 'control_treat.DESeq2.up.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(res1_down, file = 'control_treat.DESeq2.down.txt', sep = '\t', col.names = NA, quote = FALSE)

library(ggplot2)

#默认情况下，横轴展示 log2FoldChange，纵轴展示 -log10 转化后的 padj
p <- ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
geom_point(size = 1) +  #绘制散点图
scale_color_manual(values = c('red', 'gray', 'green'), limits = c('up', 'none', 'down')) +  #自定义点的颜色
labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = 'control vs treat', color = '') +  #坐标轴标题
theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), #背景色、网格线、图例等主题修改
    panel.background = element_rect(color = 'black', fill = 'transparent'), 
    legend.key = element_rect(fill = 'transparent')) +
geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +  #添加阈值线
geom_hline(yintercept = 2, lty = 3, color = 'black') +
xlim(-12, 12) + ylim(0, 35)  #定义刻度边界

p