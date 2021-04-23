#Load the necessary libraries

library(tidyr);library(dplyr);library(readr);library(GenomicRanges);
library(plyranges);library("stringr");


#List the cell types to analyze, the names should be the same that appear in the columns of EDIT4.txt file, wich contains the PC-HIC data

celltypes=c("Naive_CD4","Total_CD4_MF", 
                     "Total_CD4_Activated", "Total_CD4_NonActivated", 
                     "Naive_CD8", "Total_CD8", "Naive_B", "Total_B", "Foetal_thymus")
            
            
file = "appris_data.principal_Gencode19Ensembl74.txt"
file2 = "grch37_ensembl_biomart_fixed6.txt"
Edit4<-read.delim("EDIT4_2.tsv", sep="\t")

#Declaring variables


#PARAMETERS
for (i in celltypes){

celltype=i

nodes_txt=sprintf("%s_mutations_net_nodes.txt",celltype)
nodes <- scan(nodes_txt)

for (node in nodes){


#INPUTS
Celltype_net_txt=sprintf("ego_%s_%s_edgelist_PID.txt",node,celltype)
cell_net_gene_txt=sprintf("%s_%s_net_gene.txt",node,celltype)


#OUTPUT files names
cell_net_gene_sum_txt=sprintf("%s_%s_sum_net_gene.txt",node,celltype)
cell_net_transcript_txt = sprintf("%s_%s_net_transcript.txt",node,celltype)
cell_net_transcript_sum_txt = sprintf("%s_%s_sum_net_transcript.txt",node,celltype)
gene_list_txt= sprintf("%s_%s_gene_list.txt",node,celltype)
exones_list_txt= sprintf("%s_%s_exon_list.txt",node,celltype)



#Open the network in EDIT4 an filtering bad quality interactions

Celltype_net_raw <-cbind(Edit4[["baitID"]],Edit4[["oeID"]],Edit4[[celltype]] )
Celltype_net_raw<-data.frame(Celltype_net_raw)
colnames(Celltype_net_raw)<-c("baitID","oeID", "weight")

baits<-data.frame(ID=Edit4$baitID, chr=Edit4$baitChr, start=Edit4$baitStart, end=Edit4$baitEnd)
oes<-data.frame(ID=Edit4$oeID, chr=Edit4$oeChr, start=Edit4$oeStart, end=Edit4$oeEnd)
IDS<-rbind(baits,oes)
IDS_unique<-unique(IDS)

Celltype_net_filtered<-filter(Celltype_net_raw, Celltype_net_raw$weight >= 5 )#Is >= or just >?
colnames(Celltype_net_filtered)<-c("baitID","oeID","weight")

ID_filtered<-filter(IDS_unique, (IDS_unique$ID %in% Celltype_net_filtered$oeID) | (IDS_unique$ID %in%  Celltype_net_filtered$oeID))



#Building ego network for the node and celltype of interest

ego_net<-read.delim(Celltype_net_txt, sep=" ")
colnames(ego_net)<-c("baitID","oeID", "weight")

ego_baits<-data.frame(ID=ego_net$baitID)
ego_oes<-data.frame(ID=ego_net$oeID)
ego_IDS<-rbind(ego_baits,ego_oes)
ego_IDS_unique<-unique(ego_IDS)

ID_filtered_ego<-filter(ID_filtered, ID_filtered$ID %in% ego_IDS_unique$ID)


#Building granges object for the ego network

gr_celltype <- GRanges(
    seqnames = ID_filtered_ego$chr,
    ranges = IRanges(ID_filtered_ego$start, end = ID_filtered_ego$end, names =ID_filtered_ego$ID),
    fragment_ID=ID_filtered_ego$ID,
    fragment_start=ID_filtered_ego$start,
    fragment_end=ID_filtered_ego$end,
    fragment_size=ID_filtered_ego$end-ID_filtered_ego$start)



#Reading annotation files for gene and building granges object with them

Annot_gene_structure<-read.delim(file2)
colnames(Annot_gene_structure)

genes<-data.frame(chr=Annot_gene_structure$Chromosome_scaffold_name, 
                  Gene_start =Annot_gene_structure$Gene_start,
                  Gene_end=Annot_gene_structure$Gene_end, 
                  Gene_name=Annot_gene_structure$Gene_name,
                  Gene_stable_ID_version=Annot_gene_structure$Gene_stable_ID_version)

genes<-unique(genes)

genes$Chr <- paste0("chr", genes$chr)

gr_gene<- GRanges(
    seqnames =genes$Chr,
    ranges = IRanges(genes$Gene_start, end = genes$Gene_end,
                     names = genes$Gene_stable_ID_version),
    gene_name=genes$Gene_name,
    gene_ID=genes$Gene_stable_ID_version,
    gene_start= genes$Gene_start,
    gene_end=genes$Gene_end,
    gene_size=genes$Gene_end-genes$Gene_start)



#Annotation of principal isoforms

# Annotation of principal isoforms
Principal_isoforms<-read.delim(file)

colnames(Principal_isoforms)=c("gene_name", "gene_ID","transcript_ID", 
                               "column_4", "isoform_type")

Principal_isoform_genes= filter(Principal_isoforms, 
                                str_detect(isoform_type, "PRINCIPAL"))

Annot_gene_structure_principal=filter(Annot_gene_structure,
              Transcript_stable_ID %in% Principal_isoform_genes$transcript_ID)

Annot_gene_structure_principal$Chr <- paste0("chr",
                       Annot_gene_structure_principal$Chromosome_scaffold_name)



#building granges object with transcripts principal issoforms

gr_Annot_gene_structure<- GRanges(
    seqnames =Annot_gene_structure_principal$Chr,
    ranges = IRanges(Annot_gene_structure_principal$Exon_region_start_bp, 
                     end = Annot_gene_structure_principal$Exon_region_end_bp, 
                     names = Annot_gene_structure_principal$Exon_stable_ID),
    exon_rank=Annot_gene_structure_principal$Exon_rank_in_transcript,
    transcript_ID=Annot_gene_structure_principal$Transcript_stable_ID_version,
    gene_ID=Annot_gene_structure_principal$Gene_stable_ID_version,
    gene_name=Annot_gene_structure_principal$Gene_name,
    exon_ID=Annot_gene_structure_principal$Exon_stable_ID,
    exon_start= Annot_gene_structure_principal$Exon_region_start_bp,
    exon_end=Annot_gene_structure_principal$Exon_region_end_bp,
    exon_wide=Annot_gene_structure_principal$Exon_region_end_bp-Annot_gene_structure_principal$Exon_region_start_bp)




#Annotating fragments in network with coordinates of gene and transcripts

annotated_frags_genes<-join_overlap_inner(gr_celltype,gr_gene)

annotated<- data.frame(seqnames(annotated_frags_genes),ranges(annotated_frags_genes),values(annotated_frags_genes))

summarized<-aggregate(gene_ID ~ fragment_ID, annotated, paste, collapse = ", ")
colnames(summarized)<-c("ID","gene_IDs")
annotated_sum<-left_join(summarized,ID_filtered_ego,by="ID")

write_tsv(annotated, cell_net_gene_txt)
write_tsv(annotated_sum, cell_net_gene_sum_txt)

list_genes<-data.frame(unique(cbind(gene_name=annotated[["gene_name"]],gene_ID=annotated[["gene_ID"]])))
write_tsv(list_genes,gene_list_txt)


annotated_frags_transcript<-join_overlap_inner(gr_celltype,gr_Annot_gene_structure)

annotated_transcript<- data.frame(seqnames(annotated_frags_transcript),ranges(annotated_frags_transcript),values(annotated_frags_transcript))

summarized2<-aggregate(transcript_ID ~ fragment_ID, annotated_transcript, paste, collapse = ", ")
colnames(summarized2)<-c("ID","transcript_IDs")

annotated_transcript_sum<-left_join(summarized2,ID_filtered_ego,by="ID")

write_tsv(annotated_transcript_sum, cell_net_transcript_sum_txt)
write_tsv(annotated_transcript, cell_net_transcript_txt)#Replace by cell name

list_exones<-unique(annotated_transcript[,10:14])
write_tsv(list_exones,exones_list_txt)

}
}


