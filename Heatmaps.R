
#Eliminar variables que hayan quedado abiertas de una sesion anterior
rm(list=ls());



#Cargar las librerias necesarias
library(extrafont);#Para que todas las imagenes tengan times new roman
par(family = "Times New Roman")
#library(car);library(ggpubr);library(viridis);library(vegan);library(pander);library("ggfortify")
library(ggplot2);library(readxl);
library(grid);
library(pheatmap)
library(viridis);
library(gdata)
#library(gtools)




library("plyr")
library(tidyr);library(dplyr)



#List the cell types to analyze, the names should be the same that appear in the columns of EDIT4.txt file, wich contains the PC-HIC data

celltypes=c("Naive_CD4","Total_CD4_MF", "Total_CD4_Activated", "Total_CD4_NonActivated", "Naive_CD8", "Total_CD8","Naive_B","Total_B","Foetal_thymus")




{r message=FALSE}
my_data=list()
#colnames(my_data)<-c("gene","name","count","node","celltype")

counter=0
for (i in celltypes){
celltype=i
nodes_txt=sprintf("%s_mutations_net_nodes.tsv",celltype)
nodes <- scan(nodes_txt) #("all-present?")
for (node in nodes){
counter= counter+1
genes_txt=sprintf("%s_%s_gene_list.txt",node, celltype)
genes_read<- read.delim(genes_txt, sep ="\t")
genes<-genes_read$gene_ID
genes_name<-genes_read$gene_name
counts<-data.frame(gene=genes, name=genes_name, count=1)
counts$node<-node
counts$cell<-celltype
my_data[[counter]]<-counts
}
}

Alldata<-do.call("rbind",my_data)

all_geneID<-unique(Alldata$gene)#137
all_gene_name<-unique(Alldata$name)#135 (hay dos ID con el mismo nombre)#no se como identificarlos, computacionalmente
all_nodes<-unique(Alldata$node)#7 nodes

#name<-Alldata$name
#sort(name)



heatmap<-matrix(0,length(all_gene_name), length(celltypes))
colnames(heatmap)<-celltypes
rownames(heatmap)<-all_gene_name



#In gene presence all tissues
for (each in celltypes){
test<-filter(Alldata, Alldata$cell==each)
genes<-unique(test$name)
  for (gene in all_gene_name){
    if(gene %in% genes){
      heatmap[gene,each]=1
      }
      }
}


{r include=FALSE}
gene_in_tissue<-heatmap
png("heatmap_tissue_presence.png", width = 15, height = 45, units = "cm", res = 600, pointsize = 10)
pheatmap(gene_in_tissue, color = c("beige","orange"))
dev.off()

#Presence of the gene in at least one of the network mutations by celltype, 



#In gene presence all tissues, numb of mutations
for (each in celltypes){
test<-filter(Alldata, Alldata$cell==each)
genes<-unique(test$name)
frequency<-data.frame(count(test, name))
  for (gene in all_gene_name){
    if(gene %in% genes){
      heatmap[gene,each]=frequency[frequency$name==gene ,2]
      }
      }
}




{r echo=TRUE}
gene_in_tissue_quant<-heatmap

png("heatmap_tissue_quant.png", width = 15, height = 45, units = "cm", res = 600, pointsize = 10)
pheatmap(gene_in_tissue_quant, color=c("beige", "grey", "yellow","orange", "red","blue"))
dev.off()




heatmap2<-matrix(0,length(all_gene_name), length(all_nodes))
rownames(heatmap2)<-all_gene_name
colnames(heatmap2)<-all_nodes



#In gene presence all nodes
for (each in all_nodes){
test<-filter(Alldata, Alldata$node==each)
genes<-unique(test$name)
for (gene in all_gene_name){
  if(gene %in% genes){
    node<-sprintf("%d",each)
    heatmap2[gene,node]=1
    }
    }
}

mutations_presence<-heatmap2
colnames(mutations_presence)<-c("chr14:20343467_A_C","chr14:20180243_T_A","chr14:19815670_G_GT", "chr14:20725578_C_T","chr18:51130565_TACATATATACATACATATATAC_*","chr14:20772874_A_G", "chr14:20762812_G_A")

png("heatmap2_gene_in_mutation.png", width = 15, height = 50, units = "cm", res = 600, pointsize = 10)
pheatmap(mutations_presence, color = c("beige","orange"))
dev.off()





#In gene presence all mutations, numb of celltypes
for (each in all_nodes){
test<-filter(Alldata, Alldata$node==each)
genes<-unique(test$name)
frequency<-length(unique(test$cell))
print(frequency)
for (gene in all_gene_name){
 if(gene %in% genes){
  node<-sprintf("%d",each)
  heatmap2[gene,node]=frequency
    }
}
}




mutations_freq_tissue<-heatmap2
colnames(mutations_freq_tissue)<-c("chr14:20343467_A_C","chr14:20180243_T_A","chr14:19815670_G_GT", "chr14:20725578_C_T","chr18:51130565_TACATATATACATACATATATAC_*","chr14:20772874_A_G", "chr14:20762812_G_A")


png("heatmap2_mutations_freq_tissue.png", width = 15, height = 50, units = "cm", res = 600, pointsize = 10)
pheatmap(mutations_freq_tissue, color=c("beige", "grey", "yellow","orange", "red","blue","brown","green","pink"))
dev.off()


