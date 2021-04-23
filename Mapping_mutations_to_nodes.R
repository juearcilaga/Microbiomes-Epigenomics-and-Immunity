#This script is useful for mapping a set of mutations to the nodes in
#the network of celltypes of interest

library(tidyr);library(dplyr);library(readr);library(GenomicRanges);
library(plyranges);library("stringr");



cells<-c("Naive_CD4","Total_CD4_MF", 
         "Total_CD4_Activated", "Total_CD4_NonActivated", 
         "Naive_CD8", "Total_CD8", "Naive_B", "Total_B", "Foetal_thymus")



for (i in cells){
  
#PARAMETERS  
celltype=i
mutations_net_nodes_txt=sprintf("%s_mutations_net_nodes.txt",i)
mutations_to_nodes_txt=sprintf("%s_mutations_to_nodes_coord.txt",i)
Celltype_net_txt=sprintf("%s_net.txt", i)
#INPUTS
Edit4 <-read.delim("EDIT4_2.tsv")


Celltype_net_raw <-cbind(Edit4[["baitID"]],Edit4[["oeID"]],Edit4[[celltype]] )

Celltype_net_raw<-data.frame(Celltype_net_raw)

colnames(Celltype_net_raw)<-c("baitID","oeID", "weight")

baits<-data.frame(ID=Edit4$baitID, chr=Edit4$baitChr, start=Edit4$baitStart, 
                  end=Edit4$baitEnd)

oes<-data.frame(ID=Edit4$oeID, chr=Edit4$oeChr, start=Edit4$oeStart, 
                end=Edit4$oeEnd)

IDS<-rbind(baits,oes)

IDS_unique<-unique(IDS)



Celltype_net_filtered<-filter(Celltype_net_raw, Celltype_net_raw$weight >= 5 )
#Is >= or just >?
colnames(Celltype_net_filtered)<-c("baitID","oeID","weight")

ID_filtered<-filter(IDS_unique, 
                    (IDS_unique$ID %in% Celltype_net_filtered$oeID) | 
                        (IDS_unique$ID %in%  Celltype_net_filtered$oeID))

write_tsv(Celltype_net_filtered, Celltype_net_txt)

gr_celltype <- GRanges(
    seqnames = ID_filtered$chr,
    ranges = IRanges(ID_filtered$start, end = ID_filtered$end, 
                     names =ID_filtered$ID),
    fragment_ID=ID_filtered$ID,
    fragment_start=ID_filtered$start,
    fragment_end=ID_filtered$end,
    fragment_size=ID_filtered$end-ID_filtered$start)

gr_celltype



mutations <-read.delim("mutations.txt")
mutations_hg37 <-read.delim("mutations_hg37.txt")

#La mala noticia esque esto esta anotado con el g37 no con el g38
colnames(mutations)<-c("chr","pos", "ID","REF","ALT")

#La mala noticia esque esto esta anotado con el g37 no con el g38
#Ahora si cambiadas las coordenadas con liftover g37
mutations_hg37$ID<-mutations$ID
mutations_hg37$REF<-mutations$REF
mutations_hg37$ALT<-mutations$ALT




gr_mutations<- GRanges(
    seqnames = mutations_hg37$chr,
    ranges = IRanges(mutations_hg37$start, end = mutations_hg37$end, 
                     names =mutations_hg37$ID),
    REF=mutations_hg37$REF,
    ALT=mutations_hg37$ALT,
    mutation_names =mutations_hg37$ID)
gr_mutations



mutations_to_frags<-join_overlap_inner(gr_celltype, gr_mutations)
mapped2<- data.frame(seqnames(mutations_to_frags),ranges(mutations_to_frags),
                     values(mutations_to_frags))

nodes2<-mapped2$fragment_ID
nodes2<-unique(nodes2)



write(nodes2, mutations_net_nodes_txt)
write_tsv(mapped2, mutations_to_nodes_txt)
}


