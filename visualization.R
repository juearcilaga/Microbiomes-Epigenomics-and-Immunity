
library(dplyr)
library(readr)
library('igraph')

#INPUT
#Ego edgelist
#Simplified nodes transcript annotation (merged fields)

celltypes=c("Naive_B","Total_CD4_Activated",
          "Total_CD4_NonActivated","Total_CD4_MF","Total_CD8","Naive_CD4","Total_B")

#PARAMETERS
for (i in celltypes){

celltype=i

#celltype=celltypes[1]

nodes_txt=sprintf("%s_mutations_net_nodes.txt",celltype)
nodes <- scan(nodes_txt)

for (node in nodes){

#node=nodes[1]

#Reading ego networks for each node
ego_net_txt=sprintf("ego_%s_%s_edgelist_PID.txt",node, celltype)
ego_net=read.delim(ego_net_txt,sep=" ")
colnames(ego_net)<-c("bait","oe","weight")



#OUTPUT file names

#For exon highlights
net_name<-sprintf("net_%s_exon_%s.graphml",celltype,node)
net_exone_details<-sprintf("net_%s_exon_details_%s.graphml",celltype,node)
plot_exon_name<-sprintf("%s_exon_%s.png",celltype,node)
main_exon<-sprintf("Exons Highlighted in %s:%s",celltype,node)
plot_exone_details<-sprintf("%s_exon_details_%s.png",celltype,node)
net_exon_atributes_txt<-sprintf("%s_exon_atributes_%s.tsv",celltype,node)



#For gene highlights
net_gene_name<-sprintf("net_%s_gene_%s.graphml",celltype,node)
net_gene_details<-sprintf("net_%s_gene_details_%s.graphml",celltype,node)
plot_gene_name<-sprintf("%s_gene_%s.png",celltype,node)
plot_gene_details<-sprintf("%s_gene_details_%s.png",celltype,node)
main_genes<-sprintf("Genes Highlighted in %s:%s",celltype,node)
net_gene_atributes_txt<-sprintf("%s_gene_atributes_%s.tsv",celltype,node)



#Reading annotations for the ego networks

#For exon highlights
exon_annot_txt=sprintf("%s_%s_sum_net_transcript.txt",node, celltype)
exon_annot <- read_table2(exon_annot_txt)
with_exon<-data.frame(exon_annot)

#For gene highlights
gene_annot_txt=sprintf("%s_%s_sum_net_gene.txt",node, celltype)
gene_annot<- read.delim(gene_annot_txt, sep="\t")
with_gene<-gene_annot





#fill with NAs for nodes without exons
baits<-data.frame(ID=ego_net$bait)
oes<-data.frame(ID=ego_net$oe)
IDS<-rbind(baits,oes)
IDS_unique<-unique(IDS)


'%notin%' <- Negate('%in%')

#Select those nodes in the network that do not contain transcripts
ID_filtered<-filter(IDS_unique, IDS_unique$ID %notin% with_exon$ID) 

nones=data.frame()
for (i in ID_filtered){
  j<-seq_along(i)
  nones[j,1]<-data.frame(ID=i)
  nones[j,2]<-data.frame(transcript_IDs="None") 
  nones[j,3]<-data.frame(chr="None") 
  nones[j,4]<-data.frame(start="None")
  nones[j,5]<-data.frame(end="None")
  }

mixed<-rbind(with_exon,nones)



#Build network
net <-graph_from_data_frame(d=ego_net, vertices=mixed, directed=F)

c <-layout_in_circle(net)

V(net)[transcript_IDs=="None"]$color <- "gray50"
V(net)[transcript_IDs!="None"]$color <- "tomato"
V(net)[name==node]$color <-"gold"

V(net)$size <- ifelse(V(net)$name <100, 6, 4)

#Export network
write_graph(net, net_name, format =  "graphml")

par(mar=c(4,0,1,0))

#Plotting network
plot(net, main= main_exon, edge.arrow.size=.4,vertex.label=NA, layout=c)

#It helps to add a legend explaining the meaning of the colors we used:
legend(x=0.1,y=-1, c("without_exons","with_exons", "with_mutation"),xjust=0.5, pch=21,col="#777777", pt.bg=c("gray50", "tomato", "gold"), pt.cex=1, cex=.8, bty="n", ncol=1)



#l <-layout_on_sphere(net)
#f <-layout_with_fr(net)
#plot(net, edge.arrow.size=.4,vertex.label=NA, layout=f)
#plot(net, edge.arrow.size=.4,vertex.label=NA, layout=l)



#Export plot
png(plot_exon_name, width = 12, height = 12, units = "cm", res = 600, pointsize = 10)
V(net)$size <-ifelse(V(net)$name <100, 6, 4)

plot(net, main= main_exon, edge.arrow.size=.4,vertex.label=NA, layout=c)

legend(x=0, y=-1, xjust = 0.5, c("with_exons","without_exons", "with_mutation"), pch=21,col="#777777", pt.bg=c("gray50", "tomato", "gold"), pt.cex=1, cex=.8, bty="n", ncol=1)

dev.off()







#Graph with exon IDS


library(ggplot2);library(scales)

radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2* pi), range(x)))
}

n <- length(V(net)$name)

lab.locs <- radian.rescale(x=1:n, direction=-1,  start=0)

library(Polychrome)

# create your own color palette based on `seedcolors`

P100= createPalette(100,  c("#ff0000", "#00ff00", "#0000ff"))

colores<-P100


V(net)[transcript_IDs=="None"]$color <- "gray50"
V(net)[transcript_IDs!="None"]$color <-colores[V(net)[transcript_IDs!="None"]]
V(net)[name==node]$color <-"gold"

V(net)$size <- ifelse(V(net)$name <100, 6, 4)

V(net)$idex <-seq(1,length(V(net)$name))
  
vlabel_tr<-ifelse(V(net)$name %notin% V(net)[transcript_IDs=="None"]$name, V(net)$idex, NA)

par(mar=c(4,0,1,0))

if(V(net)$name <100) {
  vdist=1
} else if(V(net)$name<150) {
  vdist=1.5
} else {
  vdist=2.5
}


vsize=ifelse(V(net)$name <100, .7, .6)

vlegendatr<-ifelse(V(net)$color %in% V(net)$color[1:length(V(net)[transcript_IDs!="None"]$name)], V(net)$idex, NA)

plot(net,main= main_exon, edge.arrow.size=.4, vertex.label = vlabel_tr, vertex.label.color="black", vertex.label.cex=0.7, vertex.label.dist=1.5, vertex.label.degree=lab.locs, layout=c)

legend(x=0, y=-1, xjust = 0.5, vlegendatr[1:length(V(net)[transcript_IDs!="None"]$name)], pch=21,col="#777777", pt.bg=V(net)[transcript_IDs!="None"]$color, pt.cex=1, cex=.7, bty="n", ncol=13)


#Export network
write_graph(net,net_exone_details, format =  "graphml")



net_exon_atributes<-data.frame(vertex_attr(net))
write_tsv(net_exon_atributes, net_exon_atributes_txt)




png(plot_exone_details, width = 12, height = 12, units = "cm", res = 600, pointsize = 10)

plot(net, main= main_exon,edge.arrow.size=.4, vertex.label = vlabel_tr, vertex.label.color="black", vertex.label.cex=vsize, vertex.label.dist=vdist, vertex.label.degree=lab.locs, layout=c)

legend(x=0, y=-1, xjust = 0.5, vlegendatr[1:length(V(net)[transcript_IDs!="None"]$name)], pch=21,col="#777777", pt.bg=V(net)[transcript_IDs!="None"]$color, pt.cex=1, cex=.7, bty="n", ncol=13)

dev.off()













#Analysis for genes


#Select those nodes in the net_geneswork that do not contain genes
ID_filtered_genes<-filter(IDS_unique, IDS_unique$ID %notin% with_gene$ID) 

nones_genes=data.frame()
for (i in ID_filtered_genes){
  j<-seq_along(i)
  nones_genes[j,1]<-data.frame(ID=i)
  nones_genes[j,2]<-data.frame(gene_IDs="None") 
  nones_genes[j,3]<-data.frame(chr="None") 
  nones_genes[j,4]<-data.frame(start="None")
  nones_genes[j,5]<-data.frame(end="None")
  }

mixed_genes<-rbind(with_gene,nones_genes)



net_genes <-graph_from_data_frame(d=ego_net, vertices=mixed_genes, directed=F)

c <-layout_in_circle(net_genes)

V(net_genes)[gene_IDs=="None"]$color <- "gray50"
V(net_genes)[gene_IDs!="None"]$color <- "tomato"
V(net_genes)[name==node]$color <-"gold"
V(net_genes)$size <- ifelse(V(net_genes)$name <100, 6, 4)

par(mar=c(4,0,1,0))

plot(net_genes, main= main_genes, edge.arrow.size=.4,vertex.label=NA, layout=c)

legend(x=0, y=-1, xjust = 0.5, c("without_genes","with_genes", "with_mutation"), pch=21,col="#777777", pt.bg=c("gray50", "tomato", "gold"), pt.cex=1, cex=.8, bty="n", ncol=1)

#Export network
write_graph(net_genes,net_gene_name, format =  "graphml")



png(plot_gene_name, width = 12, height = 12, units = "cm", res = 600, pointsize = 10)

plot(net_genes,main= main_genes, edge.arrow.size=.4,vertex.label=NA, layout=c)

legend(x=0, y=-1, xjust = 0.5, c("without_genes","with_genes", "with_mutation"), pch=21,col="#777777", pt.bg=c("gray50", "tomato", "gold"), pt.cex=1, cex=.8, bty="n", ncol=1)

dev.off()




colores[V(net_genes)[gene_IDs!="None"]]

V(net_genes)[gene_IDs=="None"]$color <- "gray50"
V(net_genes)[gene_IDs!="None"]$color <-colores[V(net_genes)[gene_IDs!="None"]]
V(net_genes)[name==node]$color <-"gold"

V(net_genes)$size <- ifelse(V(net_genes)$name <100, 6, 4)

V(net_genes)$idex <-seq(1,length(V(net_genes)$name))
  
vlabel<-ifelse(V(net_genes)$name %notin% V(net_genes)[gene_IDs=="None"]$name, V(net_genes)$idex, NA)

par(mar=c(4,0,1,0))

vsizeg=ifelse(V(net_genes)$name <100, .7, .6)

if(V(net_genes)$name <100) {
  vdistg=1
} else if(V(net_genes)$name<150) {
  vdistg=1.5
} else {
  vdistg=2.5
}

vlegenda<-ifelse(V(net_genes)$color %in% V(net_genes)$color[1:length(V(net_genes)[gene_IDs!="None"]$name)], V(net_genes)$idex, NA)

plot(net_genes, main=main_genes, edge.arrow.size=.4, vertex.label = vlabel, vertex.label.color="black", vertex.label.cex=0.7, vertex.label.dist=1.5, vertex.label.degree=lab.locs, layout=c)

legend(x=0, y=-1, xjust = 0.5, vlegenda[1:length(V(net_genes)[gene_IDs!="None"]$name)], pch=21,col="#777777", pt.bg=V(net_genes)[gene_IDs!="None"]$color, pt.cex=1, cex=.7, bty="n", ncol=13)


#Export network
write_graph(net_genes,net_gene_details, format =  "graphml")



net_gene_atributes<-data.frame(vertex_attr(net_genes))

write_tsv(net_gene_atributes, net_gene_atributes_txt)







png(plot_gene_details, width = 12, height = 12, units = "cm", res = 600, pointsize = 10)

plot(net_genes, main=main_genes, edge.arrow.size=.4, vertex.label = vlabel, vertex.label.color="black", vertex.label.cex=vsizeg, vertex.label.dist=vdistg, vertex.label.degree=lab.locs, layout=c)

legend(x=0, y=-1, xjust = 0.5, vlegenda[1:length(V(net_genes)[gene_IDs!="None"]$name)], pch=21,col="#777777", pt.bg=V(net_genes)[gene_IDs!="None"]$color, pt.cex=1, cex=.7, bty="n", ncol=13)

dev.off()

}
}


