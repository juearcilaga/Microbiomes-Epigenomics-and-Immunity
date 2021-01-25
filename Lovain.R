Mono_score<-read.delim("Monocytes_CS_eQTL_FANTOM_hg37Genes")
str(Mono_score)
#'data.frame':	96442 obs. of  20 variables:
Mono_score<-cbind(Mono_score[["ID"]],Mono_score[["Propagation"]])
colnames(Mono_score)<-c("ID","Propagation")
head(Mono_score)

Edit4 <-read.delim("EDIT4_2.tsv")
colnames(Edit4)
# [1] "baitChr"                "baitStart"              "baitEnd"               
# [4] "baitID"                 "baitName"               "oeChr"                 
#[7] "oeStart"                "oeEnd"                  "oeID"                  
#[10] "oeName"                 "dist"                   "Monocytes"             
#[13] "Macrophages_M0"         "Macrophages_M1"         "Macrophages_M2"        
#[16] "Neutrophils"            "Megakaryocytes"         "Endothelial_precursors"
#[19] "Erythroblasts"          "Foetal_thymus"          "Naive_CD4"             
#[22] "Total_CD4_MF"           "Total_CD4_Activated"    "Total_CD4_NonActivated"
#[25] "Naive_CD8"              "Total_CD8"              "Naive_B"               
#[28] "Total_B"

str(Edit4)
#'data.frame':	728835 obs. of  28 variables:
Edit4[,4], Edit4[,9]

str(Edit4[,4])
#int [1:728835] 218 218 218 218 218 218 219 219 220 220 ...
str(Edit4[["baitID"]])
# int [1:728835] 218 218 218 218 218 218 219 219 220 220 ...
 
Mono_net_raw <-cbind(Edit4[["baitID"]],Edit4[["oeID"]],Edit4[["Monocytes"]] )
colnames(Mono_net_raw)<-c("baitID","oeID","Monocytes")

#Tcell_net <-cbind(Edit4[["baitID"]],Edit4[["oeID"]],Edit4[["XXX"]] )
#colnames(Tcell_net)<-c("baitID","oeID","XXXX")

library(tidyr) 
Mono_net_raw<-data.frame(Mono_net_raw)
str(Mono_net_raw)
#'data.frame':	728835 obs. of  3 variables:
Mono_net_filtered<-filter(Mono_net_raw, Monocytes >= 5 )#Is >= or just >?
str(Mono_net_filtered)
#'data.frame':	171430 obs. of  3 variables:


#For Luvain 

colnames(Mono_net_filtered)<-c("ID","oeID","Monocytes")
Mono_net_propag<-left_join(Mono_net_filtered,Mono_score,by="ID",copy=TRUE)

head(Mono_net_propag)
#   ID oeID Monocytes Propagation
#1 218  220 25.546568    1.467695
#2 218  224  5.050248    1.467695
colnames(Mono_net_propag)<-c("baitID","oeID","Monocytes_Chicago","Propagation_bait")



colnames(Mono_net_propag)<-c("baitID","ID","Monocytes_Chicago","Propagation_bait")
Mono_net_propag<-left_join(Mono_net_propag,Mono_score,by="ID",copy=TRUE)
colnames(Mono_net_propag)<-c("baitID","oeID","Monocytes_Chicago","Propagation_bait", "Propagation_oe")

Mono_net_propag$weight <- rowSums(Mono_net_propag[,4:5] )

Mono_net<-cbind(Mono_net_propag[["baitID"]],Mono_net_propag[["oeID"]],Mono_net_propag[["weight"]] )
colnames(Mono_net)<-c("baitID","oeID","weight")
Mono_net<-data.frame(Mono_net)
str(Mono_net)
#'data.frame':	171430 obs. of  3 variables:

#For Neutrophils
Neut_net <-cbind(Edit4[["baitID"]],Edit4[["oeID"]],Edit4[["Neutrophiles"]] )
colnames(Neut_net)<-c("baitID","oeID","Neutrophiles")
head(Neut_net)

Neut_score<-read.delim("Neutrophils_CS_eQTL_FANTOM_hg37Genes")
str(Neut_score)
#'data.frame':	96442 obs. of  20 variables:
Neut_score<-cbind(Neut_score[["ID"]],Neut_score[["Propagation"]])
colnames(Neut_score)<-c("ID","Propagation")
head(Neut_score)






