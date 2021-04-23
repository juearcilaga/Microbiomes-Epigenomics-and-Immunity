# Microbiomes-Epigenomics-and-Immunity
This repository stores the raw scripts made during the progress of the first year of my PhD. 


##  Side Project1: Interpreting rare variants on non-coding regions.
1.	From the 17 Promoter Capture HiC (PCHI) datasets available published by Javierre et al. 2016, I selected nine datasets adaptative immune cell types (T cells and B cells and thymus)
2. I built nine networks with Javierra et al. data, one network for each o the cell types and mapped the seven mutation positions to the node coordinates in the networks using the script "Mapping_mutation_to_nodes.R"
3. I built egocentric networks (for each mutated node in each cell type  using the the script "Ego_mutations.py"
4. I annotated genes and exons coordinates in the nodes of each ego networks using data downloaded from Ensemble Biomart (appris_data.principal_Gencode19Ensembl74.txt) and filtering principal transcript isoforms are according to the information in the APRIS database (appris_data.principal_Gencode19Ensembl74.txt) using the script "Annotations_ego.R"
5. I used the scripts "Visualization.R" and "Heatmap.R" to generate the images of the networks in circular layout and heaptmaps.
