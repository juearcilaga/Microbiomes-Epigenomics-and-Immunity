#!/usr/bin/env python

#Load necessary libraries
import networkx as nx
import pandas as pd
import numpy as np


def builtGraph(dataframe):#Debe tener en la primera y segunda fila las interacciones y en la tercera pesos
    Matrix=pd.DataFrame(dataframe).to_numpy()
    rowsnumb = len(Matrix)
    g = nx.Graph()
    
    for i in list(range(rowsnumb)):
        g.add_edge(Matrix[i,0],Matrix[i,1],weight=Matrix[i,2])
    return g

def ego_graphs_gene(graph,nodes_gene,radio):#graph es la grafica grande, node is a list, radio=2
    graphs={}#Create an empty directory
    for j in nodes_list:
        ego="%s"%j
        graphs[j]= nx.ego_graph(graph, int(j), radius=radio, center=True, undirected=True, distance=None)
    return graphs

def ego_graphs_export(graphs,celltype):#argument is a directory
    for key  in graphs:
    	nx.write_weighted_edgelist(graphs[key], "ego_%s_%s_edgelist_PID.txt"%(key,celltype))
    return graphs


#Input:

celltypes=["Naive_CD4","Total_CD4_MF", 
         "Total_CD4_Activated", "Total_CD4_NonActivated", 
         "Naive_CD8", "Total_CD8", "Naive_B", "Total_B", "Foetal_thymus"]



radio=2 #the label of interactions to include, distance to the egocentric node.

celltype= [f for f in celltypes]

for f in celltypes:
    
    #Input:
    file="%s_net.txt"%f
    file3="%s_mutations_net_nodes.txt"%f
    
    #Open each weighted_edge_list as dataframe
    dataframe = pd.read_csv(file, sep='\t') 

    #Build a graph with the dataframe
    graph= builtGraph(dataframe)
    
    #Open each file of nodes list as a python list
    open_file= open(file3, "r")
    read_file= open_file.read()
    nodes_list = read_file.rstrip('\n').split(" ")
    open_file.close()
    
    #Generate an ego graph for each node that contains the mutation of interest
    graphs= ego_graphs_gene(graph,nodes_list,radio)

    #Export the list of nodes in a tab sep file to use later in the vizualization step. 
    ego_graphs_export(graphs,f)

