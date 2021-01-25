#sudo apt install -y python3-pip
#source my_env/bin/activate
#Note: Within the virtual environment, you can use the command python instead of python3, and pip instead of pip3 if you would prefer. If you use Python 3 on your machine outside of an environment, you will need to use the python3 and pip3 commands exclusively.

#source Manninder/bin/activate
#pip install networkx

import networkx as nx
import pandas as pd

Mononet = pd.read_csv('Mono_net.csv')
print (Mononet)

Graphtype=nx.Graph()



#FUNCTIONS
def from_pandas_dataframe(df, source, target, edge_attr=None):
        
	import networkx as nx
        
	g = nx.Graph()

    # Index of source and target
	src_i = df.columns.get_loc(source)
	tar_i = df.columns.get_loc(target)
	
	if edge_attr:
        # If all additional columns requested, build up a list of tuples
        # [(name, index),...]
		if edge_attr is True:
            # Create a list of all columns indices, ignore nodes
			edge_i = []
			for i, col in enumerate(df.columns):
				if col is not source and col is not target:	
					edge_i.append((col, i))
        				# If a list or tuple of name is requested
		elif isinstance(edge_attr, (list, tuple)):
            		edge_i = [(i, df.columns.get_loc(i)) for i in edge_attr]
        		# If a string or int is passed
		else:
            		edge_i = [(edge_attr, df.columns.get_loc(edge_attr)),]

        		# Iteration on values returns the rows as Numpy arrays
		for row in df.values:
            		g.add_edge(row[src_i], row[tar_i], {i:row[j] for i, j in edge_i})
    
    			# If no column names are given, then just return the edges.
	else:
		for row in df.values:
			g.add_edge(row[src_i], row[tar_i])

	return g
	
	
import numpy as np
G=from_pandas_dataframe(Mononet,'baitID', 'oeID')
#G['E']['C']['weight']
#    10
#    >>> G['E']['C']['cost']
#    9

import dill
filename = 'globalsave.pkl'
dill.dump_session(filename)


import community as community_louvain
import matplotlib.cm as cm
import matplotlib.pyplot as plt


# compute the best partition
partition = community_louvain.best_partition(G)

# draw the graph
pos = nx.spring_layout(G)
# color the nodes according to their partition
cmap = cm.get_cmap('viridis', max(partition.values()) + 1)
nx.draw_networkx_nodes(G, pos, partition.keys(), node_size=40,
                       cmap=cmap, node_color=list(partition.values()))
nx.draw_networkx_edges(G, pos, alpha=0.5)
plt.show()


import dill
filename = 'globalsave.pkl'
dill.dump_session(filename)

# and to load the session again:
#dill.load_session(filename)


