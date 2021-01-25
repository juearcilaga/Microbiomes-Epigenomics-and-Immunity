#sudo apt install -y python3-pip
#source my_env/bin/activate
#Note: Within the virtual environment, you can use the command python instead of python3, and pip instead of pip3 if you would prefer. If you use Python 3 on your machine outside of an environment, you will need to use the python3 and pip3 commands exclusively.

#source Manninder/bin/activate
#pip install networkx

 
import networkx as nx
Mononet = pd.read_csv('Mono_net.csv')       
g = nx.Graph()

Matrix=pd.DataFrame(Mononet).to_numpy()
rowsnumb = len(Matrix)
columnsnumb=len(Matrix[0])


for i in list(range(rowsnumb)):
	print(Matrix[i,0],Matrix[i,1],Matrix[i,2])
	g.add_edge(Matrix[i,0],Matrix[i,1],weight=Matrix[i,2])

import community as community_louvain
import matplotlib.cm as cm
import matplotlib.pyplot as plt


# compute the best partition
partition = community_louvain.best_partition(g)

import dill
filename = 'globalsave.pkl'
dill.dump_session(filename)

# draw the graph
pos = nx.spring_layout(g)

# color the nodes according to their partition
cmap = cm.get_cmap('viridis', max(partition.values()) + 1)

nx.draw_networkx_nodes(g, pos, partition.keys(), node_size=40,
                       cmap=cmap, node_color=list(partition.values()))
nx.draw_networkx_edges(g, pos, alpha=0.5)
plt.show()




