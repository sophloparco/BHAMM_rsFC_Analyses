import pandas
import numpy
import math
import nilearn
import nistats
import os
import sys
import matplotlib.pyplot as plt
import matplotlib
import glob2 as glob
import re
from collections import OrderedDict
from datetime import date, datetime
import networkx
import SchaeferAtlas
import ConnectivityMatrix


def make_graph(matrix, atlas, edge_density=.05):

	G= networkx.from_numpy_matrix(matrix.to_numpy(), parallel_edges=False)

	print("\nRaw Graph:")
	print(networkx.info(G))

	for node in G.nodes():
		if G.has_edge(node, node):
			G.remove_edge(node, node)


	Diagonal_ZeroMat = networkx.to_numpy_matrix(G)
	print("\nDiagonal of matrix set to zero:")
	print(Diagonal_ZeroMat)

	#COMPUTE THRESHOLD
	thresh = numpy.percentile(numpy.absolute(Diagonal_ZeroMat), float(100)*(1-edge_density))

	#PRINT IT
	print("\nThreshold for top %s percent of edges: %s" %(str(int(edge_density*100)),str(thresh)))


	#TURN THE NUMPY MATRIX OBJECT INTO PANDAS DATAFRAME SO WE CAN THRESHOLD EASILY
	DZ_mat_df = pandas.DataFrame(Diagonal_ZeroMat)

	#THRESHOLD EDGES BASED OFF PERCENTILE WEIGHT CUTOFF (i.e. ONLY edge_density % OF EDGES HAVE A WEIGHT GREATER THAN OUR THRESH -- ALL OTHERS WILL BE SET TO ZERO)
	thresholded_df = DZ_mat_df.where(DZ_mat_df>thresh, 0)

	return G, thresh

def threshold_graph(g, thresh=0, neg_edges=True):

	from scipy import integrate
	alpha=.05
	backbone_graph = networkx.Graph()
	for node in g:
		k_n = len(g[node])
		if k_n > 1:
			sum_w = sum( g[node][neighbor]['weight'] for neighbor in g[node] )
			for neighbor in g[node]:
				edgeWeight = g[node][neighbor]['weight']
				if neg_edges==True:

					if (float(numpy.absolute(edgeWeight))>thresh)&(not(backbone_graph.has_edge(node,neighbor))):

						backbone_graph.add_edge(node,neighbor, weight = edgeWeight, spring_weight=(1 - edgeWeight))

				if neg_edges==False:
					
					if (float(edgeWeight)>thresh)&(not(backbone_graph.has_edge(node,neighbor))):
	
						backbone_graph.add_edge(node,neighbor, weight = edgeWeight, spring_weight=(1 - edgeWeight))


	return backbone_graph

def degree_threshold(g, k=1):

	from scipy import integrate
	#UNFINISHED
	backbone_graph = networkx.Graph()
	for node in g:
		k_n = len(g[node])
		if k_n > k:
			top_k = [g[node][neighbor]['weight'] for neighbor in g[node]]
			for neighbor in g[node]:
				edgeWeight = g[node][neighbor]['weight']
				if neg_edges==True:

					if (float(numpy.absolute(edgeWeight))>thresh)&(not(backbone_graph.has_edge(node,neighbor))):
						
						# pij = float(edgeWeight)/sum_w
						# if (1-pij)**(k_n-1) < alpha: # equation 2 Serrano et al. 2008 (https://www.pnas.org/content/106/16/6483)

						backbone_graph.add_edge(node,neighbor, weight = edgeWeight, spring_weight=(1 - edgeWeight))

				if neg_edges==False:
					if (float(edgeWeight)>thresh)&(not(backbone_graph.has_edge(node,neighbor))):
						
						# pij = float(edgeWeight)/sum_w
						# if (1-pij)**(k_n-1) < alpha: # equation 2 Serrano et al. 2008 (https://www.pnas.org/content/106/16/6483)

						backbone_graph.add_edge(node,neighbor, weight = edgeWeight, spring_weight=(1 - edgeWeight))


	return backbone_graph

def color_nodes(Graph, atlas):
	attributes={}

	network_num = str(len(atlas['Network'].unique()))

	print("Network res")
	print(network_num)

	syst_color = SchaeferAtlas.color_networks(atlas, "8")

	print("syst_color")
	print(syst_color)

	atlas = SchaeferAtlas.unpack_node_labels(atlas)
	for node in Graph.nodes():
		attributes[node]={'Network':atlas.loc[atlas['Node']==node, 'Network'].iloc[0], 
						'ColorCode':atlas.loc[atlas['Node']==node, 'ColorCode'].iloc[0], 
						'Color':atlas.loc[atlas['Node']==node, 'Color'].iloc[0], 
						'Coords':atlas.loc[atlas['Node']==node, 'Coords'].iloc[0]}

	networkx.set_node_attributes(Graph, attributes)


	return Graph

def get_attribute_list(Graph, attribute='Network'):
	attribute_dict = networkx.get_node_attributes(Graph, attribute)
	attribute_list = [attribute_dict[node] for node in Graph.nodes()]
	return attribute_list

def subgraph(g, attribute='Network', attribute_classes=[]):
	if len(attribute_classes) >0:

		nodelist=[]
		for node in g.nodes():

			if g.nodes[node][attribute] in attribute_classes:
				nodelist.append(node)
		sub_g = networkx.subgraph(g, nodelist)

		return sub_g
	else: 
		return None

def make_cmap(G, Networks, Names=False):
	nodelist = [node for node in G.nodes]
	colormap=[]
	colornames=[]

	for node in nodelist:
		if (G.nodes[node]['Network'] not in Networks):
			G.remove_node(node)
		else:
			colormap.append(G.nodes[node]['ColorCode'])
			colornames.append(G.nodes[node]['Color'])

	z = range(0,len(colormap)) 
	if Names==True:
		return colornames
	else:
		cm = matplotlib.colors.LinearSegmentedColormap.from_list("net_colors", colormap, N=len(colormap))
		return cm 

def plot_graph(Graph, atlas, Networks=[], title="Add a title you dingus!"):

	print("\nPlotting spring-embedded graph for the following networks:")
	print(Networks)
	
	cm = make_cmap(Graph, Networks)


	fig = plt.figure(figsize=(10,8))
	fig.patch.set_facecolor((0.5,0.5,0.5))
	fig.suptitle('%s'%(title), y=.98, size='x-large', fontweight='bold', color='white')
	ax= plt.subplot(111)
	ax.set_facecolor((0.5,0.5,0.5))
	plt.gca().axes.get_yaxis().set_visible(False)
	plt.gca().axes.get_xaxis().set_visible(False)

	
	pos= networkx.spring_layout(Graph, 6*(float(1)/len(Graph.nodes())**(.5)), weight='spring_weight')

	gray= tuple((0.5,0.5,0.5))
	networkx.draw_networkx_edges(Graph, pos, edge_color='white', alpha=.4, ax=ax)

	networkx.draw_networkx_nodes(Graph, pos, cmap=cm,
							node_color=range(len(Graph.nodes())),
							node_size=30, ax=ax)

	labels= {}
	for node in Graph.nodes():
		labels[node]= '%s'%(str(node))							

	return fig, Graph

