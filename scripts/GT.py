import pandas
import numpy
import math
import nilearn
import nistats
import os
import sys
import matplotlib.pyplot as plt
import matplotlib
import glob
import re
from collections import OrderedDict
from datetime import date, datetime
import networkx
import scipy
import Graph
from holoviews import opts, dim
import holoviews as hv

def SystemSegregation(graph, atlas):

	G = Graph.threshold_graph(graph, thresh=0, neg_edges=False)
	g = Graph.color_nodes(G, atlas)

	Systems = atlas['Network'].unique()
	
	System_dict= {}
	SS_dict={}
	for system in Systems:
		System_dict[system]={'within':0, 'between':0}
	
		within_edges_means=[]
		between_edges_means=[]
		

		for node in g.nodes():
			
			if g.nodes[node]['Network']==system:
				within_edges=[g[node][neighbor]['weight'] for neighbor in g[node] if (g.nodes[neighbor]['Network']==system)& (str(g[node][neighbor]['weight'])!='nan')]
				between_edges=[g[node][neighbor]['weight'] for neighbor in g[node] if (g.nodes[neighbor]['Network']!=system)& (str(g[node][neighbor]['weight'])!='nan')] 
				
				# print "\nSystem %s:"%system
				# print "\nWithin Edges:"
				# print within_edges
				# print "\nBetween Edges:"
				# print between_edges

				if len(within_edges)>0:
					mean_within= numpy.mean(within_edges)
					within_edges_means.append(mean_within)
				else:
					mean_within=0
				
				if len(between_edges)>0:
					mean_between= numpy.mean(between_edges)
					between_edges_means.append(mean_between)
				else:
					mean_between=0



		System_dict[system]['within']=numpy.mean(within_edges_means)
		System_dict[system]['between']=numpy.mean(between_edges_means)
		System_dict[system]['SystemSegregation']=(System_dict[system]['within']-System_dict[system]['between'])/System_dict[system]['within']
		

		SS_dict[system]= System_dict[system]
		


		# print "\n%s Sub-Network Mean Z-transformed Within-System Edges: %s"%(system, str(System_dict[system]['within']))
		# # print within_edges_means
		# print "%s Sub-Network Mean Z-transformed Between-System Edges: %s"%(system, str(System_dict[system]['between']))
		# # print between_edges_means
		# print "\n%s Sub-Network System Segregation: %s"%(system, str(System_dict[system]['SystemSegregation']))

	col_list = ['Network', 'Within', 'Between', 'SystemSegregation']
	
	out_dict = OrderedDict({})
	for c in col_list:
		out_dict[c] = []

	for network, measures in SS_dict.items():
		out_dict['Network'].append(network)
		out_dict['Within'].append(measures['within'])
		out_dict['Between'].append(measures['between'])
		out_dict['SystemSegregation'].append(measures['SystemSegregation'])

	# print "\nOut dict:"
	# print out_dict
	
	out_df = pandas.DataFrame.from_dict(out_dict, orient='index').T

	print(out_df)
	
	return out_df, g

def Centrality(graph, atlas, GT_df=pandas.DataFrame()):



	Networks = atlas['Network'].unique()



	BC_dict = {}
	# EC_dict = {}
	DC_dict = {}

	betweeness_dict = networkx.algorithms.centrality.betweenness_centrality(graph, weight='weight')
	# eigen_dict = networkx.algorithms.centrality.eigenvector_centrality(graph, weight='weight')
	degree_dict =  networkx.algorithms.centrality.degree_centrality(graph)

	for net in Networks:
		BC_list = []
		# EC_list = []
		DC_list = []

		for node in graph.nodes():

			if graph.nodes[node]['Network']== net:
				BC = betweeness_dict[node]
				# EC = eigen_dict[node]
				DC = degree_dict[node]

				BC_list.append(BC)
				# EC_list.append(EC)
				DC_list.append(DC)

		mean_BC = numpy.mean(BC_list)
		# mean_EC = numpy.mean(EC_list)
		mean_DC = numpy.mean(DC_list)
		BC_dict[net]=mean_BC
		# EC_dict[net]=mean_EC
		DC_dict[net] = mean_DC


	
	GT_df['Betweenness_Centrality'] = GT_df['Network'].map(BC_dict)
	# GT_df['Eigenvector_centrality'] = GT_df['Network'].map(EC_dict)
	GT_df['Degree_centrality'] = GT_df['Network'].map(DC_dict)
	
	return GT_df


		# sub_g = Graph.subgraph(graph, attribute='Network', attribute_classes=[net])

def Efficiency(graph, GT_df=pandas.DataFrame()):


	global_efficiency = networkx.algorithms.global_efficiency(graph)
	local_efficiency = networkx.algorithms.local_efficiency(graph)

	GT_df['Global_Efficiency']= global_efficiency
	GT_df['Local_Efficiency']= local_efficiency


	# efficiency = networkx.algorithms.efficiency(graph)

	return GT_df

def Smallworld(graph, GT_df=pandas.DataFrame()):
	# Networks = atlas['Network'].unique()

	sigma = networkx.algorithms.smallworld.sigma(graph, niter=100, nrand=10, seed=None)

	if len(GT_df['Network']>0):
		GT_df["SW_Sigma"]= sigma

	return GT_df

def get_network_density(result_df, labels_df):

    labels_unique=list(pandas.Series(labels_df["Network"]).drop_duplicates().to_numpy())
    N=len(labels_unique)
    nw_density_mat=numpy.zeros((N,N))



    result_df = result_df.where(numpy.triu(numpy.ones(result_df.shape)).astype(numpy.bool))
    result_df = result_df.fillna(0)
    result_df[result_df != 0 ] = 1



    thresh_links = hv.Dataset((list(labels_df['Network_id']), list(labels_df['Network_id']), result_df),
                      ['node1', 'node2'], 'value').dframe()



    thresh_links = thresh_links[thresh_links['value'] != 0].reset_index(drop=True)
    thresh_links_nw = thresh_links.groupby(["node1", "node2"]).sum()[["value"]].astype(int).reset_index()


    all_nodes_nw = numpy.triu(numpy.ones(result_df.shape))
    sum_all_nodes = all_nodes_nw.sum().sum()


    links = hv.Dataset((list(labels_df['Network_id']), list(labels_df['Network_id']), all_nodes_nw),
                      ['node1', 'node2'], 'value').dframe()

    links = links[links['value'] != 0].reset_index(drop=True)
    links_nw = links.groupby(["node1", "node2"]).sum()[["value"]].astype(int).reset_index()
    links_nw['node1_Network'] = [labels_unique[i] for i in links_nw['node1']]
    links_nw['node2_Network'] = [labels_unique[i] for i in links_nw['node2']]

    merged_df = pandas.merge(links_nw, thresh_links_nw, left_on=["node1", "node2"], right_on=["node1", "node2"], suffixes=('_allConn', '_threshConn'))

    merged_df['density']=round((merged_df['value_threshConn']/merged_df['value_allConn'])*100, 2)

    for ind, row in merged_df.iterrows():
        nw_density_mat[row['node1']][row['node2']] = round((row['value_threshConn']/row['value_allConn'])*100, 2)

    nw_density_mat = nw_density_mat + nw_density_mat.T

    nw_density_mat[numpy.diag_indices_from(nw_density_mat)]/= 2.0

    return nw_density_mat,merged_df

if __name__=='__main__':

	SS = SystemSegregation(G)
