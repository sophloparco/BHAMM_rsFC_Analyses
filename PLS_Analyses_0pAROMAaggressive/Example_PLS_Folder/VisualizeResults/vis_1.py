import pandas
import numpy
import nilearn
import nistats
import os
import sys
import matplotlib.pyplot as plt
import matplotlib
import glob
import re
from datetime import date, datetime
import SchaeferAtlas
import ConnectivityMatrix
from holoviews import opts, dim
import holoviews as hv
import Graph
import networkx


def add_labels(corrmtx, atlas_df, sort_rois, ax=None, title="Set a title you dingus!"):
	"""
	plot correlation matrix after reordering
	"""
	label_locs= range(0, len(atlas_df['Network'].values))
	networks=atlas_df['Network'].values

	ax.set_yticks(label_locs)
	ax.set_yticklabels(networks, color='black', fontsize=20)
	ax.set_xticks(label_locs)
	ax.set_xticklabels(networks, rotation=75, color='black', visible=True, fontsize=20)
	ax.set_title(title, color='black')



	return corrmtx, ax

def get_atlas_file(PLS_version, PLS_dir):
	if PLS_version=="PLS_WBHC":
		atlas = pandas.read_csv(os.path.join(PLS_dir, 'atlas-schaefer_space-MNI_parcels-200_voxelres-2mm_yeonetworks-7_SubNets-6.csv'))

	return atlas

def plot_BSR(df, atlas, PLS_version, plot_dir, lv=''):

	fig = ConnectivityMatrix.plot(df.to_numpy(), atlas, vmax=3, vmin=-3, title="%s BSR Matrix"%lv)
	f = "%s_%s_BSRMat"%(PLS_version, lv)
	fig.savefig(os.path.join(plot_dir, "%s.png"%f), facecolor=fig.get_facecolor())

def plot_density(matrix, atlas_df, sort_rois=False, vmax=0, vmin=0, cmap=matplotlib.cm.get_cmap('jet'), title="Add a title you dingus!"):
	from nilearn import plotting
	
	
	colorbar_position= [.85, .05, .005, .85]

		
	fig = plt.figure(figsize=(16,12))
	fig.patch.set_facecolor((1.0,1.0,1.0))
	fig.suptitle('%s'%(title), y=.98, size='x-large', fontweight='bold', color='black')
	ax= plt.subplot(111)


	if ((vmax==0)&(vmin==0)):
		vmax = numpy.max(numpy.abs(matrix.to_numpy()))+(1./8)*(numpy.max(numpy.abs(matrix.to_numpy())))
		vmin = 0

	title= ""
	
	reordered, ax= add_labels(matrix, atlas_df, sort_rois, ax=ax, title=title)
	

	im = plotting.plot_matrix(matrix, vmin=vmin,vmax=vmax, axes=ax, cmap=cmap, colorbar=False)

	label_locs= range(0, len(atlas_df['Network'].values))
	networks=atlas_df['Network'].values
	
	ax.set_yticks(label_locs)
	ax.set_yticklabels(networks, color='black', fontsize=20)
	ax.set_xticks(label_locs)
	ax.set_xticklabels(networks, rotation=75, color='black', visible=True, fontsize=20)
	ax.set_title(title, color='black')


	plt.subplots_adjust(left=0.1, right=.99, bottom=0.2, top=0.9)
	norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
	cbar_ax= fig.add_axes(colorbar_position)
	cbar = plt.colorbar(im, cmap=cmap, cax=cbar_ax, fraction=0.1, pad=0.1)
	cbar.solids.set_edgecolor("face")
	cbar.ax.tick_params(labelsize=20)


	return fig

def get_pos_threshold_matrix(matrix, threshold=0.0):
    pos_threshold_mat = matrix[matrix>threshold]
    pos_threshold_mat = pos_threshold_mat.fillna(0)
    pos_threshold_mat = pos_threshold_mat.round(decimals=2)

    return pos_threshold_mat

def get_neg_threshold_matrix(matrix, threshold=0.0):
    neg_threshold_mat = matrix[matrix<-1*threshold]
    neg_threshold_mat = neg_threshold_mat.fillna(0)
    neg_threshold_mat = neg_threshold_mat.round(decimals=2)

    return neg_threshold_mat

def get_thresholded_matrix(matrix,  threshold=0.0):

    pos_thresh_mat = get_pos_threshold_matrix(matrix,threshold) 
    neg_thresh_mat = get_neg_threshold_matrix(matrix,threshold)


    return pos_thresh_mat, neg_thresh_mat

def get_network_density(result_df, atlas):

	Networks=list(pandas.Series(atlas["Network"]).drop_duplicates().to_numpy())
	print('Networks in atlas**:')
	print(Networks)

	N=len(Networks)
	nw_density_mat=numpy.zeros((N,N))

	G, thresh_irrelevant = Graph.make_graph(result_df, atlas)
	g_thresh = Graph.color_nodes(G, atlas)

	Net_dict= {}
	Density_dict={}
	i=0
	j=0
	for Network in Networks:
		Net_dict[Network]={}
		i = i%N
		j = j%N
		for Other_Network in Networks:
			
			if Network != Other_Network:
				sub_g = Graph.subgraph(g_thresh, attribute_classes=[Network, Other_Network])
				sub_g_1 = Graph.subgraph(g_thresh, attribute_classes=[Network])
				sub_g_2 = Graph.subgraph(g_thresh, attribute_classes=[Other_Network])

				print("sub_g nodes:")
				print(sub_g.nodes)
				print("sub_g.edges:")
				print(sub_g.edges)

				between_edges = [e for e in sub_g.edges if e not in (list(sub_g_1.edges) + list(sub_g_2.edges))] 

				# edges = networkx.edge_boundary(sub_g, sub_g_1.nodes, sub_g_2.nodes)
				print("between edges:")
				print(between_edges)

				
				m = len(between_edges)
				n1 = len(sub_g_1.nodes)
				print("N1 (%s):%s" %(Network, str(n1)))
				n2= len(sub_g_2.nodes)
				print("N2 (%s):%s" %(Other_Network, str(n2)))
				density = (m)/(n1*n2)

			else:
				sub_g = Graph.subgraph(g_thresh, attribute_classes=[Network])
				m = len(sub_g.edges)
				print("\nm: %s"%(str(m)))
				n = len(sub_g.nodes)
				print("\nn: %s"%(str(n)))

				density = (2.0*m)/(n*(n - 1))

			
			Net_dict[Network].update({Other_Network:density*100})
			nw_density_mat[i,j]= density*100


			j+=1
		i+=1
	
	nw_density_mat = numpy.array(nw_density_mat)

	max_density = numpy.max(nw_density_mat)

	density_df = pandas.DataFrame(data=nw_density_mat[:,:], # values
 			index=range(0,len(nw_density_mat[:,0])),    # 1st column as index
			columns=range(0,len(nw_density_mat[0,:])))
	


	return density_df, Net_dict, max_density

def main():
	#DEFINE DIRECTORIES
	PLS_dir = os.path.join('Path', 'to','PLS_Analyses_0pAROMAaggressive')
	wd = os.path.join(PLS_dir, 'PLS_Age_20-60_peri_1-Group','VisualizeResults' )
	plot_dir = os.path.join(wd, 'plots')

	#LOAD ATLAS
	atlas_name = "Schaefer"
	space="MNI"
	voxel_res="2mm"
	parcel_res="200"
	yeo_network_res="7"
	

	PLS_version="PLS_WBHC"

	atlas = get_atlas_file(PLS_version, PLS_dir)

	LV1=True
	LV2=False
	LV3=False
	LV4=False

	#BSR RESULT FILES WE WANT TO PLOT
	if LV1 ==True:
		PLS_LV1_file = '%s_result_rawBSR_LV1.csv'%PLS_version
	
	if LV2==True:
		PLS_LV2_file = '%s_result_rawBSR_LV2.csv'%PLS_version

	if LV3==True:
		PLS_LV3_file = '%s_result_rawBSR_LV3.csv'%PLS_version

	if LV4==True:
		PLS_LV4_file = '%s_result_rawBSR_LV4.csv'%PLS_version
	
	#LOAD THEM & PLOT THEM
	if LV1==True:
		df1 = pandas.read_csv(os.path.join(wd,PLS_LV1_file), header=None)

		plot_BSR(df1, atlas, PLS_version, plot_dir, lv='LV1')

	if LV2==True:
		df2 = pandas.read_csv(os.path.join(wd,PLS_LV2_file), header=None)

		plot_BSR(df2, atlas, PLS_version, plot_dir, lv='LV2')

	if LV3==True:
		df3 = pandas.read_csv(os.path.join(wd,PLS_LV3_file), header=None)

		plot_BSR(df3, atlas, PLS_version, plot_dir, lv='LV3')

	if LV4==True:
		df4 = pandas.read_csv(os.path.join(wd,PLS_LV4_file), header=None)

		plot_BSR(df4, atlas, PLS_version, plot_dir, lv='LV4')




	thresh_file = open(os.path.join(plot_dir, "%s_DensityThresholds.txt"%PLS_version), "w")


	#MAKE DENSITY MATS LV1
	if LV1==True:
		thresh = numpy.percentile(df1.abs().to_numpy(), 95)
		print("LV1 Top 5% Thresh:")
		print(thresh)

		thresh_file.write("LV1 BSR Tresh %f" %(thresh))


		pos_thresh_df, neg_thresh_df = get_thresholded_matrix(df1, thresh)

		#POS DENSITY MATRIX
		density_df_pos, density_dict_pos, max_density_pos = get_network_density(pos_thresh_df, atlas)

		#NEG DENSITY MATRIX
		density_df_neg, density_dict_neg, max_density_neg = get_network_density(neg_thresh_df, atlas)

		max_density= max([max_density_pos, max_density_neg])
		max_density=9
		#PLOT POS DENSITY MAT
		Networks=pandas.DataFrame({'Network':pandas.Series(atlas["Network"]).drop_duplicates().to_numpy()})

		print("Networks in Atlas: ")
		print(Networks)

		fig = plot_density(density_df_pos, Networks, vmax=max_density, vmin=0, cmap=matplotlib.cm.get_cmap('Reds'), title="Positive Density Matrix")
		f = "%s_resultLV1_PosDensityMat"%PLS_version
		fig.savefig(os.path.join(plot_dir, "%s.png"%f), facecolor=fig.get_facecolor())


		#PLOT NEG DENSITY MAT

		Networks=pandas.DataFrame({'Network':pandas.Series(atlas["Network"]).drop_duplicates().to_numpy()})

		print('Networks in atlas:')
		print(Networks)

		fig = plot_density(density_df_neg, Networks, vmax=max_density, vmin=0, cmap=matplotlib.cm.get_cmap('Blues') ,title="Negative Density Matrix")
		f = "%s_resultLV1_NegDensityMat"%PLS_version
		fig.savefig(os.path.join(plot_dir, "%s.png"%f), facecolor=fig.get_facecolor())


	#MAKE DENSITY MATS LV2
	if LV2==True:
		thresh = numpy.percentile(df2.abs().to_numpy(), 95)
		print("LV2 Top 5% Thresh:")
		print(thresh)

		thresh_file.write("LV2 BSR Tresh %f\n" %(thresh))
		

		pos_thresh_df, neg_thresh_df = get_thresholded_matrix(df2, thresh)

		#POS DENSITY MATRIX
		density_df_pos, density_dict_pos, max_density = get_network_density(pos_thresh_df, atlas)

		#NEG DENSITY MATRIX
		density_df_neg, density_dict_neg, max_density_neg = get_network_density(neg_thresh_df, atlas)

		max_density= max([max_density_pos, max_density_neg])

		#PLOT POS DENSITY MATRIX

		Networks=pandas.DataFrame({'Network':pandas.Series(atlas["Network"]).drop_duplicates().to_numpy()})

		print(Networks)

		fig = plot_density(density_df_pos, Networks, vmax=max_density, vmin=0, cmap=matplotlib.cm.get_cmap('Reds'), title="Positive Density Matrix")
		f = "%s_resultLV2_PosDensityMat"%PLS_version
		fig.savefig(os.path.join(plot_dir, "%s.png"%f), facecolor=fig.get_facecolor())

		#PLOT NEG DENSITY MATRIX
		Networks=pandas.DataFrame({'Network':pandas.Series(atlas["Network"]).drop_duplicates().to_numpy()})

		print(Networks)

		fig = plot_density(density_df_neg, Networks, vmax=max_density, vmin=0, cmap=matplotlib.cm.get_cmap('Blues'), title="Negative Density Matrix")
		f = "%s_resultLV2_NegDensityMat"%PLS_version
		fig.savefig(os.path.join(plot_dir, "%s.png"%f), facecolor=fig.get_facecolor())
	
	#MAKE DENSITY MATS LV3
	if LV3==True:
		thresh = numpy.percentile(df3.abs().to_numpy(), 95)
		print("LV3 Top 5% Thresh:")
		print(thresh)

		thresh_file.write("LV3 BSR Tresh %f\n" %(thresh))
		

		pos_thresh_df, neg_thresh_df = get_thresholded_matrix(df3, thresh)

		#POS DENSITY MATRIX
		density_df, density_dict, max_density = get_network_density(pos_thresh_df, atlas)


		Networks=pandas.DataFrame({'Network':pandas.Series(atlas["Network"]).drop_duplicates().to_numpy()})

		print(Networks)

		fig = plot_density(density_df, Networks, vmax=max_density, vmin=0, cmap=matplotlib.cm.get_cmap('Reds'), title="Positive Density Matrix")
		f = "%s_resultLV3_PosDensityMat"%PLS_version
		fig.savefig(os.path.join(plot_dir, "%s.png"%f), facecolor=fig.get_facecolor())


		#NEG DENSITY MATRIX
		density_df, density_dict, max_density = get_network_density(neg_thresh_df, atlas)


		Networks=pandas.DataFrame({'Network':pandas.Series(atlas["Network"]).drop_duplicates().to_numpy()})

		print(Networks)

		fig = plot_density(density_df, Networks, vmax=max_density, vmin=0, cmap=matplotlib.cm.get_cmap('Blues'), title="Negative Density Matrix")
		f = "%s_resultLV3_NegDensityMat"%PLS_version
		fig.savefig(os.path.join(plot_dir, "%s.png"%f), facecolor=fig.get_facecolor())

	#MAKE DENSITY MATS LV4
	if LV4==True:
		thresh = numpy.percentile(df4.abs().to_numpy(), 95)
		print("LV4 Top 5% Thresh:")
		print(thresh)

		thresh_file.write("LV4 BSR Tresh %f\n" %(thresh))
		

		pos_thresh_df, neg_thresh_df = get_thresholded_matrix(df4, thresh)

		#POS DENSITY MATRIX
		density_df, density_dict, max_density = get_network_density(pos_thresh_df, atlas)


		Networks=pandas.DataFrame({'Network':pandas.Series(atlas["Network"]).drop_duplicates().to_numpy()})

		print(Networks)

		fig = plot_density(density_df, Networks, vmax=max_density, vmin=0, cmap=matplotlib.cm.get_cmap('Reds'), title="Positive Density Matrix")
		f = "%s_resultLV4_PosDensityMat"%PLS_version
		fig.savefig(os.path.join(plot_dir, "%s.png"%f), facecolor=fig.get_facecolor())


		#NEG DENSITY MATRIX
		density_df, density_dict, max_density = get_network_density(neg_thresh_df, atlas)


		Networks=pandas.DataFrame({'Network':pandas.Series(atlas["Network"]).drop_duplicates().to_numpy()})

		print(Networks)

		fig = plot_density(density_df, Networks, vmax=max_density, vmin=0, cmap=matplotlib.cm.get_cmap('Blues'), title="Negative Density Matrix")
		f = "%s_resultLV4_NegDensityMat"%PLS_version
		fig.savefig(os.path.join(plot_dir, "%s.png"%f), facecolor=fig.get_facecolor())

	thresh_file.close()

if __name__ == '__main__':
	

	main()