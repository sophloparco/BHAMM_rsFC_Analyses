import os
import sys
import glob2 as glob
import re
import pandas
import numpy
import math
from collections import OrderedDict
import nilearn
import nibabel


def get(space, voxel_res, parcel_res, yeo_network_res, centroid=True, img_affine=None):

	#THIS SERVES AS THE MAIN FUNCTION IN THE SCHAEFERATLAS.PY SCRIPT -- TAKES INPUT PARAMETERS TO GRAB THE APPROPRIATE FILE IN SCHAEFER'S FORMAT AND TRANSFORM IT TO FORMAT FOR THIS ANALYSIS

	atlas_dir= "path/to/atlas-Schaefer_space-MNI/"

	if (centroid==True)&(space=="MNI"): 
		atlas_file= os.path.join(atlas_dir, "Centroid_coordinates","Schaefer2018_%sParcels_%sNetworks_order_FSLMNI152_%s.Centroid_RAS.csv" %(parcel_res, yeo_network_res, voxel_res))


	atlas_df = pandas.read_csv(atlas_file)

	
	labeled_atlas_df = unpack_node_labels(atlas_df)

	std_df = color_networks(labeled_atlas_df, yeo_network_res)

	return std_df

def unpack_node_labels(in_atlas):

	#TRANSFORM THE INPUT ATLAS TO MYLE STANDARD FORMAT USED BY ANALYSIS WORKFLOW 
	
	out_atlas = in_atlas
	# gapminder['gdpPercap_ind'] = gapminder.gdpPercap.apply(lambda x: 1 if x >= 1000 else 0)


	out_atlas['Hemisphere'] = out_atlas["ROI Name"].apply(lambda x: "LH" if "LH" in x else "RH")
	out_atlas['Network'] = out_atlas["ROI Name"].apply(lambda x: x.split("_")[len(x.split("_"))-3] if len(x.split('_'))==5 else x.split("_")[len(x.split("_"))-2])
	out_atlas['Region'] = out_atlas["ROI Name"].apply(lambda x: x.split("_")[len(x.split("_"))-2])
	out_atlas['Member'] = out_atlas["ROI Name"].apply(lambda x: x.split("_")[len(x.split("_"))-1])
	out_atlas['Coords'] = list(zip(out_atlas['R'], out_atlas['A'], out_atlas['S']))
	out_atlas = out_atlas.sort_values(by='Network')
	out_atlas['Node'] = range(0,len(out_atlas['Region'].values))



	return out_atlas


def color_networks(atlas, yeo_network_res):

	if yeo_network_res == '7':
		color_dict = {'Default':'red', 'Cont':'yellow', 'DorsAttn':'green', 'SalVentAttn':'black', 'Limbic':'cyan', 'SomMot':'orange', 'Vis':'blue'}

		colors= {'Limbic':(0,255,255), 'SomMot':(255,128,0),
	 	'Default':(255,0,0), 'Vis':(0, 0, 255),
	 	'Cont':(255, 255, 0), 'SalVentAttn':(0, 0, 0), 
	 	'DorsAttn':(0, 204, 0)}

	elif yeo_network_res == '17':

		color_dict = {'DefaultA':'red', 'DefaultB':'red','DefaultC':'red','ContA':'yellow','ContB':'yellow','ContC':'yellow', 'DorsAttnA':'green','DorsAttnB':'green', 'SalVentAttnA':'black','SalVentAttnB':'black', 'LimbicA':'cyan','LimbicB':'cyan', 'TempPar':'purple', 'SomMotA':'orange','SomMotB':'orange', 'VisCent':'blue','VisPeri':'blue'}
	
		colors= {'LimbicA':(0,255,255),'LimbicB':(0,255,255), 
		'SomMotA':(255,128,0),'SomMotB':(255,128,0), 
		'TempPar':(127,0,255),
	 	'DefaultA':(255,0,0),'DefaultB':(255,0,0), 'DefaultC':(255,0,0), 
	 	'VisCent':(0, 0, 255),'VisPeri':(0, 0, 255),
	 	'ContA':(255, 255, 0), 'ContB':(255, 255, 0),'ContC':(255, 255, 0), 
	 	'SalVentAttnA':(0, 0, 0),'SalVentAttnB':(0, 0, 0), 
	 	'DorsAttnA':(0, 204, 0),'DorsAttnB':(0, 204, 0)}

	syst_color= {}
	for system in atlas['Network'].unique():
		r, g, b = colors[system]
		if r!=0:
			r1=float(r)/255
		else:
			r1=0
		if g!=0:
			g1=float(g)/255
		else:
			g1=0
		if b!=0:
			b1=float(b)/255
		else:
			b1=0
		syst_color[system]= (r1, g1, b1)

	atlas['Color'] = atlas['Network'].map(color_dict)
	atlas['ColorCode'] = atlas['Network'].map(syst_color)

	return atlas


def get_network_coords(atlas, network=None):

	if network==None:

		coords = list(atlas.apply(lambda x: (x['R'], x['A'], x['S']), axis=1))
	
	else: 
		network_rois = atlas[atlas['Network'].str.contains(network)]

		coords = list(network_rois.apply(lambda x: (x['R'], x['A'], x['S']), axis=1))



	return coords


def get_node_colors(atlas, network=None):

	if network==None:

		color_list = list(atlas['Color'])

	else:

		sub_network = atlas[atlas['Network'].str.contains(network)]
		color_list = list(sub_network['Color'])

	return color_list


def inspect_coords(atlas, network=None, coord_system='RAS'):
	from nilearn import plotting


	coords = get_network_coords(atlas, network=network)#list(network_rois.apply(lambda x: (x['R'], x['A'], x['S']), axis=1))
	color_list = get_node_colors(atlas, network=network)


	view = view_coords(coords, color_list, )


	return view

def view_coords(coords, colors):

	from nilearn import plotting

	view = plotting.view_markers(coords, colors, marker_size=5)

	return view

if __name__ =="__main__":


	space="MNI"
	voxel_res="2mm"
	parcel_res="200"
	yeo_network_res="7"


	df = get(space, voxel_res, parcel_res, yeo_network_res)

	df.to_csv(os.path.join("atlas-schaefer_space-%s_parcels-%s_voxelres-%s_yeonetworks-%s.csv"%(space, parcel_res, voxel_res, yeo_network_res)), index=False)
