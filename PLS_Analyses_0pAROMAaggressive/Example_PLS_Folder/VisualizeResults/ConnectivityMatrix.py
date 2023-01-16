#Author: Sophia LoParco 
#Email: soph.loparco@gmail.com

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
import TimeSeries



def make(timeseries, confounds=None, conn_measure='correlation', mat_type='individual', save_mat=False, filename=''):

	
	matrix = to_conn_mat(timeseries, confounds=confounds, conn_measure=conn_measure, mat_type=mat_type)

	if (save_mat==True) & (len(filename)>0):
		save(matrix, filename)

	elif (save_mat==True) & (len(filename)==0):
		print("\nWARNING: No filename provided to save connectivity matrix\nConnectivity matrix will not be saved...")

	return matrix

def to_conn_mat(timeseries, confounds=None, conn_measure='correlation', mat_type='individual'):

	from nilearn.connectome import ConnectivityMeasure

	correlation_measure = ConnectivityMeasure(kind=conn_measure)


	if mat_type=="individual":
		if confounds!=None:
			matrix = correlation_measure.fit_transform([timeseries], confounds)[0]
		else: 
			matrix = correlation_measure.fit_transform([timeseries])[0]

	elif mat_type=='group':


		if confounds!=None:
			correlation_measure.fit_transform(timeseries, confounds)
		else: 
			correlation_measure.fit_transform(timeseries)

		matrix = correlation_measure.mean_


	return matrix

def save(matrix, filename=""):


	print("\nSaving matrix to file:")
	print(filename)
	
	matrix_df= pandas.DataFrame(data=matrix[:,:], # values
				index=range(1,len(matrix[0,:])+1),    # 1st column as index
				columns=range(1,len(matrix[0,:])+1))

	matrix_df.to_csv(filename, mode='w', header=True,index=False, na_rep="NaN")

def reorder_corrs(corrmtx, atlas_df):
	"""
	reorder correlation matrix according to network labels
	"""

	idx=numpy.argsort(atlas_df['Network'])
	tmp=corrmtx[:,idx]
	return(tmp[idx,:],atlas_df.iloc[idx,:])

def plot_reordered_corrs(corrmtx, atlas_df, sort_rois, ax=None, title="Set a title you dingus!"):
	"""
	plot correlation matrix after reordering
	"""
	
	if sort_rois == True:
		corr_reord, parceldata_reord=reorder_corrs(corrmtx, atlas_df)
	else:
		corr_reord = corrmtx
		parceldata_reord = atlas_df

	breaks=numpy.array([int(not i) for i in (parceldata_reord['Network'].values[:-1]==parceldata_reord['Network'].values[1:])])
	breaks[-1]= 1
	breaklocs=numpy.where(breaks)[0]

	for b in breaklocs:
		if b!=breaklocs[-1]:
			plt.plot([0,corrmtx.shape[0]-1],[b,b],color='white',linewidth=1.25)
			plt.plot([b,b],[0,corrmtx.shape[0]-1],color='white',linewidth=1.25)
	# find label locations
	# add a zero to help find label locations 
	breaklocs2=numpy.hstack(([0],breaklocs))
	label_locs=numpy.mean(numpy.vstack((breaklocs,breaklocs2[:-1])),0)
	networks=parceldata_reord['Network'].values[breaklocs]
	
	ax.set_yticks(label_locs)
	ax.set_yticklabels(networks, color='white', fontsize=14)
	ax.set_xticks(label_locs)
	ax.set_xticklabels(networks, rotation=90, color='white', fontsize=15)
	ax.set_title(title, color='white')

	
	return corr_reord, ax

def plot(matrix, atlas_df, sort_rois=True, vmax=0, vmin=0, cmap=matplotlib.cm.get_cmap('jet'), title="Add a title you dingus!"):
	from nilearn import plotting
	
	
	colorbar_position= [.85, .2, .005, .75]

		
	
	fig = plt.figure(figsize=(12,8))
	fig.patch.set_facecolor((0.5,0.5,0.5))
	fig.suptitle('%s'%(title), y=.98, size='x-large', fontweight='bold', color='white')
	ax= plt.subplot(111)


	if ((vmax==0)&(vmin==0)):
		vmax = numpy.max(numpy.abs(matrix.to_numpy()))+(1./8)*(numpy.max(numpy.abs(matrix.to_numpy())))
		vmin = 0

	title= ""
	
	im = plotting.plot_matrix(matrix, vmin=vmin,vmax=vmax, axes=ax, cmap=cmap, colorbar=False)

	reordered, ax= plot_reordered_corrs(matrix, atlas_df, sort_rois, ax=ax, title=title)

	norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
	cbar_ax= fig.add_axes(colorbar_position)
	cbar = plt.colorbar(im, cmap=cmap, cax=cbar_ax, fraction=0.1, pad=0.1)
	cbar.solids.set_edgecolor("face")
	plt.subplots_adjust(left=0.1, right=.99, bottom=0.2, top=0.95)

	return fig

def plot_rois_only(matrix, labels, sort_rois=True, vmax=0, vmin=0, cmap=matplotlib.cm.get_cmap('jet'), title="Add a title you dingus!"):
	from nilearn import plotting

	colorbar_position= [.85, .05, .005, .85]

		
	
	fig = plt.figure(figsize=(12,10))
	fig.patch.set_facecolor((0.5,0.5,0.5))
	fig.suptitle('%s'%(title), y=.98, size='x-large', fontweight='bold', color='white')
	ax= plt.subplot(111)

	

	if ((vmax==0)&(vmin==0)):
		vmax = numpy.max(numpy.abs(matrix.to_numpy()))+(1./8)*(numpy.max(numpy.abs(matrix.to_numpy())))
		vmin = 0

	title= ""
	
	
	label_locs = numpy.arange(0, len(labels)+1, step=1)
	ax.set_yticks(label_locs)
	ax.set_yticklabels(labels, color='white')
	ax.set_xticks(label_locs)
	ax.set_xticklabels(labels, rotation=75, color='white')


	im = plotting.plot_matrix(matrix, vmin=vmin,vmax=vmax, axes=ax, cmap=cmap, colorbar=False)
	plt.subplots_adjust(left=0.1, right=.9, bottom=0.15, top=0.95)
	norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
	cbar_ax= fig.add_axes(colorbar_position)
	cbar = plt.colorbar(im, cmap=cmap, cax=cbar_ax, fraction=0.1, pad=0.1)
	cbar.solids.set_edgecolor("face")

	return fig




