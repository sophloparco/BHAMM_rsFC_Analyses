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
import nistats



def get(img, ROIs, roi_type = 'sphere', radius=1, confounds=None, save_ts=False, filename=''):
	from nilearn import plotting

	timeseries = extract_signal(img, ROIs, roi_type, radius, confounds)


	if (save==True) & (len(filename)>0):
		save(timeseries, filename)
	elif (save==True) & (len(filename)==0):
		print("\nWARNING: No filename provided to save timeseries data\nTimeseries will not be saved...")

	return timeseries


def extract_signal(img, ROIs, roi_type='sphere', radius=1, confounds=None):

	if roi_type == 'sphere':
		from nilearn.input_data import NiftiSpheresMasker

		masker = NiftiSpheresMasker(ROIs, radius=radius, allow_overlap=False, standardize=False, detrend=False, t_r=2)



	elif roi_type == 'parcel':
		from nilearn.input_data import NiftiLabelsMasker
		masker = NiftiLabelsMasker(labels_img=ROIs, standardize=False,
			low_pass=0.08, high_pass=0.003, t_r=2)

	# Here we go from nifti files to the signal time series in a numpy
	# array. Note how we give confounds to be regressed out during signal
	# extraction
	timeseries = masker.fit_transform(img)
		
									
	return timeseries

def save(timeseries, filename=""):


	print("\nSaving timeseries to file:")
	print(filename)
	
	timeseries.to_csv(filename, mode='w', header=True,index=False, na_rep="NaN")

def to_conn_mat(timeseries, confound_file=None, conn_measure='correlation'):

	from nilearn.connectome import ConnectivityMeasure
	correlation_measure = ConnectivityMeasure(kind=conn_measure)
	
	if confound_file!=None:
		correlation_matrix = correlation_measure.fit_transform([timeseries], confound_file)[0]
	else: 
		correlation_matrix = correlation_measure.fit_transform([timeseries])[0]

	
	return correlation_matrix


