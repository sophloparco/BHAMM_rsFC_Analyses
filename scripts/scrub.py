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

def get_confounds(ID, task):
	
	folder = os.path.join(preproc_folder, 'sub-%s' %(str(ID)), "func")
	confound_file= glob.glob("%s/*%s*regressors.tsv"%(folder,task))[0]
	confound_df = pandas.read_csv(os.path.join(folder, confound_file), 
		sep="\t")
	return confound_df

def Diff(li1, li2): 
	#simple function for getting the set difference of two lists
	li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2] 
	return li_dif

def unique(li1, li2):
	#give the unique values across two lists
	li_uniq = list(set(li1+li2))
	return li_uniq

def add_buffers(df, buffer_vols=0):
	"""
	take a vector of Trues and Falses (representing outlier volumes to remove from a timeseries)
	Add the number of volumes equal to buffer_vols parameter preceding and following these outlier volumes
	e.g. if buffer_vols = 1, 1 volume preceding and 1 following the outlier volume will be included for scrubbing
	returns the indices of volumes to scrub from the timeseries data as a list.
	
	"""
	 

	indices = df.index[df['outlier_volumes']==True].tolist()
	
	buffer_indices=[]
	for i in indices:
		#add buffers for each outlier index
		for b in range(1, buffer_vols+1):
			#add indices of volumes to buffer 1 by 1 
			buffer_indices.append(i-1)
			buffer_indices.append(i+1)

	#filter any indices that may fall outside the range of the timeseries df 
	buffers = [b for b in buffer_indices if ((b <= len(df.index)) & (b>=0))]

	#get only the unique values for original outliers and buffers 
	vols2scrub = unique(indices, buffers)

	return vols2scrub

def scrub(TS, outlier_vector, additional_volumes=0):
	global project_path, preproc_folder
	project_path = "path/to/project-folder-enclosing-BIDS"
	preproc_folder = os.path.join(project_path, "BIDS", "derivatives", "fmriprep20.2.0", "fmriprep")
	to_remove = add_buffers(outlier_vector, additional_volumes)


	
	TS_scrubbed = TS[~TS.index.isin(to_remove)]

	scrubbed_volume_vector = outlier_vector.index.isin(to_remove) 

	return TS_scrubbed, scrubbed_volume_vector

if __name__=='__main__':

	scrub()