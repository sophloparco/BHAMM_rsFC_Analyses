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

def binarize_vector(criteria, outlier_value, confound_df, thresh_inclusive=False):
	"""
	Take a key from the outlier criteria dictionary and make a binarized vector of values 
	that satisfy the specified outlier_value for that criteria
	outlier_value should be a string, formatted accordingly: [opperator]_[value]
	choices for opperator: >, <, =

	"""
	l = outlier_value.split('_', 1)
	opperator= l[0]
	# print 'Opperator: %s' %opperator
	val = l[1]
	# print 'Val: %s' %val
	if opperator == '>':
		if thresh_inclusive==True:
			v = confound_df[criteria]>=(float(val))*1
		else:
			v = confound_df[criteria]>(float(val))*1

	if opperator == '<':
		if thresh_inclusive==True:
			v = confound_df[criteria]<=(float(val))*1
		else:
			v = confound_df[criteria]<(float(val))*1

	if opperator == '=':
			v = confound_df[criteria]==(float(val))*1

	return v 

def main(project, fmriprep_version, confound_df, criteria={}):
	"""
	Take a df of confound time-series for a specific functional image 
	criteria is a dictionary with keys = confounds (column headers in the df) and 
	values are thresholds to flag for outlier volumes -- e.g. if one would like to flag 
	all volumes where there is greater than .5 mm framewise_displacement:
	criteria = {'framewise_displacement':'>_.5'} (by default threshold values are inclusive)

	"""

	global project_path, preproc_folder
	project_path = project
	preproc_folder = os.path.join(project_path, "BIDS", "derivatives", fmriprep_version, "fmriprep")
	df=pandas.DataFrame()
	i=0
	for k,v in criteria.items():

		outlier_vect = binarize_vector(k, v, confound_df)

		
		if i==0:
			old_outlier_vect = outlier_vect
		else:
			new_outlier_vect = outlier_vect | old_outlier_vect
			old_outlier_vect = new_outlier_vect
		
		df[k] = confound_df[k]
		i+=1

	df['outlier_volumes'] = new_outlier_vect



	return df

if __name__=='__main__':

	#can take full path to csv file as input confounds df but not recommended

	file = sys.argv[1]

	df = pandas.read_csv(file)



	main(df)



