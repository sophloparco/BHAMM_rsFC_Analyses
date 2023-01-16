import pandas
import numpy
import nilearn
from nilearn import image
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
import ConnectivityMatrix
import find_outlier_vols
import scrub
import Graph
import GT
import networkx

def get_confounds(ID, task, preproc_folder):
	print("ID: ")
	print(ID)
	
	folder = os.path.join(preproc_folder, 'sub-%s' %(str(ID)), "func")

	file_pattern = "%s/sub-%s_task-%s_run-*_desc-confounds_timeseries.tsv"%(folder,ID,task)

	print("\nLooking for confound file:")
	print(file_pattern)

	confound_file= glob.glob(file_pattern)[0]

	confound_df = pandas.read_csv(os.path.join(folder, confound_file), 
		sep="\t")
	return confound_df

def get_group(group_key, denoise_type):
	analysis_folder = 'path/to/analysis_folder'

	file_pattern = os.path.join(analysis_folder, "Denoised_Images/pre-fmriprep20.2.0_task-rest_smooth-6FWHM_denoise-%s/sub-%s_task-rest_smooth-6FWHM_denoise-%s_bold.nii.gz"%(denoise_type, group_key, denoise_type))

	print("\nLooking for files that match the following pattern:")
	print(file_pattern)

	files = glob.glob(file_pattern)

	files.sort()

	print("\nThe following denoised image data files will be included in the mean matrix:")
	print(files)



	IDs=[]
	img_files=[]
	for f in files:


		print("\nf:")
		print(f)

		f_name = f.rsplit('/', 1)[1]
		print("\nf_name:")
		print(f_name)

		ID = f_name.rsplit('_')[0].split('-')[1]
		print("\nID:")
		print(ID)

		IDs.append(ID)
		img_files.append(f)


	data= {'ID':IDs, 'File':img_files}
	

	return data

if __name__=='__main__':

	fmriprep_version = 'fmriprep20.2.0'

	project_folder = 'path/to/analysis_folder'
	preproc_folder = os.path.join('path', 'to','BIDS', 'derivatives',fmriprep_version, 'fmriprep')

	task='rest'

	atlas_name = "Schaefer"
	space="MNI"
	voxel_res="2mm"
	parcel_res_list=["200"]
	yeo_networks=['7']
	density = .04

	FD_Thresh = 0.5
	DVARS_Thresh = 2.0

	radii=[3]
	scrubbing=True
	scrub_criteria = {'framewise_displacement': '>_%s'%(str(FD_Thresh)), 'std_dvars':'>_%s'%(str(DVARS_Thresh))}
	denoise_types=['0pAROMAaggressive']

	for denoise_type in denoise_types:
		PLS_dir = os.path.join(project_folder, 'PLS_Analyses_%s'%denoise_type)
		if os.path.isdir(PLS_dir)==False:
			os.mkdir(PLS_dir)

		OutlierVol_dir=os.path.join(PLS_dir, 'OutlierVolumes_FD-%s_DVARS-%s'%(str(FD_Thresh),str(DVARS_Thresh)))


		if os.path.isdir(OutlierVol_dir)==False:
			os.mkdir(OutlierVol_dir)

		for yeo_network_res in yeo_networks: 
			for parcel_res in parcel_res_list:
				for radius in radii:

					out_dir = os.path.join(PLS_dir, "Schaefer-%sparcel_YeoNetwork-%s_centroid-%smmradius"%(parcel_res, yeo_network_res, str(radius)))

					if os.path.isdir(out_dir)==False:
							os.mkdir(out_dir)


					groups = {"Matrices":"[1,3,5,7]*"}
						submat_dict={'6':'Whole-Brain'}

					for group,key in groups.items():
						data = get_group(key, denoise_type)
						print("\nData:")
						print(data)
						for atlas_type, desc in submat_dict.items(): 

							IDmat_out_dir = os.path.join(out_dir, "%s_Matrices"%(desc))
							if os.path.isdir(IDmat_out_dir)==False:
								os.mkdir(IDmat_out_dir) 


							atlas = pandas.read_csv(project_folder+'/atlas-schaefer_space-MNI_parcels-200_voxelres-2mm_yeonetworks-%s_SubNets-%s.csv'%(yeo_network_res, atlas_type))
							print(atlas)

							coords = SchaeferAtlas.get_network_coords(atlas)


							TS_list = []
							for i in range(0, len(data['ID'])):
								ID = data['ID'][i]
								img_file = data['File'][i]

								print("Finding file: %s" %img_file)
								img = nilearn.image.load_img(img_file)

								print("Making TimeSeries from: %s's image:"%ID)

								timeseries = TimeSeries.get(img, coords, radius=radius)
								
								confound_df = get_confounds(ID, task, preproc_folder)

								trimmed_confound_df = confound_df.iloc[5:]
								outlier_df = find_outlier_vols.main(project_folder, fmriprep_version, trimmed_confound_df, scrub_criteria)

								outlier_vector = pandas.DataFrame({'outlier_volumes':outlier_df['outlier_volumes']})

								print("\nTime Series:")
								print(timeseries)

								ts_df= pandas.DataFrame(data=timeseries[:,:], # values
								index=range(0,len(timeseries[:,0])),    # 1st column as index
								columns=range(0,len(timeseries[0,:])))
							
								
								scrubbed_ts_df, WithBuffer = scrub.scrub(ts_df, outlier_vector, additional_volumes=1)
								
								if scrubbing==True:
									outlier_df['WithBuffer'] = WithBuffer
									outlier_df.to_csv(os.path.join(OutlierVol_dir, 'scrubbed_FD-%s_DVARS-%s_volumes_%s.csv'%(str(FD_Thresh),str(DVARS_Thresh),str(ID))))
									scrubbed_ts_df.to_csv(os.path.join(OutlierVol_dir, 'scrubbed_FD-%s_DVARS-%s_timeseries_%s.csv'%(str(FD_Thresh),str(DVARS_Thresh),str(ID))))


											
								if len(scrubbed_ts_df.index)>=180:

									if scrubbing==True: 
										matrix = ConnectivityMatrix.make(scrubbed_ts_df.to_numpy(), mat_type='individual')
										
									else:
										matrix = ConnectivityMatrix.make(ts_df.to_numpy(), mat_type='individual')
										ts_df.to_csv(os.path.join(OutlierVol_dir, 'timeseries_%s.csv'%(str(ID))))

								
									matrix_df= pandas.DataFrame(data=matrix[:,:], # values
									index=range(0,len(matrix[:,0])),    # 1st column as index
									columns=range(0,len(matrix[0,:])))
									print("\nMatrix:")
									print(matrix_df)

									
									matrix_df.to_csv(os.path.join(IDmat_out_dir, "%s_mat.txt" %(str(ID))), header=False, index=False, sep='\t')
									#SAVE AS A TSV WITH THE .TXT ENDING TO BE COMPATIBLE WITH CONN_PLS MATLAB SCRIPTS

							 


							




