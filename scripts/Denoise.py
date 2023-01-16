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

def Diff(li1, li2): 
	#simple function for getting the set difference of two lists
	li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2] 
	return li_dif

def get_IDdf():
	
	

	p_file= glob.glob(os.path.join(analysis_folder, 'PLS_Analyses', "to_denoise.csv"))[0]
	print(p_file)

	participants= pandas.read_csv(p_file, sep=",")

	return participants

def filter_IDs(ID_df):

	derivatives = derivatives_dir 
	
	ID_df['ID']= ID_df['ID'].astype('int')
	has_data = [ID for ID in ID_df['ID'].values if os.path.isdir(os.path.join(derivatives,("sub-%s"%(str(ID))), 'func'))]


	
	no_data = ID_df[~ID_df['ID'].isin(has_data)]

	not_updated_in_participant_file = [ID for ID in has_data if ID not in ID_df['ID'].values]
	ID_df_filtered= ID_df[ID_df['ID'].isin(has_data)|ID_df['ID'].isin(not_updated_in_participant_file)]
	if no_data.empty==False:
		print("\nWarning: The following participants are listed as preprocessed in the participants.tsv file\nhowever, they do not seem to have data in the BIDS/derivatives folder for the project.")
		missing_data=no_data
		print('Derivatives Folder: %s'%(derivatives))

		print(missing_data['ID'].values)
	else:
		print("\nAll participants listed as preprocessed in the participants.tsv file have data in the BIDS/derivatives folder.")

	if len(not_updated_in_participant_file)!=0:
		print("\nWarning: The following particpiants have data in the BIDS/derivatives folder\nhowever they do not seem to be listed as preprocessed in the participants.tsv file.")
		print(not_updated_in_participant_file)
		print("\nPlease update their 'preprocessed' values to TRUE in the participants.tsv file.")

	else:
		print("\nAll participants with data in the BIDS/derivatives folder are listed as preprocessed in the participants.tsv file. ")
	

	return ID_df_filtered

def get_tasks(ID):

	func_folder= os.path.join(derivatives_dir, 'sub-%s'%(str(ID)), 'func')

	if os.path.isdir(func_folder):
		files = os.listdir(func_folder)

	else:

		print("\nWarning: Could not find a func data folder for participant %s" %(str(ID)))
		print("The following folder does not exist:")
		print(func_folder)

		files=[]
		return files

	tasks=[]
	for file in files:

		task = [t.split("-", 1)[1] for t in file.split("_") if "task" in t ]

		if (len(task)>0):
			if task[0] not in tasks:
				tasks.append(task[0])


	full_task_list= ['easy1', 'easy2', 'easy3', 'easy4', 'hard1', 'hard2', 'hard3', 'hard4', 'rest']

	if len(tasks)<len(full_task_list):
		print("\nThe following tasks were found for participant %s:"%(str(ID)))
		print(tasks)

		print("\nParticipant %s is missing the following tasks:"%(str(ID)))
		missing = Diff(tasks, full_task_list)
		print(missing)

	return tasks

def get_IDs():
	args= sys.argv[1:]

	ID_df= get_IDdf()


	if len(args)>0:
		IDs=args
		ID_df = ID_df[ID_df['ID'].isin(IDs)]
		if ID_df.empty:
			print("\nError: None of the IDs provided exist in the participants.tsv file.\nAborting script.")
			sys.exit()

		print("\nInput ID's detected, CorMats will be made for the following participants:")
		QC_df_option = 'a'
	
	elif len(args)==0:
		print("\nNo input ID's detected")
		QC_df_option = 'w'
	

	df = filter_IDs(ID_df)
	print('\nCorMats will be made for all available participants:')
	return df['ID'].values

def get_runs(ID, task, space_label, USE_AROMA, denoise, params, regressors):
	task_label = task#'run-%s' %(str(run))
	# space_label = 'MNI152NLin2009cAsym'
	space_label = 'MNI152NLin2009cAsym_desc'
	if USE_AROMA:
		denoise = 'aggressive'
		descriptor = 'preproc'
	else: 
		descriptor = 'preproc'
	mask = 'brain_mask'
	folder = os.path.join(derivatives_dir, 'sub-%s' %(str(ID)), "func")

	print("\nLooking for data in participant %s's functional folder: "%(str(ID)))
	print(folder)
	
	file_list= glob.glob("%s/*%s*%s*%s*.nii.gz"%(folder,task, space_label, descriptor))

	if len(file_list)==0:
		print("\nWarning: Unable to locate image for participant %s"%(str(ID)))
		img_file="None"
	else:
		img_file=file_list[0]

	print("\nCurrent Image:")

	print(img_file)

	file_list = glob.glob("%s/*%s*timeseries.tsv"%(folder,task))

	if len(file_list)==0:
		print("\nWarning: Unable to locate confound file for participant %s"%(str(ID)))
		confound_file="None"

	else:
		confound_file=file_list[0]


	print("\nConfound File: ")
	print(confound_file)

	if confound_file!='None':
		confound_df = pandas.read_csv(os.path.join(folder, confound_file), 
			sep="\t")
	else: 
		print("\nAborting script, due to missing data for participant %s. \nPlease verify at %s"%(str(ID), folder))
		sys.exit()
	if descriptor != "smoothAROMAnonaggr":
		

		
		if params!="AllConfounds":

			if denoise=="aggressive":
				compsfile = glob.glob("%s/*%s*AROMAnoiseICs.csv"%(folder,task))[0]
				comps_list= list(pandas.read_csv(os.path.join(folder, compsfile), sep=",").columns.values)
				for comp in comps_list: 
					regressors.append("aroma_motion_%s" %(comp.zfill(2)))

		elif params=="AllConfounds": 
			regressors= [c for c in confound_df.columns.values if not("aroma" in c)]
		elif (params=="AllConfounds"&USE_AROMA==True):
			regressors= confound_df.columns.values

		print('\nTask: %s; Regressors:\n' %(str(task)))
		if ((params=="0p")&(USE_AROMA==False)):
			empty_list=[]
			print(empty_list)
		else:
			print(regressors)

		confounds= confound_df[regressors].copy()

		

		confounds= confounds.fillna(value=0)
		
	else:
		regressors=["csf","white_matter", "global_signal"]
		confounds= None



	from nilearn import image
	img = nilearn.image.load_img(img_file)

	mask_file= glob.glob("%s/*%s*%s*%s*.nii.gz"%(folder,task, space_label, mask))[0]
	mask_img = nilearn.image.load_img(mask_file)


	return img, confounds, mask_img, os.path.join(folder, confound_file)

def get_surfs(ID, task):
	task_label = task
	space_label = 'fsaverage5'


	folder = os.path.join(dir_path, 'fmriprep_out', 'withFS', 'sub-%s' %(str(ID)), "func")


	right_img_file= glob.glob("%s/*%s*%s*hemi-R.func.gii"%(folder,task, space_label))

	left_img_file= glob.glob("%s/*%s*%s*hemi-L.func.gii"%(folder,task, space_label))




	print("\nCurrent Images:")

	print(left_img_file)
	print(right_img_file)



	surfaces = {'func_left':left_img_file, 'func_right':right_img_file}

	return surfaces

def get_anat(ID, space="MNI152NLin2009cAsym"):

	anat_folder = os.path.join(dir_path, 'fmriprep_out', 'withFS', 'sub-%s' %(str(ID)), "anat")

	space_label = space

	descriptor = "preproc_T1w"

	anat_file= glob.glob("%s/*%s*%s*.nii.gz"%(anat_folder, space_label, descriptor))
	anat_img= nilearn.image.load_img(anat_file)

	descriptor= "brain_mask"
	mask_file= glob.glob("%s/*%s*%s*.nii.gz"%(anat_folder, space_label, descriptor))
	mask_img= nilearn.image.load_img(mask_file)

	img = nilearn.image.math_img("img1 * img2", img1=anat_img, img2=mask_img)

	return img

def get_base_regs(params, USE_AROMA):

	if params=="28p":
		base_regs=["csf","white_matter", "global_signal","trans_x", "trans_y", "trans_z",
			 "rot_x", "rot_y", "rot_z", "framewise_displacement", "trans_x_derivative1","trans_y_derivative1", 
			 "trans_z_derivative1", "trans_x_derivative1_power2", "trans_y_derivative1_power2", "trans_z_derivative1_power2",
			 "rot_x_derivative1", "rot_y_derivative1", "rot_z_derivative1", "rot_x_derivative1_power2", "rot_y_derivative1_power2",
			 "rot_z_derivative1_power2", "rot_x_power2", "rot_y_power2", "rot_z_power2", "trans_x_power2", "trans_y_power2",
			 "trans_z_power2"]
	elif params=='36p':
		base_regs=["csf", 'csf_derivative1', 'csf_derivative1_power2', 'csf_power2',"white_matter", 
			'white_matter_derivative1', 'white_matter_derivative1_power2',
			'white_matter_power2', "global_signal", 'global_signal_derivative1', 
			'global_signal_derivative1_power2', 'global_signal_power2',"trans_x", "trans_y", 
			"trans_z", "rot_x", "rot_y", "rot_z", "trans_x_derivative1","trans_y_derivative1", 
			"trans_z_derivative1", "trans_x_derivative1_power2", "trans_y_derivative1_power2", 
			"trans_z_derivative1_power2", "rot_x_derivative1", "rot_y_derivative1", "rot_z_derivative1", 
			"rot_x_derivative1_power2", "rot_y_derivative1_power2","rot_z_derivative1_power2", "rot_x_power2", 
			"rot_y_power2", "rot_z_power2", "trans_x_power2", "trans_y_power2", "trans_z_power2"]
	elif params=='32p':
		base_regs=["csf", 'csf_derivative1', 'csf_derivative1_power2', 'csf_power2',"white_matter", 
			'white_matter_derivative1', 'white_matter_derivative1_power2',
			'white_matter_power2', "trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z", 
			"trans_x_derivative1","trans_y_derivative1", 
			"trans_z_derivative1", "trans_x_derivative1_power2", "trans_y_derivative1_power2", 
			"trans_z_derivative1_power2", "rot_x_derivative1", "rot_y_derivative1", "rot_z_derivative1", 
			"rot_x_derivative1_power2", "rot_y_derivative1_power2","rot_z_derivative1_power2", "rot_x_power2", 
			"rot_y_power2", "rot_z_power2", "trans_x_power2", "trans_y_power2", "trans_z_power2"]


	elif params=="6p":
		base_regs=["trans_x", "trans_y", "trans_z",
			 "rot_x", "rot_y", "rot_z"]
	
	elif params=='7p':
		base_regs=["trans_x", "trans_y", "trans_z",
			 "rot_x", "rot_y", "rot_z", "global_signal"]

	elif params=='8p':
		base_regs=["trans_x", "trans_y", "trans_z",
			 "rot_x", "rot_y", "rot_z", "csf", "white_matter"]
	
	elif params=='9p':
		base_regs=["trans_x", "trans_y", "trans_z",
			 "rot_x", "rot_y", "rot_z", "csf", "white_matter", "global_signal"]
	

	elif params=="24p":
		base_regs=["trans_x", "trans_y", "trans_z","rot_x", "rot_y", "rot_z", "trans_x_derivative1","trans_y_derivative1", 
			 "trans_z_derivative1", "trans_x_derivative1_power2", "trans_y_derivative1_power2", "trans_z_derivative1_power2",
			 "rot_x_derivative1", "rot_y_derivative1", "rot_z_derivative1", "rot_x_derivative1_power2", "rot_y_derivative1_power2",
			 "rot_z_derivative1_power2", "rot_x_power2", "rot_y_power2", "rot_z_power2", "trans_x_power2", "trans_y_power2",
			 "trans_z_power2"]
	elif params=="3p":
		base_regs=["csf","white_matter", "global_signal"]
	elif params=="0p":
		base_regs=[]
	elif params=="AllConfounds":
		base_regs=[]

	return base_regs

def clean_4Dimg(img, confounds, mask_img, params=None, USE_AROMA=False, smooth=True, smooth_fwhm=6, bandpass=(0.003, 0.08)):

	from nilearn.image import concat_imgs, index_img
	sub_imgs=list()

	for i in range(5, img.shape[3]):
		if (i<=confounds.shape[0]):
			sub_imgs.append(index_img(img, i))
		else: 
			break


	print("\n Concatenating %s 3D Images... " %(str(len(sub_imgs))))
	trimmed_img= concat_imgs(sub_imgs)
	print("\nConfound TimeSeries Length:\n")
	print(len(confounds.index.values))
	print("\nTrimmed Confound TimeSeries (first 5 rows removed):\n")
	trimmed_confounds= confounds.iloc[5:]
	print(len(trimmed_confounds.index.values))


	if smooth==True:
		in_img= nilearn.image.smooth_img(trimmed_img, fwhm=smooth_fwhm)
	else:
		in_img= trimmed_img

	#highpass-lowpass frequencies respectively for our bandpass filter of the signal
	min_freq, max_freq = bandpass

	if params!="0p":
		clean_img= nilearn.image.clean_img(in_img, confounds=trimmed_confounds.values, mask_img=mask_img, high_pass=min_freq, low_pass=max_freq, t_r=2.0)
	elif ((params=="0p")&(USE_AROMA==True)): 
		clean_img= nilearn.image.clean_img(in_img, confounds=trimmed_confounds.values, mask_img=mask_img, high_pass=min_freq, low_pass=max_freq, t_r=2.0)
	else:
		clean_img= nilearn.image.clean_img(in_img, mask_img=mask_img, high_pass=min_freq, low_pass=max_freq, t_r=2.0)


	

	
	print("\n Length of cleaned image: %s\n" %(str(clean_img.shape[3])))



	return clean_img

def plot_FD(ID,confounds):
	import matplotlib.pyplot as plt
	
	t= confounds.index
	FD= confounds['framewise_displacement']
	fig, ax = plt.subplots()
	ax.plot(t, FD)
	ax.grid()
	ax.set(xlabel='Frame (TR= 2 seconds)', ylabel='FD (mm)',
       title='%s Framewise Displacement'%(str(ID)))
	return fig

def plot_QCFC(QC_dict, denoise_type=""):
	import matplotlib.pyplot as plt
	
	fig, ax = plt.subplots()

	FCs= QC_dict['mFC']
	FDs= QC_dict['rmsFD']

	from scipy.stats.stats import pearsonr
	cor= pearsonr(FCs, FDs)
	

	ax.scatter(FDs, FCs)
	ax.set(xlabel='rms FD(individual)', ylabel='mean FC(individual)',
       title='FD ~ FC Correlation: %s'%(str(cor)))

	return fig

def get_meanConn(connectivity_matrix):

	CM=connectivity_matrix
	conn_dict= {}
	mask = numpy.ones(CM.shape,dtype='bool')

	mask[numpy.tril_indices(len(CM))] = False
	numpy.fill_diagonal(mask, val=True)
	
	ltri= CM.mask(mask)
	new_df= ltri.fillna(0)

	nonzero_array= new_df.to_numpy()[numpy.nonzero(new_df.to_numpy())]
	print("\nnon-zero array:\n")
	print(nonzero_array)
	
	mean= numpy.sum(nonzero_array)/len(nonzero_array)#numpy.sum(numpy.absolute(new_df.to_numpy()))/(len(new_df.values)*len(new_df.values)/2)
	return mean
	
def get_TimeSeries(ID, atlas, atlas_type, img):


	if atlas_type== "Volumetric Parcellation":

		from nilearn.input_data import NiftiMapsMasker
		masker = NiftiMapsMasker(maps_img=atlas.maps, standardize=True)
		

		timeseries = masker.fit_transform(img)


		print("\nDimensions of %s timeseries:\n" %(atlas_type))
		print(timeseries.shape)

		return timeseries

	elif atlas_type== "ROI Coordinate Map":
		coords= atlas['ROIs']
		from nilearn.input_data import NiftiSpheresMasker
		masker = NiftiSpheresMasker(coords, radius=5, standardize=True,
		low_pass=0.08, high_pass=0.01, t_r=2)


		timeseries = masker.fit_transform(img)
									
		return timeseries

	elif atlas_type== "Surface Parcellation":
		from nilearn import surface

		surfaces = get_surfs(ID, task)

		surf_L = surface.load_surf_data(surfaces['func_left'][0])
		surf_R = surface.load_surf_data(surfaces['func_right'][0])

		
		hr_surf= nilearn.surface.vol_to_surf(clean_img, big_fsaverage['pial_left'])

		clean_surfaces= {
			'all':hr_surf
		}

		nilearn.signal.clean(timeseries_L, sessions=None, confounds=confounds, t_r=2.0)
		nilearn.signal.clean(timeseries_R, sessions=None, confounds=confounds, t_r=2.0)

		clean_surfaces['all']

def mkCorMat(ID, atlas, atlas_type, img=None, confound_file=None, conn_measure="correlation", save=False):

	timeseries = get_TimeSeries(ID, atlas, atlas_type, img)
	


	from nilearn.connectome import ConnectivityMeasure
	correlation_measure = ConnectivityMeasure(kind=conn_measure)
	
	if confound_file!=None:
		correlation_matrix = correlation_measure.fit_transform([timeseries], confound_file)[0]
	else: 
		correlation_matrix = correlation_measure.fit_transform([timeseries])[0]

	
	return correlation_matrix, timeseries

def mkCorMats(IDs, atlas, atlas_type, imgs=None, confound_file=None, conn_measure="correlation", save=False):

	TSs=[]
	for i in range(0,len(IDs)):
		ID= IDs[i]
		img= imgs[i] 
		timeseries = get_TimeSeries(ID, atlas, atlas_type, img)
		TSs.append(timeseries)



	from nilearn.connectome import ConnectivityMeasure
	correlation_measure = ConnectivityMeasure(kind=conn_measure)
	
	
	correlation_matrix = correlation_measure.fit_transform(TSs)

	mean_matrix = correlation_measure.mean_

	
	return mean_matrix

def plot_CM(CM_df, vmax=0, vmin=0, cmap=matplotlib.cm.get_cmap('jet'), title="Add a title you dingus!"):
	from nilearn import plotting
	

	fig = plt.figure(figsize=(10,8))
	fig.patch.set_facecolor((0.5,0.5,0.5))
	fig.suptitle('%s'%(title), y=.98, size='x-large', fontweight='bold', color='white')
	ax= plt.subplot(111)

	if ((vmax==0)&(vmin==0)):
		vmax = numpy.max(numpy.abs(CM_df.to_numpy()))+(1./8)*(numpy.max(numpy.abs(CM_df.to_numpy())))
		vmin = 0

	plotting.plot_matrix(CM_df, vmin=vmin,vmax=vmax, axes=ax, cmap=cmap,colorbar=True)
		
	return fig

def main(project_path, analysis_folder, fmriprep_version, task, params, input_IDs = [], space_label = 'MNI152NLin6Asym', USE_AROMA=False, save=False, smooth=True, smooth_fwhm=6, save_denoised_niftis=False, bandpass=(0.003,0.08)):
	"""
	This script will make denoised image files for each individual with 
	specified confounds to regress, smooth True/False, smooth_fwhm specified (default is 6), bandpass filtered (default is 0.003 - 0.08 Hz range)

	params = the parameter type for regressing confounds in the images -- '36p' is recommended 
	but specific assortments of regressors can be found in the get_base_regs() function of this script
	
	"""

	#Set up our global variables for absolute paths to input images/files

	global dir_path, today, project_folder, BIDS_folder, derivatives_dir

	dir_path= project_path
	today = date.today()
	
	project_folder = project_path
	BIDS_folder= os.path.join(project_folder, "BIDS")
	derivatives_dir = os.path.join(BIDS_folder, 'derivatives', fmriprep_version,'fmriprep')


	if len(input_IDs) == 0:
		IDs_raw = get_IDs()
		
		IDs = [ID for ID in IDs_raw if 'rest' in get_tasks(ID)]

		missing_rest=Diff(list(IDs_raw), list(IDs))

		if len(missing_rest)!=0:
			print('\nWarning, the following participants do not have resting state data:')
			print(missing_rest)


	else: 
		IDs = input_IDs
	
	print('\nData will be run on the following %s participants:' %(str(len(IDs))))
	print(IDs)


	################################################################################
	

	descriptor=fmriprep_version

	if smooth==True: 
		smooth_desc= "%sFWHM"%(str(smooth_fwhm))
	else: 
		smooth_desc= "unsmoothed"
	

	if USE_AROMA==True: 
		denoising= "AROMAaggressive"
	elif USE_AROMA==False:
		denoising= ""
		
	if USE_AROMA==False:
		denoise="non-aggressive"

	else: 
		denoise="aggressive"

	from nilearn import plotting 
	output_dir = os.path.join(analysis_folder, 'Denoised_Images')

	if save==True:
		
		if os.path.isdir(output_dir)==False:
			print("\nMaking output_dir")
			print(output_dir)
			os.mkdir(output_dir)
		
		job_dir = os.path.join(output_dir,"pre-%s_task-%s_smooth-%s_denoise-%s"%(descriptor, str(task), str(smooth_desc), params+denoising))
		if os.path.isdir(job_dir):
			now = datetime.now()
			job_dir= job_dir+"_%s" %(str(now.strftime("_time-%H-%M-%S")))
		

		os.mkdir(job_dir)


	QC_dict= {'ID':IDs}
	QC_dict['rmsFD']=[]
	QC_dict['mFC']=[]
	cleaned_imgs=[]

	for ID in IDs:

		regs= get_base_regs(params, USE_AROMA)
		img, confounds, mask_img, confound_file= get_runs(ID, task, space_label, USE_AROMA, denoise, params, regs)
		
 		
		clean_img= clean_4Dimg(img, confounds, mask_img=mask_img, params=params, USE_AROMA=USE_AROMA, smooth=smooth, smooth_fwhm=smooth_fwhm, bandpass=bandpass)

		#WILL SAVE THE DENOISED NIFTI IF WE WANT TO 
		if save_denoised_niftis==True:
			clean_img_filename= os.path.join(job_dir, 'sub-%s_task-%s_smooth-%s_denoise-%s%s_bold.nii.gz'%(str(ID), str(task), str(smooth_desc),str(params), denoising))
			clean_img.to_filename(clean_img_filename)

		#ALWAYS RECORD THE TIME-SERIRES OF THE CONFOUNDS WE REGRESSED OUT OF THE IMAGE
		regressed_confounds_filename=os.path.join(job_dir,'sub-%s_task-%s_smooth-%s_denoise-%s_regressedconfounds.csv'%(str(ID), str(task), str(smooth_desc), str(params)))
		confounds.to_csv(regressed_confounds_filename, index=False)

		cleaned_imgs.append(clean_img)

	return job_dir
	

if __name__=="__main__":
	
	project_folder = os.path.join('path', 'to', 'directory', 'enclosing', 'BIDS', 'Folder')
	fmriprep_version = 'fmriprep20.2.0'
	analysis_folder= os.path.join('path', 'to','workspace')


	input_IDs = sys.argv[1:]
	print("\nInput IDs:")
	print(input_IDs)



	task = 'rest'
	param_list= ['0p'] 
	make=True
	if make==True:
		for params in param_list:

			main(project_folder, analysis_folder, fmriprep_version, task, params, space_label = 'MNI152NLin2009cAsym_desc', USE_AROMA=True, save=True,  save_denoised_niftis=True)




