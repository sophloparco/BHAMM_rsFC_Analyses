AUTHOR: Sophia LoParco

EMAIL: soph.loparco@gmail.com

1) Setup python virtual environment by following the instructions in README_python_environment_setup.txt

2) create a to_denoise.csv file with the list of participants to denoise listed in a single column with the header "ID"

3) Denoise fmriprep outputs:
    Summary: This step will take images and confounds produced by fmriprep's initial preprocessing and perform denoising. 
    Outputs: Denoised Images will be made in a folder titled Denoised_Images

    i. Open the scripts/Denoise.py script: 
        Edit path variables initialized under line 648: "if __name__=='__main__':"

        project_folder --> path to folder containing your BIDS folder
        fmriprpe_version --> the version of fmriprep used (e.g., fmriprep20.2.0)
        analysis_folder --> the path to this analysis workspace

4) Create the Correlation Matrices: 
    Summary: Once the nifti images are denoised you can create the correlation matrices for your participants. 
    Outputs: Matrices that will be used by PLS connectivity and thus must be saved as [participant ID].txt files in the tab-delimited format 
 
    i. Open the scripts/wf.py script:
        Edit the path variables initialized under line 99: "if __name__=='__main__':"

        project_folder --> path to folder containing your BIDS folder
        fmriprpe_version --> the version of fmriprep used (e.g., fmriprep20.2.0)
        analysis_folder --> the path to this analysis workspace

    ii. This script wf.py will extract time-series from the denoised images in step 3 based on pre-specified centroids in the Schaefer 2018 atlas (see the previously used example atlas file: atlas-schaefer_space-MNI_parcels-200_voxelres-2mm_yeonetworks-7_SubNets-6.csv)
    
        You may select different options such as the number of centroids, radius of the centroids, etc. -- the variables for which are initialized in the script just below the path variables
    
        Time-series by default are scrubbed based on motion outlier threshold variables (defaults -- framewise-displacement: 0.5; DVARS: 2)
    
        The script will create a directory titled PLS_Analyses_[DenoiseType] and save the correlation matrices within it


5) PLS Connectivity Setup: 
    Summary: PLS Connectivity will be run in the PLS_Analyses_[DenoiseType] directory created in the previous step. 
    This folder will need to include:

        i. The folders of different types of PLS's you will run
        ii. The atlas file you are using (should be in the workspace folder)
        iii. Your matrices created in the previous step
        iiii. A folder of time-series titled OutlierVolumes_FD-0.5_DVARS-2.0 created in the previous step

    Within a PLS Analysis folder (e.g., PLS_Age_20-60_peri_1-Group) there should be:

        i. The scripts folder: this folder should contain your matlab scripts for running PLS connectivity. mainPLS_rsWholeBrainHC.m is the main script for conducting the analyses

            -displayPLSresults.m
            -create_colourbar.m
            -mainPLS_rs_WholeBrainHC.m
            -saveOutputFiles.m
            -stackPLSdatamat.m

	ii. VisualResults folder: should contain vis_1.py for visualizing BSR matrices specified within the same enclosed folder and r scripts for plotting correlation profiles etc. 