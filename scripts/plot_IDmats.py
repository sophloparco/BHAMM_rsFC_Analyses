import pandas
import ConnectivityMatrix
import glob
import os
import sys


def main():
    project_folder = 'Path/to/analysis/folder'
    denoise_types=['0pAROMAaggressive']
    yeo_network_res='7'
    atlas_type='6' 


    for denoise_type in denoise_types:
        PLS_dir = os.path.join(project_folder, 'PLS_Analyses_%s'%denoise_type)
        out_dir= os.path.join(PLS_dir, "Matrix_Plots")
        if os.path.isdir(out_dir)==False:
            os.mkdir(out_dir)
        mat_dir = os.path.join(PLS_dir, 'Schaefer-200parcel_YeoNetwork-7_centroid-3mmradius', 'Whole-Brain_Matrices')
        print("\nmat_dir:")
        print(mat_dir)
        f_list = glob.glob(os.path.join(mat_dir, '[1,3,5,7]*'))
        # print('\nfiles:')
        # print(f_list)
        for f in f_list:
            atlas=pandas.read_csv(PLS_dir+'/atlas-schaefer_space-MNI_parcels-200_voxelres-2mm_yeonetworks-%s_SubNets-%s.csv'%(yeo_network_res, atlas_type))

            matrix = pandas.read_csv(f, sep='\t')

            ID=f.rsplit('/', 1)[-1].split('_', 1)[0]
            # print('\nID:')
            # print(ID)
            fig = ConnectivityMatrix.plot(matrix, atlas, sort_rois=False, vmin=-.4, vmax=.9, title="ID:%s"%(str(ID)))
            fig.savefig(os.path.join(out_dir, "%s_mat.png"%(str(ID))), facecolor=fig.get_facecolor())


        



if __name__=="__main__":

    main()