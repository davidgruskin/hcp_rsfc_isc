#!/bin/python
import os
import numpy as np
import scipy.io as io
import pandas as pd
import nibabel as nb
import nilearn
import subprocess
import sys
import argparse

def main():
    # Set up parser
    parser = argparse.ArgumentParser(description='This script denoises HCP 7T data')
    parser.add_argument("-first_id", "--first_id", required=False, type=str, help="string indicating index of the first subject (lowest possible value = 1) on which this script will run")
    parser.add_argument("-last_id", "--last_id", required=False, type=str, help="string indicating index of the last subject (highest possible value = 184) on which this script will run")
    parser.add_argument("-tasktype", "--tasktype", required=False, type=str, help="string indicating whether to process REST or MOVIE runs")
    parser.add_argument("-startswith", "--startswith", required=False, type=str, help="string indicating subject IDs that start with X (e.g. '2' to run on all subjects whose ID starts with '2') on which this script will run")
    parser.add_argument("-proj_dir", "--proj_dir", required=False, default='/data/data7/HCP_7T', type=str, help="string indicating location of data")
    parser.add_argument("-regress_motion", "--regress_motion", required=True, type=int, help="binary indicating whether or not to perform spike regression")
    parser.add_argument("-wb", "--wb", required=False, default='/Applications/workbench/bin_macosx64/wb_command', type=str, help="string indicating path to wb_command")
    args = vars(parser.parse_args())

    # Set up subject list
    wb = args['wb']
    proj_dir = args['proj_dir']
    hcp_pheno = [x for x in os.listdir(proj_dir) if '.py' not in x]
    if args['first_id']:
        sub_list = hcp_pheno[int(args['first_id'])-1:int(args['last_id'])-1]
    elif args['startswith']:
        sub_list = [x for x in hcp_pheno if x.startswith(args['startswith'])]
    else:
        sub_list = hcp_pheno

    # Loop over subjects
    for sub in sub_list:
        print(sub)
        
        sub_dir = os.path.join(proj_dir, str(sub), 'MNINonLinear/Results')
        scan_list = ['rfMRI_REST1_7T_PA', 'rfMRI_REST2_7T_AP', 'rfMRI_REST3_7T_PA', 'rfMRI_REST4_7T_AP','tfMRI_MOVIE1_7T_AP', 'tfMRI_MOVIE2_7T_PA', 'tfMRI_MOVIE3_7T_PA', 'tfMRI_MOVIE4_7T_AP']
        if args['tasktype']:
            scan_list = [x for x in scan_list if args['tasktype'] in x]
        for scantype in scan_list:
            print('\t' + scantype)

            # Find data and make filenames
            cii_file = os.path.join(sub_dir, scantype, scantype + '_Atlas_hp2000_clean.dtseries.nii')
            fake_nii = cii_file.replace(".dtseries", "fake") 
            clean_cii = fake_nii.replace('fake.nii','_nilearn.dtseries.nii')
            if os.path.exists(cii_file):

                # Make a fake nifti file from cifti data 
                subprocess.call(['/Applications/workbench/bin_macosx64/wb_command', '-cifti-convert', '-to-nifti', cii_file, fake_nii])
                gsr_text = os.path.join(sub_dir, scantype, 'GSR.txt')
                gsr_cmd = '/gpfs/milgram/apps/hpc.rhel7/software/AFNI/2018.08.28/3dmaskave -quiet '
                gsr_cmd = gsr_cmd + cii_file + ' > ' + gsr_text

                # Calculate global signal using AFNI
                with open(gsr_text, 'w') as f:
                    subprocess.call(['/Users/admin/abin/3dmaskave', '-quiet', fake_nii], stdout=f)

                # Load AFNI global signal output (and high FD frames to be regressed out of data if relevant)
                gsr = pd.read_csv(gsr_text, header=None)
                gsr_deriv = pd.Series(np.gradient(gsr[0].values))
                motion_file = os.path.join(sub_dir,scantype,'motion_regressors1.csv')
                if os.path.exists(motion_file) and args['regress_motion'] == 1:
                    motion_regressors = pd.read_csv(motion_file, header=None)
                    confounds = pd.concat([gsr, gsr_deriv,motion_regressors], axis=1)
                else:
                    confounds = pd.concat([gsr, gsr_deriv], axis=1)
                dat = nilearn.image.load_img(fake_nii)

                # Run denoising
                clean_dat = nilearn.image.clean_img(dat, t_r=1.000, confounds=np.mat(confounds))

                # Export to new file
                clean_out = fake_nii.replace('fake.nii','_nilearn.nii')
                clean_dat.to_filename(clean_out)

                # Convert back to nifti
                combine_cii = 'wb_command -cifti-convert -from-nifti ' + clean_out + ' ' + cii_file + ' ' + clean_cii
                subprocess.call(['/Applications/workbench/bin_macosx64/wb_command', '-cifti-convert', '-from-nifti', clean_out, cii_file, clean_cii])
                
                # Delete intermediary files
                os.remove(fake_nii) 
                os.remove(clean_out)

if __name__ == "__main__":
    #sub = sys.argv[0]
    main()










