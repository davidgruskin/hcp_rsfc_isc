#!/bin/python
import os
import subprocess
import sys 
import argparse

def main():
    # sSet up parser
    parser = argparse.ArgumentParser(description='This script parcellates HCP 7T data')
    parser.add_argument("-first_id", "--first_id", required=False, type=str, help="string indicating index of the first subject (lowest possible value = 1) on which this script will run")
    parser.add_argument("-last_id", "--last_id", required=False, type=str, help="string indicating index of the last subject (highest possible value = 184) on which this script will run")
    parser.add_argument("-tasktype", "--tasktype", required=False, type=str, help="string indicating whether to process REST or MOVIE runs")
    parser.add_argument("-startswith", "--startswith", required=False, type=str, help="string indicating subject IDs that start with X (e.g. '2' to run on all subjects whose ID starts with '2') on which this script will run")
    parser.add_argument("-proj_dir", "--proj_dir", required=False, default='/data/data7/HCP_7T', type=str, help="string indicating location of data")
    parser.add_argument("-regress_motion", "--regress_motion", required=True,type=int,help="binary indicating whether or not to perform spike regression")
    parser.add_argument("-wb", "--wb", required=False, default='/Applications/workbench/bin_macosx64/wb_command', type=str, help="string indicating path to wb_command")
    args = vars(parser.parse_args())
    
    # Set up subject list
    wb = args['wb']
    proj_dir = args['proj_dir']
    hcp_pheno = [x for x in os.listdir(proj_dir) if '.py' not in x]
    args = vars(parser.parse_args())
    hcp_pheno = [x for x in os.listdir(proj_dir) if '.py' not in x]
    if args['first_id']:
        sub_list = hcp_pheno[int(args['first_id'])-1:int(args['last_id'])]
    elif args['startswith']:
        sub_list = [x for x in hcp_pheno if x.startswith(args['startswith'])]
    else:
        sub_list = hcp_pheno

    # Loop over subjects and parcellate cifti files
    parc_files = ['/data/data7/restsync/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR_use.dlabel.nii']#['/data/data1/atlases/Glasser_et_al_2016_HCP_MMP1.0_RVVG/HCP_PhaseTwo/Q1-Q6_RelatedValidation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii']
    for sub in sub_list:
        print(sub)
        sub_dir = os.path.join(proj_dir,sub)

        tasks = [x for x in os.listdir(os.path.join(sub_dir,'MNINonLinear/Results')) if 'DS' not in x and '.txt' not in x]
        if args['tasktype']:
            tasks = [x for x in tasks if args['tasktype'] in x]

        for task in tasks:
            print('\t' + task)
            dtseries_file = os.path.join(sub_dir,'MNINonLinear/Results',task, task + '_Atlas_hp2000_clean_nilearn.dtseries.nii')
            for parc_file in parc_files:
                if parc_file == '/data/data1/atlases/Glasser_et_al_2016_HCP_MMP1.0_RVVG/HCP_PhaseTwo/Q1-Q6_RelatedValidation210/MNINonLinear/fsaverage_LR32k/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii':
                    parc_name = 'glasser360v'
                elif parc_file == '/data/data7/restsync/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR_use.dlabel.nii':
                    parc_name = 'cole360'
                ptseries_file = dtseries_file.replace('.dtseries.nii', '_' + parc_name + '.ptseries.nii')
                corr_cmd = wb + ' -cifti-parcellate' + ' ' + dtseries_file + ' ' + parc_file + ' COLUMN ' + ptseries_file + ' -only-numeric' + ''
                if os.path.exists(dtseries_file):
                    os.system(corr_cmd)
                    if args['regress_motion'] == 1:
                        output = ptseries_file.replace("ptseries.nii", "csv")
                    elif args['regress_motion'] == 0:
                        output = ptseries_file.replace(".ptseries.nii", "GSR_only.csv")
                    conv_cmd = wb + ' -cifti-convert -to-text' + ' ' + ptseries_file + ' ' + output + ' -col-delim ' + ','
                    print('\t' + 'created: ' + output)
                    if os.path.exists(ptseries_file):

                    	# Delete intermediate files
                        os.system(conv_cmd)
                        os.remove(ptseries_file)



if __name__ == "__main__":
    main()