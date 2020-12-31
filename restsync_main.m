%%%%%%%%
% Author: David Gruskin
% Contact: dcg2153@cumc.columbia.edu
% Last updated: 12/2020
% Project: Brain connectivity at rest predicts individual differences in normative activity during movie watching 
% Description: This is the main script for the HCP 7T RSFC-ISC project. 

% Dependencies: 
% [calc_FD.m, make_friston_regressors.m] from https://github.com/DCAN-Labs/dcan_bold_processing
% [kfold_family.m, computeLambdaMax.m, find_lambda.m] from https://github.com/YaleMRRC/CPM/
% [fd_censoring.m, ridgeCPM_bagging.m, hcp_vis.m] from https://github.com/davidgruskin/restsync
% [rotate_parcellation.m, perm_sphere_p.m] from https://github.com/frantisekvasa/rotate_parcellation
%%%%%%%%

%% Define variables
% Set these variables before running through the code. If running from
% scratch, create_filtered_fd_files, run_setup, and load_data should be set
% to 1. Otherwise, configure as necessary.

% Set to 1 to calculate FC/ISC/correlation between FC and ISC
run_setup = 1;

% Load pre-generated FC,ISC, and FC-ISC outputs
if run_setup == 0
    load('/Users/davidgruskin/HCP_Restsync/data/imaging/restsync_results_231220.mat')
end

% Name of parcellation being used
parc_name = 'glasser360v';

% Set to 1 to generate filtered FD traces
create_filtered_fd_files = 0;

% Exclude scans if > this fraction of frames are censored
fd_perc = .5;

% Set some motion censoring values here
fd_threshold   = 0.2; % remove all frames with fd greater than this value
num_contiguous = 5; % if a run of uncensored frames is less than this value, censor them

% Set to 1 to load raw data from csvs
load_data = 1;

% Set to 1 to create scalars for visualization
create_scalars = 0;

% Set to 1 if on kraepelin/mountcastle
kraepelin = 1;

% Path to workbench
wb_path = '/Users/davidgruskin/Downloads/workbench/bin_macosx64/wb_command';


%% Add paths and load setup files

if kraepelin == 1
    addpath('/data/data7/restsync/')
    addpath('/data/data1/bin/')
    addpath('/home/davidgruskin/Desktop')
    addpath('/data/data1/bin/matlab/n42_analyses/')
    addpath('/data/data7/restsync/spin_dir/')
    proj_dir = '/data/data7/HCP_7T';
    load('/data/data7/restsync/final_code/netassignments_glasser.mat');
    
    % Flip the Cole/Anticevic network labels to account for R/L swap
    netassignments_new(1:180,1) = netassignments(181:360);
    netassignments_new(181:360,1) = netassignments(1:180);
else
    addpath('/Users/davidgruskin/Downloads')
    addpath('/Users/davidgruskin/HCP_Restsync')
    addpath('/Users/davidgruskin/Downloads/RainCloudPlots-master/')
    addpath('/Users/davidgruskin/Downloads/Robust_Statistical_Toolbox-master/')
    
    % Set up visualization
    load('/Users/davidgruskin/HCP_Restsync/miscellaneous/visualization_mats/netassignments_coleanticevic.mat')
   
    % Flip the Cole/Anticevic network labels to account for R/L swap
    netassignments_new(1:180,1) = netassignments(181:360);
    netassignments_new(181:360,1) = netassignments(1:180);

    save_stem        = '/Users/davidgruskin/HCP_Restsync/';
    wb_path          = '~/Downloads/workbench/bin_macosx64/wb_command';
    scalar_template  = '/Users/davidgruskin/HCP_Restsync/workbench_files/workbench_setup/glasser_360v.pscalar.nii';
    pconn_template   = '/Users/davidgruskin/HCP_Restsync/workbench_files/workbench_setup/archive/glasser.pconn.nii';
end

% Set up parcellation/RSN variables
num_node                   = 360;
num_edges                  = 64620;
network_names              = {'Visual1','Visual2','Somatomotor','Cingulo-Opercular','Dorsal-attention','Language','Frontoparietal','Auditory','Default','Posterior-Multimodal','Ventral-Multimodal','Orbito-Affective'};
coleanticevic_rgb_networks = [0, 0, 225;100, 0, 255; 0, 255, 255; 153, 0 ,153; 0, 255, 0; 0, 154, 154; 255, 255, 0; 249, 61, 251; 255, 0, 0; 177, 89, 40; 255, 156, 0; 65, 124, 0];
coleanticevic_rgb_nodes    = zeros(num_node,3);
for node = 1:num_node
    coleanticevic_rgb_nodes(node,:) = coleanticevic_rgb_networks(netassignments(node,1),:);
end

%% Data loading

% Prepare pheno tables
if kraepelin == 0
    load('/Users/davidgruskin/inputs/hcp_mztwin.mat')
    load('/Users/davidgruskin/inputs/hcp_dztwin.mat')
    
    pheno_table      = readtable('/Users/davidgruskin/inputs/hcp_7t_pheno.csv');
    pheno_restricted = readtable('/Users/davidgruskin/inputs/RESTRICTED_dgruskin_3_8_2019_15_47_21.csv');
    
else
    load('/data/data7/restsync/hcp_mztwin.mat')
    load('/data/data7/restsync/hcp_dztwin.mat')
    
    pheno_table      = readtable('/data/data7/restsync/hcp_7t_pheno.csv');
    pheno_restricted = readtable('/data/data7/restsync/RESTRICTED_dgruskin_3_8_2019_15_47_21.csv');
end

pheno_restricted = [pheno_restricted(:,1:5) array2table(MZTWIN) array2table(DZTWIN) pheno_restricted(:,6:end) pheno_table(:,2:end)];
IDs              = pheno_restricted.Subject;

% Set basic variables
num_subj  = 184;
num_scan  = 8;
num_net   = 12;
num_days  = 2;

% These subjects are missing some amount of rest/movie data, so they
% will be excluded from analysis
all_excludes = [8,66,124,126,130,135,179,183];
all_used               = 1:num_subj;
all_used(all_excludes) = [];

% Get demographics
age_mean = mean(pheno_restricted.Age_in_Yrs(all_used,1));
age_std  = std(pheno_restricted.Age_in_Yrs(all_used,1));

num_female = sum(grp2idx(pheno_table.Gender(all_used,1))==2);

race_percentages = countcats(categorical(pheno_restricted.Race(all_used,1)))/length(all_used);

ethnicity_percentages = countcats(categorical(pheno_restricted.Ethnicity(all_used,1)))/length(all_used);

% Set scan names for data loading
scan_names    = cell(8,1);
scan_names{1} = 'rfMRI_REST1_7T_PA';scan_names{2} = 'rfMRI_REST2_7T_AP';scan_names{3} = 'rfMRI_REST3_7T_PA';scan_names{4} = 'rfMRI_REST4_7T_AP';
scan_names{5} = 'tfMRI_MOVIE1_7T_AP';scan_names{6} = 'tfMRI_MOVIE2_7T_PA';scan_names{7} = 'tfMRI_MOVIE3_7T_PA';scan_names{8} = 'tfMRI_MOVIE4_7T_AP';

% Set family IDs
family_ids(:,1) = 1:num_subj;
family_ids(:,2) = grp2idx(pheno_restricted.Family_ID);

%% Create filtered fd files (this section adapted from https://github.com/DCAN-Labs/dcan_bold_processing)
if create_filtered_fd_files == 1
    head_ratio_cm = 5;
    TR            = 1;
    LP_freq_min   = 12;
    order         = 4;
    
    % filter design
    hr_min = LP_freq_min; % recasted to reuse code
    hr     = hr_min/60;
    fs     = 1/TR;
    fNy    = fs/2;
    fa     = abs(hr - floor((hr + fNy) / fs) * fs);
    
    % cutting frequency normalized between 0 and nyquist
    Wn = min(fa)/fNy;
    if ~isempty(order)
        b_filt = fir1(order, Wn, 'low');
        a_filt = 1;
    end
    num_f_apply = 0;
    
    % Read individual movement regressors files
    if create_filtered_fd_files == 1
        for subj = 1:num_subj
            path_mov_reg  = strcat('/data/data7/HCP_7T/',num2str(IDs(subj,1)),'/MNINonLinear/Results');
            pathstring    = [path_mov_reg filesep '*' filesep 'Movement_Regressors.txt'];
            path_contents = dir(pathstring);
            
            n = size(path_contents,1);
            for i = 1:n
                % Read motion numbers
                file_mov_reg = [path_contents(i).folder filesep path_contents(i).name]; %-- for MATLAB USE
                MR           = dlmread(file_mov_reg);
                MR_ld        = make_friston_regressors(MR);%% Using this function to only get the linear displacements
                MR_ld        = MR_ld(:,1:6);
                split_name   = split(file_mov_reg,'/');
                scan_id      = split_name{end-1};
                
                %'FiltFilt_all';
                MR_filt = filtfilt(b_filt,a_filt,MR_ld);
                for i = 1:num_f_apply-1
                    MR_filt = filtfilt(b_filt,a_filt,MR_filt);
                end
                
                
                hd_mm     = 10 * head_ratio_cm; % cm to mm conversion
                MR_backed = MR_filt;
                
                MR_backed(:, 4:end) = 180*MR_backed(:,4:end)/(pi*hd_mm);
                
                %% The last 6 movement regressors are derivatives of the first 6. Needed for task fMRI, but not used in the Fair Lab
                second_derivs = cat(1, [0 0 0 0 0 0], diff(MR_backed(:,1:6)));
                MR_backed     = [MR_backed second_derivs];
                MR_new        = MR_backed(:, [1 2 3 4 5 6]);
                MR_new(:,4:end) = MR_new(:,4:end)*pi*50/180; % Calculate length of arc in mm
                
                dX = diff(MR_new); % calculate derivatives
                FD = (sum(abs(dX),2)); % calculate FD
                
                filename_out = strcat(proj_dir,'/',num2str(IDs(subj,1)),'/MNINonLinear/Results/',scan_id,'/FD1.mat');
                save(filename_out,'FD')
            end
        end
        %% Get indices of censored volumes
        
        if create_filtered_fd_files == 1
            for scan_id = 1:num_scan
                for subj = 1:length(IDs)
                    filename = strcat(proj_dir,'/',num2str(IDs(subj,1)),'/MNINonLinear/Results/',scan_names{scan_id},'/FD1.mat');
                    if exist(filename,'file') == 2
                        FD = load(strcat(proj_dir,'/',num2str(IDs(subj,1)),'/MNINonLinear/Results/',scan_names{scan_id},'/FD1.mat'));
                        vols_to_censor = fd_censoring(FD.FD,fd_threshold,num_contiguous);
                        if ~isempty(vols_to_censor)
                            vols_mat = zeros(length(FD.FD)+1,size(vols_to_censor,1));
                            for row = 1:size(vols_to_censor,1)
                                vols_mat(vols_to_censor(row,1),row) = 1;
                            end
                            if create_filtered_fd_files == 1
                                writematrix(vols_mat,strcat(proj_dir,'/',num2str(IDs(subj,1)),'/MNINonLinear/Results/',scan_names{scan_id},'/motion_regressors1.csv'));
                                clear vols_mat FD
                            end
                        end
                    end
                end
            end
        end
    end
end

%% GSR/Spike Regression and Parcellation

% After running all of the code up until this point, a file called
% "motion_regressors1.csv" should have been created for all scans in which
% there are frames that exceed the specified FD threshold. To censor these
% frames and remove the global signal (and its derivative) from the data,
% run "restsync_denoise.py" followed by "restsync_parcellate.py," which
% will generate csvs with denoised, parcellated data to be loaded in the
% next section. Instructions for running the python scripts can be found in
% the header of each script.

%% Load imaging data

% Initialize matrices
sub_data           = cell(num_scan,1); % this variable will contain the parcel-wise raw timecourses
perc_mat           = zeros(num_subj,num_scan); % this variable will contain the % of frames censored per individual/scan
mean_FD_uncensored = NaN(num_subj,num_scan); % this variable will contain the mean FD of the uncensored FD trace per individual/scan
mean_FD_unfiltered = NaN(num_subj,num_scan); % ' ' but for the unfiltered FD trace
mean_FD_censored   = NaN(num_subj,num_scan); % ' ' but for the censored FD trace
full_fd_mat        = cell(num_scan,1);

if load_data == 1
    for scan_id = 1:num_scan
        scan_id
        % Regress motion only for rest scans
        if scan_id < 5
            regress_motion = 1;
            exemplar_subdata = csvread(strcat(proj_dir,'/','100610','/MNINonLinear/Results/',scan_names{scan_id},'/',scan_names{scan_id},sprintf('_Atlas_hp2000_clean_nilearn_%s.csv',parc_name)));
        else
            regress_motion = 0;
            exemplar_subdata = csvread(strcat(proj_dir,'/','100610','/MNINonLinear/Results/',scan_names{scan_id},'/',scan_names{scan_id},sprintf('_Atlas_hp2000_clean_nilearn_%sGSR_only.csv',parc_name)));
        end
        for subj = 1:length(IDs)
            if regress_motion == 1
                filename = strcat(proj_dir,'/',num2str(IDs(subj,1)),'/MNINonLinear/Results/',scan_names{scan_id},'/',scan_names{scan_id},sprintf('_Atlas_hp2000_clean_nilearn_%s.csv',parc_name));
            elseif regress_motion == 0
                filename = strcat(proj_dir,'/',num2str(IDs(subj,1)),'/MNINonLinear/Results/',scan_names{scan_id},'/',scan_names{scan_id},sprintf('_Atlas_hp2000_clean_nilearn_%sGSR_only.csv',parc_name));
            end
            if exist(filename,'file') == 2
                try
                    csvdata = csvread(filename);
                    if size(csvdata,2) == size(exemplar_subdata,2)
                        sub_data{scan_id}(:,:,subj)      = csvdata;
                        FD                               = load(strcat(proj_dir,'/',num2str(IDs(subj,1)),'/MNINonLinear/Results/',scan_names{scan_id},'/FD1.mat'));
                        full_fd_mat{scan_id}(subj,:)     = FD.FD;
                        mean_FD_uncensored(subj,scan_id) = nanmean(FD.FD);
                        vols_to_censor                   = fd_censoring(FD.FD,fd_threshold,num_contiguous);
                        perc_mat(subj,scan_id)           = length(vols_to_censor)/(length(FD.FD));
                        hcp_motion_numbers               = calcFD(dlmread(strcat(proj_dir,'/',num2str(IDs(subj,1)),'/MNINonLinear/Results/',scan_names{scan_id},'/Movement_Regressors_dt.txt')));
                        mean_FD_unfiltered(subj,scan_id) = mean(hcp_motion_numbers);
                        
                        if ~isempty(vols_to_censor) && regress_motion == 1
                            vols_mat = zeros(length(FD.FD)+1,size(vols_to_censor,1));
                            for row = 1:size(vols_to_censor,1)
                                vols_mat(vols_to_censor(row,1),row) = 1;
                            end
                            
                            sub_data{scan_id}(:,vols_to_censor,subj) = NaN;
                            FD.FD(vols_to_censor,:) = NaN;
                            mean_FD_censored(subj,scan_id) = nanmean(FD.FD);
                        end
                        
                    else
                        sub_data{scan_id}(:,:,subj) = NaN;
                    end
                    
                catch
                    sub_data{scan_id}(:,:,subj) = NaN;
                end
                
            else
                sub_data{scan_id}(:,:,subj) = NaN;
            end
        end
    end
end

%% Clean up movie timecourses

% Each movie clip is preceded by 20 seconds of rest. Here, we remove those
% 20 seconds as well as the first 20 seconds of each clip to account for
% potential onset transients.

% Initialize matrices
rest_trs   = cell(4,1);
movie_cut  = cell(4,1);
movie_mask = cell(4,1);

% Rest blocks- Run 1
starts = [0,264.0833,505.75,713.7917,797.5833,901];
starts = floor(starts+1);
ends   = [19.9583,284.0417,525.7083,733.75,817.5417,920.9583];
ends   = floor(ends+1);

for i = 1:length(starts)
    if i == 1
        rest_trs{1,1} = [rest_trs{1,1} starts(i):ends(i)];
    else
        rest_trs{1,1} = [rest_trs{1,1} starts(i)+5:ends(i)];
    end
end

% Rest blocks- Run 2
starts = [0,246.75,525.375,794.625,898];
starts = floor(starts+1);
ends   = [19.9583,266.7083,545.3333,814.5417,917.9583];
ends   = floor(ends+1);

for i = 1:length(starts)
    if i == 1
        rest_trs{2,1} = [rest_trs{2,1} starts(i):ends(i)];
    else
        rest_trs{2,1} = [rest_trs{2,1} starts(i)+5:ends(i)];
    end
end

% Rest blocks- Run 3
starts = [0,200.5833,405.125,629.25,791.7917,895];
starts = floor(starts+1);
ends   = [19.9583,220.5417,425.0833,649.2083,811.5417,914.9583];
ends   = floor(ends+1);

for i = 1:length(starts)
    if i == 1
        rest_trs{3,1} = [rest_trs{3,1} starts(i):ends(i)];
    else
        rest_trs{3,1} = [rest_trs{3,1} starts(i)+5:ends(i)];
    end
end

% Rest blocks- Run 4
starts = [0,252.3333,502.2083,777.4167,881];
starts = floor(starts+1);
ends   = [19.9583,272.2917,522.1667,797.5417,900.9583];
ends   = floor(ends+1);

for i = 1:length(starts)
    if i == 1
        rest_trs{4,1} = [rest_trs{4,1} starts(i):ends(i)];
    else
        rest_trs{4,1} = [rest_trs{4,1} starts(i)+5:ends(i)];
    end
end

% Movie blocks- Run 1
starts = [20,284.0833,525.75,733.7917,817.5833];
ends   = [264.0417,505.7083,713.75,797.5417,900.9583];
starts = ceil(starts+1);
ends   = floor(ends);

for i = 1:length(starts)
    movie_cut{1,1} = [movie_cut{1,1} starts(i):starts(i)+19];
    movie_mask{1,1}{i} = [starts(i)+19 ends(i)];
end

% Movie blocks- Run 2
starts = [20,266.75,545.375,814.5833];
ends   = [246.7083,525.3333,794.5833,897.9583];
ends   = floor(ends);
starts = ceil(starts+1);

for i = 1:length(starts)
    movie_cut{2,1} = [movie_cut{2,1} starts(i):starts(i)+19];
    movie_mask{2,1}{i} = [starts(i)+19 ends(i)];
end

% Movie blocks- Run 3
starts = [20,220.5833,425.125,649.25,811.5833];
ends   = [200.5417,405.0833,629.2083,791.75,894.9583];
ends   = floor(ends);
starts = ceil(starts+1);

for i = 1:length(starts)
    movie_cut{3,1} = [movie_cut{3,1} starts(i):starts(i)+19];
    movie_mask{3,1}{i} = [starts(i)+19 ends(i)];
end

% Movie blocks- Run 4
starts = [20,272.3333,522.2083,797.5833];
ends   = [252.2917,502.1667,777.375,880.9583];
ends   = floor(ends);
starts = ceil(starts+1);

for i = 1:length(starts)
    movie_cut{4,1} = [movie_cut{4,1} starts(i):starts(i)+19];
    movie_mask{4,1}{i} = [starts(i)+19 ends(i)];
end

% Get indices of all frames that need to be removed
all_cut      = cell(4,1);
all_cut{1,1} = unique(cat(2,rest_trs{1,1},movie_cut{1,1}));
all_cut{2,1} = unique(cat(2,rest_trs{2,1},movie_cut{2,1}));
all_cut{3,1} = unique(cat(2,rest_trs{3,1},movie_cut{3,1}));
all_cut{4,1} = unique(cat(2,rest_trs{4,1},movie_cut{4,1}));

% Extract movies by clips (not used here but might be helpful to have)
movie_data = cell(4,1);
for scan_id = 5:8
    for clip_id = 1:size(movie_mask{scan_id-4,1},2)
        movie_data{scan_id-4,clip_id} = sub_data{scan_id,1}(:,(movie_mask{scan_id-4,1}{clip_id}(1,1):movie_mask{scan_id-4,1}{clip_id}(1,2)),:);
    end
end

% Remove rest/onset frames
for scan_id = 5:8
    sub_data{scan_id,1}(:,all_cut{scan_id-4,1},:) = [];
end

% Get mean FD in concatenated rest and movie scans
full_fd_cat{1,1} = cat(2,full_fd_mat{1},full_fd_mat{2});
full_fd_cat{2,1} = cat(2,full_fd_mat{3},full_fd_mat{4});
full_fd_cat{3,1} = cat(2,full_fd_mat{5},full_fd_mat{6});
full_fd_cat{4,1} = cat(2,full_fd_mat{6},full_fd_mat{8});

mean_FD_uncensored_cat = NaN(num_subj,4);
for scan_id = 1:4
    mean_FD_uncensored_cat(:,scan_id) = mean(full_fd_cat{scan_id,1},2);
end

% Define function for z-scoring with NaNs
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

% Z-score all timecourses
for scan_id = 1:8
    for subj = 1:length(IDs)
        sub_data{scan_id,1}(:,:,subj) = zscor_xnan(squeeze(sub_data{scan_id,1}(:,:,subj))')';
    end
end

% Concatenate rest and movie scans from the same day
% Order is: Rest Day 1, Rest Day 2, Movie Day 1, Movie Day 2
sub_data_cat{1,1}  = cat(2,sub_data{1},sub_data{2});
sub_data_cat{2,1}  = cat(2,sub_data{3},sub_data{4});
sub_data_cat{3,1}  = cat(2,sub_data{5},sub_data{6});
sub_data_cat{4,1}  = cat(2,sub_data{7},sub_data{8});

% Get percentage of frames censored
perc_mat = zeros(num_subj,4);
for scan_id = 1:4
    for subj = 1:num_subj
        perc_mat(subj,scan_id) = length(find(isnan(sub_data_cat{scan_id,1}(1,:,subj))))/length(sub_data_cat{scan_id,1});
    end
end

%% Remove subjects who have incomplete data

for scan_num = 1:4
    sub_data_cat{scan_num,1}(:,:,all_excludes) = NaN;
end

%% Calculate FC
hcp_fc   = cell(4,1);
num_node = size(sub_data_cat{1,1},1);

for scan_id = 1:4
    fc_out = zeros(num_node,num_node,num_subj);
    parfor subj = 1:num_subj
        fc_out(:,:,subj) = corr(squeeze(sub_data_cat{scan_id,1}(:,:,subj))',squeeze(sub_data_cat{scan_id,1}(:,:,subj))','rows','complete','type','pearson');  %#ok<*PFBNS>
    end
    hcp_fc{scan_id,1} = fc_out;
end

parfor scan_id = 1:4
    for subj = 1:num_subj
        for node = 1:num_node
            hcp_fc{scan_id,1}(node,node,subj) = NaN;
        end
    end
end

%% Calculate LOO ISC

hcp_isc = zeros(num_subj,num_node,2);
for scan_id = 1:2
    sliced = sub_data_cat{scan_id+2,1};
    for subj = 1:num_subj
        tmp           = sliced;
        targ          = squeeze(sliced(:,:,subj))';
        tmp(:,:,subj) = [];
        tmp_mean      = nanmean(tmp,3)';
        for node = 1:size(sliced,1)
            hcp_isc(subj,node,scan_id) = corr(targ(:,node),tmp_mean(:,node),'rows','complete','type','pearson');
        end
    end
end

%% Calculate relationships between FC and ISC

% Initialize matrices
X        = zeros(num_node);
edges    = (num_node*(num_node-1))/2;
reshaped = zeros(edges, num_subj);

output_matrix_r = zeros(num_node,edges,2);
output_matrix_p = zeros(num_node,edges,2);
square_matrix_r = NaN(num_node,num_node,num_node,4);

% Get square indices
[i,j] = find(tril(ones(num_node), -1));

counter = 1;
for isc_scan = 1:2
    for fc_scan = 1:2
        
        % Reshape square matrix to 1D vector
        reshaped = zeros(edges, num_subj);
        for subj = 1:num_subj
            X                = squeeze(hcp_fc{fc_scan,1}(:,:,subj));
            reshaped(:,subj) = X(logical(tril(ones(size(X)),-1)));
        end
        
        % Calculate correlation between RSFC and ISC across participants
        fd_cov = mean_FD_uncensored_cat(:,[fc_scan,isc_scan+2]);
        parfor node = 1:num_node
            [output_matrix_r(node,:,fc_scan),output_matrix_p(node,:,fc_scan)] = partialcorr(squeeze(hcp_isc(:,node,isc_scan)),reshaped',fd_cov,'rows','complete','type','spearman') ;
        end
        
        % Reshape 1D vector to square matrix
        for node = 1:num_node
            for index = 1:length(i)
                square_matrix_r(node,i(index),j(index),counter) = output_matrix_r(node,index,fc_scan);
                square_matrix_r(node,j(index),i(index),counter) = output_matrix_r(node,index,fc_scan);
            end
        end
        
        counter = counter + 1;
    end
end

%% Get gISC test-retest reliability

% Get mean ISC across participants
mean_isc = conv_z2r(squeeze(nanmean(conv_r2z(hcp_isc),2)));

% Get mean ISC across parcels
parcel_mean_isc = conv_z2r(squeeze(nanmean(conv_r2z(hcp_isc))));

% How consistent are ISC maps across both days?
parcel_mean_trt = corr(parcel_mean_isc,'type','spearman');

if kraepelin == 1  
    perm_id            = rotate_parcellation('/data/data7/spin_dir/left_out.txt','/data/data7/spin_dir/right_out.txt',10000);
    parcel_mean_trt_p  = perm_sphere_p(parcel_mean_isc(:,1),parcel_mean_isc(:,2),perm_id,'spearman');
end

% Make scalars
if create_scalars == 1
    group_name = 'hcp1';
    file_name  = strcat(parc_name,'_',group_name,'_parcel_meanISC');
    write_vec  = parcel_mean_isc(:,1);
    hcp_vis(write_vec,scalar_template,file_name);
    
    group_name = 'hcp2';
    file_name  = strcat(parc_name,'_',group_name,'_parcel_meanISC');
    write_vec  = parcel_mean_isc(:,2);
    hcp_vis(write_vec,scalar_template,file_name);
end

% Residualize ranks with motion (FD) traces
mdl1  = fitlm(cat(2,tiedrank(mean_FD_uncensored_cat(:,3)),tiedrank(mean_FD_uncensored_cat(:,4))),tiedrank(mean_isc(:,1)));
mdl2  = fitlm(cat(2,tiedrank(mean_FD_uncensored_cat(:,3)),tiedrank(mean_FD_uncensored_cat(:,4))),tiedrank(mean_isc(:,2)));
[r,p] = partialcorr(mean_isc(:,1),mean_isc(:,2),mean_FD_uncensored_cat(:,3:4),'rows','complete','type','spearman');

% Create scatter of day1/day2 gISC
figure
scatter(mdl1.Residuals.Raw,mdl2.Residuals.Raw,'filled'); lsline

xlim([-100 150]); ylim([-100 150])
title('ISC Test-Retest Reliability')
xlabel('Day 1 gISC'); ylabel('Day 2 gISC')
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')
rtxt = sprintf('Spearman Rho = %s \np-value= %s',num2str(r), num2str(p));
t1 = TextLocation(rtxt,'Location','Best'); t1.FontSize = 14;

% Create scatter of day1/day2 gISC (raw)
figure
scatter(mean_isc(:,1),mean_isc(:,2),'filled'); lsline

title('ISC Test-Retest Reliability')
xlabel('Day 1 gISC'); ylabel('Day 2 gISC')
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')
rtxt = sprintf('Spearman Rho = %s \np-value= %s',num2str(r), num2str(p));
t1 = TextLocation(rtxt,'Location','Best'); t1.FontSize = 14;

%% Ridge CPM

% First, run 100 iterations of 10-fold CV to train model on Day 1 data and
% test using day 2 data. This will take a while to run (~hours).
[mdl_bagged, ridgeCPMoutput] = ridgeCPM_bagging(hcp_fc, hcp_isc, all_used, sub_data_cat, mean_FD_uncensored_cat, family_ids, 100, 0);

% Next, run 10,000 iterations of the same model with shuffled mean ISC
% scores to generate null distribution for CV comparison. This will take a
% while to run (~days).
[mdl_bagged_perm, ridgeCPMoutput_perm] = ridgeCPM_bagging(hcp_fc, hcp_isc, all_used, sub_data_cat, mean_FD_uncensored_cat, family_ids, 10000, 1);

%% Ridge CPM visualizations

% Get permutation p_value
perm_p = (length(find(ridgeCPMoutput_perm.spearman_out > median(ridgeCPMoutput.spearman_out)))+1)/10001;

% Get order of parcels by RSN
counter         = 1;
indices_reshape = zeros(num_node,1);

for net = 1:num_net
    indices = find(netassignments_new(:,1)==net);
    indices_reshape((counter:counter+length(indices)-1),1) = indices;
    
    counter = counter + length(indices);
end

% Find average ridge coefficient for each FC edge
all_betas = zeros(num_edges,1);
for x = 1:length(mdl_bagged.edge_idx)
    all_betas(mdl_bagged.edge_idx(x),1) = mdl_bagged.betas(x,1);
end

% Reshape average ridge coefficients to square
square_matrix_cpm = zeros(num_node,num_node);
[i,j] = find(tril(ones(num_node), -1));
for index = 1:length(i)
    square_matrix_cpm(i(index),j(index)) = all_betas(index);
    square_matrix_cpm(j(index),i(index)) = all_betas(index);
end

% Separate parcel IDs by RSN affiliation
net_edges = cell(num_net,1);
for net = 1:num_net
    net_edges{net,1} = find(netassignments_new(:,1) == net);
end

% Get within network edge IDs
net_coords = cell(num_net,1);
for net = 1:num_net
    net_coords{net,1} = combnk(net_edges{net,1},2);
end

% Get between network edge IDs
net_coords_between = cell(num_net,num_net);
for net1 = 1:num_net
    for net2 = 1:num_net
        row_counter = 1;
        for row1 = 1:(size(net_edges{net1,1},1))
            for row2 = 1:(size(net_edges{net2,1},1))
                net_coords_between{net1,net2}(row_counter,1) = net_edges{net1,1}(row1,1);
                net_coords_between{net1,net2}(row_counter,2) = net_edges{net2,1}(row2,1);
                
                row_counter = row_counter +1;
            end
        end
    end
end

% Combine within + between network edge IDs
net_coords_all = net_coords_between;
for net = 1:num_net
    net_coords_all{net,net} = net_coords{net,1};
end

% Create heatmaps to show ridge coefficient represenation by RSN
cont_pos = zeros(num_node);
cont_neg = zeros(num_node);

for node1 = 1:num_node
    for node2 = 1:num_node
        if square_matrix_cpm(node1,node2)>0
            cont_pos(node1,node2) = square_matrix_cpm(node1,node2);
        elseif square_matrix_cpm(node1,node2)<0
            cont_neg(node1,node2) = square_matrix_cpm(node1,node2);
        else
        end
    end
end

% Negative coefficients
neg_net_out = zeros(num_net,num_net);
for net1 = 1:num_net
    for net2 = 1:num_net
        if neg_net_out(net1,net2) ~=0
            neg_net_out(net2,net1) = 0;
        else
            tmp = 0;
            counter = 0;
            for row = 1:size(net_coords_all{net1,net2},1)
                tmp = tmp + cont_neg(net_coords_all{net1,net2}(row,1),net_coords_all{net1,net2}(row,2));
                if cont_neg(net_coords_all{net1,net2}(row,1),net_coords_all{net1,net2}(row,2)) ~=0
                    counter = counter + 1;
                end
            end
            neg_net_out(net2,net1) = tmp/counter;
        end
    end
end
neg_net_out(isnan(neg_net_out)) = 0; neg_net_out(isinf(neg_net_out)) = 0;
neg_net_out(:,13) = 0;
neg_net_out(13,:) = 0;

ax(1) = gca;
pcolor(neg_net_out); pbaspect([1 1 1])
caxis([min(min(neg_net_out)),0]); colormin = [50 50 200]; colormax = [255 255 255]; colorneutral = [152 152 227];
n    = 1000;

Rlow  = linspace(colormin(1)/255,colorneutral(1)/255,n); Rhigh = linspace(colorneutral(1)/255,colormax(1)/255,n);
Blow  = linspace(colormin(3)/255,colorneutral(3)/255,n); Bhigh = linspace(colorneutral(3)/255,colormax(3)/255,n);
Glow  = linspace(colormin(2)/255,colorneutral(2)/255,n); Ghigh = linspace(colorneutral(2)/255,colormax(2)/255,n);

colormap(ax(1), [cat(1,Rlow(:),Rhigh(:)), cat(1,Glow(:),Ghigh(:)), cat(1,Blow(:),Bhigh(:))] );  %// create colormap
colorbar;
set(gca, 'YDir','reverse')

axis image

%xticks(1:num_net); yticks(1:num_net)
%xticklabels(network_names); yticklabels(network_names)
%xtickangle(45); set(gca,'TickLength',[0 0])

% Positive coefficients
pos_net_out = zeros(num_net,num_net);
for net1 = 1:num_net
    for net2 = 1:num_net
        if pos_net_out(net1,net2) ~=0
            pos_net_out(net2,net1) = 0;
        else
            tmp = 0;
            counter = 0;
            for row = 1:size(net_coords_all{net1,net2},1)
                tmp = tmp + cont_pos(net_coords_all{net1,net2}(row,1),net_coords_all{net1,net2}(row,2));
                if cont_pos(net_coords_all{net1,net2}(row,1),net_coords_all{net1,net2}(row,2)) ~=0
                    counter = counter + 1;
                end
            end
            pos_net_out(net2,net1) = tmp/counter;
        end
    end
end
pos_net_out(isnan(pos_net_out)) = 0; pos_net_out(isinf(pos_net_out)) = 0;

pos_net_out(:,13) = 0;
pos_net_out(13,:) = 0;
ax(2) = gca; s=pcolor(pos_net_out); pbaspect([1 1 1]);
caxis([0,max(max(pos_net_out))]); colormin = [255 255 255]; colormax = [200 50 50]; colorneutral = [227 152 152];

Rlow  = linspace(colormin(1)/255,colorneutral(1)/255,n); Rhigh = linspace(colorneutral(1)/255,colormax(1)/255,n);
Blow  = linspace(colormin(3)/255,colorneutral(3)/255,n); Bhigh = linspace(colorneutral(3)/255,colormax(3)/255,n);
Glow  = linspace(colormin(2)/255,colorneutral(2)/255,n); Ghigh = linspace(colorneutral(2)/255,colormax(2)/255,n);

colormap(ax(2),[cat(1,Rlow(:),Rhigh(:)), cat(1,Glow(:),Ghigh(:)), cat(1,Blow(:),Bhigh(:))] );  %// create colormap
c = colorbar;
xticks(1:num_net); yticks(1:num_net)
s.LineWidth = 1;
set(gca, 'YDir','reverse')
axis image

%xticks(1:num_net); yticks(1:num_net)
%xticklabels(network_names); yticklabels(network_names)
%xtickangle(45); set(gca,'TickLength',[0 0])


% Make scatter for day 1 -> day 2 prediction (ranks)
[r_rank, p_rank] = partialcorr(ridgeCPMoutput.predicted_ISC',ridgeCPMoutput.observed_ISC(:,1),ridgeCPMoutput.observed_ISC(:,2:end),'rows','complete','type','spearman');

mdl1 = fitlm(tiedrank(ridgeCPMoutput.observed_ISC(:,2:end)),tiedrank(ridgeCPMoutput.predicted_ISC'));
mdl2 = fitlm(tiedrank(ridgeCPMoutput.observed_ISC(:,2:end)),tiedrank(ridgeCPMoutput.observed_ISC(:,1)));

figure
scatter(mdl1.Residuals.Raw,mdl2.Residuals.Raw,'filled'); lsline

title('Ridge CPM Performance: Actual gISC vs. Predicted gISC')
xlabel('Predicted gISC'); ylabel('Actual gISC'); xlim([-100 150])
rtxt = strcat(['Spearman Rho = ' num2str(r_rank) ' p-value = ' num2str(p_rank)]);
t1   = TextLocation(rtxt,'Location','Best'); t1.FontSize = 14;
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')

% Make scatter for day 1 -> day 2 prediction (OG values)
figure
scatter(conv_z2r(ridgeCPMoutput.predicted_ISC'),conv_z2r(ridgeCPMoutput.observed_ISC(:,1)),'filled'); lsline

title('Ridge CPM Performance: Actual gISC vs. Predicted gISC')
xlabel('Predicted gISC'); ylabel('Actual gISC');
rtxt = strcat(['Spearman Rho = ' num2str(r_rank) ' p-value = ' num2str(p_rank)]);
t1   = TextLocation(rtxt,'Location','Best'); t1.FontSize = 14;
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')


% Compare null distribution from permutation tests with actual results
X          = ridgeCPMoutput_perm.spearman_out;
color      = [.5 .5 .5];
[f, Xi, ~] = ksdensity(X, 'bandwidth', []);
h{1}       = area(Xi, f); hold on

set(h{1}, 'FaceColor', color);
set(h{1}, 'EdgeColor', [0 0 0]);
set(h{1}, 'LineWidth', 2);
set(h{1}, 'FaceAlpha', 1);
xline(median(ridgeCPMoutput.spearman_out),'LineWidth',2,'Color',[1,0,0])


yl = get(gca, 'YLim');
set(gca, 'YLim', [-yl(2)*1 yl(2)]);

% Width of boxplot
wdth = yl(2) * 0.25;

% Jitter for raindrops
jit = (rand(size(X)) - 0.5) * wdth;

% Info for making boxplot
quartiles   = quantile(X, [0.25 0.75 0.5]);
iqr         = quartiles(2) - quartiles(1);
Xs          = sort(X);
whiskers(1) = min(Xs(Xs > (quartiles(1) - (1.5 * iqr))));
whiskers(2) = max(Xs(Xs < (quartiles(2) + (1.5 * iqr))));
Y           = [quartiles whiskers];

% Raindrops
drops_pos = jit - yl(2) / 2;

h{2}                 = scatter(X, drops_pos);
h{2}.SizeData        = 10;
h{2}.MarkerFaceColor = color;
h{2}.MarkerEdgeColor = 'none';

% Mean line
box_pos = [Y(1) -yl(2)/2-(wdth * 0.5) Y(2)-Y(1) wdth];

h{4} = line([Y(3) Y(3)], [-yl(2)/2-(wdth * 0.5) -yl(2) / 2 + (wdth * 0.5)], 'col', [0 0 0], 'LineWidth', 2);

h{3} = rectangle('Position', box_pos);
set(h{3}, 'EdgeColor', [0 0 0])
set(h{3}, 'LineWidth', 2);

% Whiskers
h{5} = line([Y(2) Y(5)], [-yl(2)/2 -yl(2)/2], 'col', [0 0 0], 'LineWidth', 2);
h{6} = line([Y(1) Y(4)], [-yl(2)/2 -yl(2)/2], 'col', [0 0 0], 'LineWidth', 2);

X     = ridgeCPMoutput.spearman_out;
color = [1 0 0];

% Jitter for raindrops
jit = (rand(size(X)) - 0.5) * wdth;

% Info for making boxplot
quartiles   = quantile(X, [0.25 0.75 0.5]);
iqr         = quartiles(2) - quartiles(1);
Xs          = sort(X);
whiskers(1) = min(Xs(Xs > (quartiles(1) - (1.5 * iqr))));
whiskers(2) = max(Xs(Xs < (quartiles(2) + (1.5 * iqr))));
Y           = [quartiles whiskers];

drops_pos = jit - yl(2) / 2;

h{2} = scatter(X, drops_pos);
h{2}.SizeData = 10;
h{2}.MarkerFaceColor = color;
h{2}.MarkerEdgeColor = 'none';
box_pos = [Y(1) -yl(2)/2-(wdth * 0.5) Y(2)-Y(1) wdth];
ylim([-5 5])

rtxt = strcat(['Median rho value = ' num2str(median(ridgeCPMoutput.spearman_out)) ' permutation p = ' num2str(perm_p)]);
t1 = TextLocation(rtxt,'Location','Best'); t1.FontSize = 14;
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')

%% Create outbound and inbound matrices and pscalars

% Initialize matrices
outbound_mat    = zeros(num_node,num_node,2);
inbound_mat     = zeros(num_node,num_node,2);

% Fill matrices
for day_id = 1:num_days
    for node1 = 1:num_node
        for node2 = 1:num_node
            outbound_mat(node1,node2,day_id) = squeeze(square_matrix_r(node2,node1,node2,day_id));
            inbound_mat(node1,node2,day_id)  = squeeze(square_matrix_r(node1,node1,node2,day_id));
        end
    end
end

% Make scalars
if create_scalars == 1
    for day_id = 1:2
        group_name = strcat('hcp',string(day_id));
        for node = 1:num_node
            
            % for outbound
            file_name = strcat(parc_name,'_',group_name,'_outbound_',string(node));
            write_vec = squeeze(outbound_mat(node,:,day_id));
            hcp_vis(write_vec,scalar_template,file_name);
            
            % for inbound
            file_name = strcat(parc_name,'_',group_name,'_inbound_',string(node));
            write_vec = squeeze(inbound_mat(node,:,day_id));
            %hcp_vis(write_vec,scalar_template,file_name);
        end
    end
end

% Make outbound/inbound pconns
if create_scalars == 1
    for day_id = 1:2
        inbound_mat_pconn  = squeeze(inbound_mat(:,:,day_id));
        outbound_mat_pconn = squeeze(outbound_mat(:,:,day_id));
        
        % for inbound
        cii_dscalar       = ciftiopen(pconn_template, wb_path);
        cii_dscalar.cdata = squeeze(inbound_mat_pconn(:,:));
        
        group_name = strcat('hcp',string(day_id));
        file_name  = strcat(parc_name,'_',group_name,'_inbound','.pconn.nii');
        save_path  = sprintf('%s%s',save_stem,file_name);
        ciftisave(cii_dscalar, save_path, wb_path)
        
        % for outbound
        cii_dscalar       = ciftiopen(pconn_template, wb_path);
        cii_dscalar.cdata = squeeze(outbound_mat_pconn(:,:));
        
        group_name = strcat('hcp',string(day_id));
        file_name  = strcat(parc_name,'_',group_name,'_outbound','.pconn.nii');
        save_path  = sprintf('%s%s',save_stem,file_name);
        ciftisave(cii_dscalar, save_path, wb_path)
    end
end

% Calculate outbound/inbound means
mean_outbound_absvalue = squeeze(nanmean(abs(conv_r2z(outbound_mat)),2)); mean_inbound_absvalue = squeeze(nanmean(abs(conv_r2z(inbound_mat)),2));
mean_outbound_noabs    = squeeze(nanmean(conv_r2z(outbound_mat),2)); mean_inbound_noabs = squeeze(nanmean(conv_r2z(inbound_mat),2));

% Make mean pscalars
if create_scalars == 1
    for day_id = 1:num_days
        group_name = strcat('hcp',string(day_id));
        
        % For outbound
        file_name = strcat(parc_name,'_',group_name,'_meanoutbound_absvalue');
        write_vec = conv_z2r(mean_outbound_absvalue(:,day_id));
        hcp_vis(write_vec,scalar_template,file_name);
        
        file_name = strcat(parc_name,'_',group_name,'_meanoutbound_noabs');
        write_vec = conv_z2r(mean_outbound_noabs(:,day_id));
        hcp_vis(write_vec,scalar_template,file_name);
        
        % For inbound
        file_name = strcat(parc_name,'_',group_name,'_meaninbound_absvalue');
        write_vec = conv_z2r(mean_inbound_absvalue(:,day_id));
        hcp_vis(write_vec,scalar_template,file_name);
        
        file_name = strcat(parc_name,'_',group_name,'_meaninbound_noabs');
        write_vec = conv_z2r(mean_inbound_noabs(:,day_id));
        hcp_vis(write_vec,scalar_template,file_name);
     
    end
end

%% T-tests and RM-ANOVAs for mean inbound/outbound maps
inbound_table      = cat(2,conv_z2r(squeeze(nanmean(conv_r2z(inbound_mat),2))),netassignments_new);
inbound_table      = array2table(inbound_table,'VariableNames',{'day1','day2','nets'});
inbound_table.nets = categorical(inbound_table.nets);

% First, for inbound

% Create boxplots for visualization
inbound_bpdata = cell(1,1);
inbound_bpdata_mat = zeros(78,num_net*num_days);
counter = 1;
for net = 1:num_net
    inbound_bpdata{net} = inbound_table{grp2idx(inbound_table.nets)==net,1:2};
    inbound_bpdata{net}(size(inbound_bpdata{net},1)+1:78,:) = NaN;
    for col= 1:2
        inbound_bpdata_mat(:,counter) = inbound_bpdata{net}(:,col) ;
        counter = counter+1;
    end
end
glasser_rgb_networks_flip = flip(coleanticevic_rgb_networks);

% Calculate inbound RSN averages
mean_inbound_RSN = conv_z2r(nanmean(conv_r2z(inbound_bpdata_mat))');

% Visualize with boxplots
figure; hold on
% Overlay individual points on boxplot
x = repmat(1:num_net*num_days,length(inbound_bpdata_mat),1);
for net = 1:num_net*num_days
    scatter(x(:,net),inbound_bpdata_mat(:,net),'filled','MarkerFaceAlpha',0.25','jitter','on','jitterAmount',0.15,'MarkerFaceColor',coleanticevic_rgb_networks(ceil(net/2),:)/255);
end

boxplot(inbound_bpdata_mat,'Symbol','w+');
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');

% Set colors according to Cole/Anticevic definitions
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),glasser_rgb_networks_flip(ceil(j/2),:)/255,'FaceAlpha',.5);
end

xticks(1.5:2:num_net*num_days); xticklabels(network_names); xtickangle(45); xlabel('Network')
ylabel('Average Relationship Between RSFC and ISC (Spearman Rho)'); title('Average of Inbound RSFC-ISC Maps for Each Parcel, by Resting State Network')
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')
ylim([-.1 .1])

% Now the same thing for outbound
outbound_table      = cat(2,conv_z2r(squeeze(nanmean(conv_r2z(outbound_mat),2))),netassignments_new);
outbound_table      = array2table(outbound_table,'VariableNames',{'day1','day2','nets'});
outbound_table.nets = categorical(outbound_table.nets);

outbound_bpdata = cell(1,1);
outbound_bpdata_mat = zeros(78,num_net*num_days);
counter = 1;
for net = 1:num_net 
    outbound_bpdata{net} = outbound_table{grp2idx(outbound_table.nets)==net,1:2};
    outbound_bpdata{net}(size(outbound_bpdata{net},1)+1:78,:) = NaN;
    for col = 1:2
        outbound_bpdata_mat(:,counter) = outbound_bpdata{net}(:,col) ;
        counter = counter+1;
    end
end
glasser_rgb_networks_flip = flip(coleanticevic_rgb_networks);

% Calculate inbound RSN averages
mean_outbound_RSN = conv_z2r(nanmean(conv_r2z(outbound_bpdata_mat))');

% Visualize with boxplots
figure; hold on
x = repmat(1:num_net*num_days,length(outbound_bpdata_mat),1);
for net = 1:num_net*num_days
    scatter(x(:,net),outbound_bpdata_mat(:,net),'filled','MarkerFaceAlpha',0.25','jitter','on','jitterAmount',0.15,'MarkerFaceColor',coleanticevic_rgb_networks(ceil(net/2),:)/255);
end
boxplot(outbound_bpdata_mat,'Symbol','w+');
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),glasser_rgb_networks_flip(ceil(j/2),:)/255,'FaceAlpha',.5);
end

xticks(1.5:2:24); xticklabels(network_names); xtickangle(45); xlabel('Network');
ylabel('Average Relationship Between RSFC and ISC (Spearman Rho)'); title('Average of Outbound RSFC-ISC Maps for Each Parcel, by Resting State Network')

set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')
ylim([-.08 .08])

%% Run bootstrap resampling to test mean inbound/outbound network-level effects
num_bootstraps     = 1000;
inbound_bootstrap  = cell(num_bootstraps,1);
outbound_bootstrap = cell(num_bootstraps,1);
 
s = RandStream('mlfg6331_64','Seed',2);
options = statset('UseParallel',true, ...
    'Streams',s,'UseSubstreams',true);
 
for iter = 1:num_bootstraps
    boot = randsample(s,all_used,length(all_used),true);
    % Initialize matrices
    X                 = zeros(num_node);
    edges             = (num_node*(num_node-1))/2;
    reshaped_boostrap = zeros(edges, length(all_used));
    
    output_matrix_r_bootstrap        = zeros(num_node,edges,2);
    square_matrix_bootstrap          = NaN(num_node,num_node,num_node,2);
    mean_FD_uncensored_cat_bootstrap = mean_FD_uncensored_cat(boot,:);
    
    % Get square indices
    [i,j] = find(tril(ones(num_node), -1));
    
    for isc_scan = 1:2
        hcp_isc_bootstrap = hcp_isc(boot,:,:);
        hcp_fc_bootstrap  = hcp_fc{isc_scan,1}(:,:,boot);
        reshaped_boostrap = zeros(edges, length(all_used));
        for subj = 1:length(all_used)
            X                         = squeeze(hcp_fc_bootstrap(:,:,subj));
            reshaped_boostrap(:,subj) = X(logical(tril(ones(size(X)),-1)));
        end
        
        fd_cov = mean_FD_uncensored_cat_bootstrap(:,[fc_scan,isc_scan+2]);
        parfor node = 1:num_node
            [output_matrix_r_bootstrap(node,:,isc_scan),~] = partialcorr(squeeze(hcp_isc_bootstrap(:,node,isc_scan)),reshaped_boostrap',fd_cov,'rows','complete','type','spearman') ;
        end
    end
    
    for isc_scan = 1:2
        for node = 1:num_node
            for index = 1:length(i)
                square_matrix_bootstrap(node,i(index),j(index),isc_scan) = output_matrix_r_bootstrap(node,index,isc_scan);
                square_matrix_bootstrap(node,j(index),i(index),isc_scan) = output_matrix_r_bootstrap(node,index,isc_scan);
            end
        end
    end
    
    outbound_mat_bootstrap = zeros(num_node,num_node,2);
    inbound_mat_bootstrap  = zeros(num_node,num_node,2);
    
    % Fill matrices
    for day_id = 1:num_days
        for node1 = 1:num_node
            for node2 = 1:num_node
                outbound_mat_bootstrap(node1,node2,day_id) = squeeze(square_matrix_bootstrap(node2,node1,node2,day_id));
                inbound_mat_bootstrap(node1,node2,day_id)  = squeeze(square_matrix_bootstrap(node1,node1,node2,day_id));
            end
        end
    end
    inbound_mean_net  = zeros(num_net,num_days);
    outbound_mean_net = zeros(num_net,num_days);
    
    mean_outbound_noabs_bootstrap = conv_z2r(squeeze(nanmean(conv_r2z(outbound_mat_bootstrap),2)));
    mean_inbound_noabs_bootstrap  = conv_z2r(squeeze(nanmean(conv_r2z(inbound_mat_bootstrap),2)));
    for day_id = 1:2
        for net = 1:num_net
            inbound_mean_net(net,day_id)  = conv_z2r(nanmean(conv_r2z(mean_inbound_noabs_bootstrap(netassignments_new==net,day_id))));
            outbound_mean_net(net,day_id) = conv_z2r(nanmean(conv_r2z(mean_outbound_noabs_bootstrap(netassignments_new==net,day_id))));
        end
    end
    inbound_bootstrap{iter}  = inbound_mean_net;
    outbound_bootstrap{iter} = outbound_mean_net;
end
%% Evaluate results

% Reshape to matrix
outbound_out_mat = zeros(num_net,num_days,num_bootstraps);
inbound_out_mat  = zeros(num_net,num_days,num_bootstraps);

for net = 1:num_net
    for day_id = 1:num_days
        for iter = 1:num_bootstraps
            outbound_out_mat(net,day_id,iter) = outbound_bootstrap{iter}(net,day_id);
            inbound_out_mat(net,day_id,iter)  = inbound_bootstrap{iter}(net,day_id);
        end
    end
end

% Generate confidence intervals
conf_outbound = zeros(num_net,num_days,2);
conf_inbound  = zeros(num_net,num_days,2);

p     = (1-(.05/12)) * 100;
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));

for net = 1:num_net
    for day_id = 1:num_days
        conf_outbound(net,:,day_id) = conv_z2r(CIFcn(conv_r2z(outbound_out_mat(net,day_id,:)),p));
        conf_inbound(net,:,day_id)  = conv_z2r(CIFcn(conv_r2z(inbound_out_mat(net,day_id,:)),p));
    end
end

%% Get test/re-test reliability on outbound and inbound maps for each parcel
outbound_trt = zeros(num_node,1);
inbound_trt  = zeros(num_node,1);

for node = 1:num_node
    outbound_trt(node,1) = corr(outbound_mat(node,:,1)',outbound_mat(node,:,2)','rows','complete','type','spearman');
    inbound_trt(node,1)  = corr(inbound_mat(node,:,1)',inbound_mat(node,:,2)','rows','complete','type','spearman');
end

% Make scalars
if create_scalars == 1
    group_name = 'hcp_combined';
    
    file_name     = strcat(parc_name,'_',group_name,'_outboundTRTgsronlyages');
    write_vec     = outbound_trt;
    hcp_vis(write_vec,scalar_template,file_name);
    
    file_name     = strcat(parc_name,'_',group_name,'_inboundTRTgsronlyages');
    write_vec     = inbound_trt;
    hcp_vis(write_vec,scalar_template,file_name);
end
%% Make scatterplots for FEF/TPOJ2 outbound and inbound maps

% FEF inbound
node1 = 228;
node2 = 190;
[r,p] = partialcorr(squeeze(hcp_fc{1,1}(node1,node2,:)),squeeze(hcp_isc(:,node2,1)),mean_FD_uncensored_cat(:,[1,3]),'rows','complete','type','spearman');

mdl1 = fitlm(tiedrank(mean_FD_uncensored_cat(:,[1,3])),tiedrank(squeeze(hcp_fc{1,1}(node1,node2,:))));
mdl2 = fitlm(tiedrank(mean_FD_uncensored_cat(:,[1,3])),tiedrank(squeeze(hcp_isc(:,node2,1))));

% Ranked scatter
scatter(mdl1.Residuals.Raw,mdl2.Residuals.Raw,'filled'); lsline

rtxt = strcat(['Spearman Rho = ' num2str(r) ' P = ' num2str((p))]);
t1 = TextLocation(rtxt,'Location','Best'); t1.FontSize = 14;
title('Relationship between lFEF ISC and lFEF-lLIPv RSFC')
xlabel('lFEF-lLIPv RSFC | Mean FD (Residuals)'); ylabel('lFEF ISC | Mean FD (Residuals)')
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')

% Raw scatter
figure
scatter(squeeze(hcp_fc{1,1}(node1,node2,:)),squeeze(hcp_isc(:,node2,1)),'filled'); lsline

rtxt = strcat(['Spearman Rho = ' num2str(r) ' P = ' num2str((p))]);
t1 = TextLocation(rtxt,'Location','Best'); t1.FontSize = 14;
title('Relationship between lFEF ISC and lFEF-lLIPv RSFC')
xlabel('lFEF-lLIPv RSFC (Pearson r)'); ylabel('lFEF ISC (Pearson r)')
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')

% FEF outbound
node1 = 190;
node2 = 259;
[r,p] = partialcorr(squeeze(hcp_fc{1,1}(node1,node2,:)),squeeze(hcp_isc(:,node2,1)),mean_FD_uncensored_cat(:,[1,3]),'rows','complete','type','spearman');

mdl1 = fitlm(tiedrank(mean_FD_uncensored_cat(:,[1,3])),tiedrank(squeeze(hcp_fc{1,1}(node1,node2,:))));
mdl2 = fitlm(tiedrank(mean_FD_uncensored_cat(:,[1,3])),tiedrank(squeeze(hcp_isc(:,node2,1))));

% Ranked scatter
figure
scatter(mdl1.Residuals.Raw,mdl2.Residuals.Raw,'filled'); lsline

rtxt = strcat(['Spearman Rho = ' num2str(r) ' P = ' num2str((p))]);
t1 = TextLocation(rtxt,'Location','Best'); t1.FontSize = 14;
title('Relationship between lIFJa ISC and lFEF-lIFJa RSFC')
xlabel('lFEF-lIFJa RSFC | Mean FD (Residuals)'); ylabel('lIFJa ISC | Mean FD (Residuals)')
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')

% Raw scatter
figure
scatter(squeeze(hcp_fc{1,1}(node1,node2,:)),squeeze(hcp_isc(:,node2,1)),'filled'); lsline

rtxt = strcat(['Spearman Rho = ' num2str(r) ' P = ' num2str((p))]);
t1 = TextLocation(rtxt,'Location','Best'); t1.FontSize = 14;
title('Relationship between lIFJa ISC and lFEF-lIFJa RSFC')
xlabel('lFEF-lIFJa RSFC (Pearson r)'); ylabel('lIFJa ISC (Pearson r)')
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')

% TPOJ2 inbound
node1 = 220;
node2 = 320;
[r,p] = partialcorr(squeeze(hcp_fc{1,1}(node1,node2,:)),squeeze(hcp_isc(:,node2,1)),mean_FD_uncensored_cat(:,[1,3]),'rows','complete','type','spearman');

mdl1 = fitlm(tiedrank(mean_FD_uncensored_cat(:,[1,3])),tiedrank(squeeze(hcp_fc{1,1}(node1,node2,:))));
mdl2 = fitlm(tiedrank(mean_FD_uncensored_cat(:,[1,3])),tiedrank(squeeze(hcp_isc(:,node2,1))));

% Ranked scatter
figure
scatter(mdl1.Residuals.Raw,mdl2.Residuals.Raw,'filled'); lsline

rtxt = strcat(['Spearman Rho = ' num2str(r) ' P = ' num2str((p))]);
t1 = TextLocation(rtxt,'Location','Best'); t1.FontSize = 14;
title('Relationship between lTPOJ2 ISC and lTPOJ2-l24dd RSFC')
xlabel('lTPOJ2-l24dd RSFC | Mean FD (Residuals)'); ylabel('lTPOJ2 ISC | Mean FD (Residuals)')
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')

% Raw scatter
figure
scatter(squeeze(hcp_fc{1,1}(node1,node2,:)),squeeze(hcp_isc(:,node2,1)),'filled'); lsline

rtxt = strcat(['Spearman Rho = ' num2str(r) ' P = ' num2str((p))]);
t1 = TextLocation(rtxt,'Location','Best'); t1.FontSize = 14;
title('Relationship between lTPOJ2 ISC and lTPOJ2-l24dd RSFC')
xlabel('lTPOJ2-l24dd RSFC (Pearson r)'); ylabel('lTPOJ2 ISC (Pearson r)')
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')

% TPOJ2 outbound
node1 = 320;
node2 = 259;
[r,p] = partialcorr(squeeze(hcp_fc{1,1}(node1,node2,:)),squeeze(hcp_isc(:,node2,1)),mean_FD_uncensored_cat(:,[1,3]),'rows','complete','type','spearman');

mdl1 = fitlm(tiedrank(mean_FD_uncensored_cat(:,[1,3])),tiedrank(squeeze(hcp_fc{1,1}(node1,node2,:))));
mdl2 = fitlm(tiedrank(mean_FD_uncensored_cat(:,[1,3])),tiedrank(squeeze(hcp_isc(:,node2,1))));

% Ranked scatter
figure
scatter(mdl1.Residuals.Raw,mdl2.Residuals.Raw,'filled'); lsline

rtxt = strcat(['Spearman Rho = ' num2str(r) ' P = ' num2str((p))]);
t1 = TextLocation(rtxt,'Location','Best'); t1.FontSize = 14;
title('Relationship between lTE2p ISC and lTPOJ2-lIFJa RSFC')
xlabel('lTPOJ2-lIFJa RSFC | Mean FD (Residuals)'); ylabel('lIFJa ISC | Mean FD (Residuals)')
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')


% Raw scatter
figure
scatter(squeeze(hcp_fc{1,1}(node1,node2,:)),squeeze(hcp_isc(:,node2,1)),'filled'); lsline

rtxt = strcat(['Spearman Rho = ' num2str(r) ' P = ' num2str((p))]);
t1 = TextLocation(rtxt,'Location','Best'); t1.FontSize = 14;
title('Relationship between lIFJa ISC and lTPOJ2-lIFJa RSFC')
xlabel('lTPOJ2-lIFJa RSFC (Pearson r)'); ylabel('lIFJa ISC (Pearson r)')
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')

%% Check day-to-day reliability of FEF/TPOJ2 inbound/outbound maps

if kraepelin == 1  
    perm_id           = rotate_parcellation('/data/data7/spin_dir/left_out.txt','/data/data7/spin_dir/right_out.txt',10000);
    fef_inbound_trt_r = inbound_trt(190);
    fef_inbound_trt_p = perm_sphere_p(squeeze(inbound_mat(10,:,1))',squeeze(inbound_mat(10,:,2))',perm_id,'spearman');
    
    perm_id            = rotate_parcellation('/data/data7/spin_dir/left_out.txt','/data/data7/spin_dir/right_out.txt',10000);
    fef_outbound_trt_r = outbound_trt(190);
    fef_outbound_trt_p = perm_sphere_p(squeeze(outbound_mat(10,:,1))',squeeze(outbound_mat(10,:,2))',perm_id,'spearman');
    
    perm_id            = rotate_parcellation('/data/data7/spin_dir/left_out.txt','/data/data7/spin_dir/right_out.txt',10000);
    tpoj_inbound_trt_r = inbound_trt(320);
    tpoj_inbound_trt_p = perm_sphere_p(squeeze(inbound_mat(140,:,1))',squeeze(inbound_mat(140,:,2))',perm_id,'spearman');
    
    perm_id             = rotate_parcellation('/data/data7/spin_dir/left_out.txt','/data/data7/spin_dir/right_out.txt',10000);
    tpoj_outbound_trt_r = outbound_trt(320);
    tpoj_outbound_trt_p = perm_sphere_p(squeeze(outbound_mat(140,:,1))',squeeze(outbound_mat(140,:,2))',perm_id,'spearman');
end

%% Do PCA

% Zero out diagonals
for node = 1:num_node
    outbound_mat(node,node,:) = 0;
    inbound_mat(node,node,:)  = 0;
end

% Initialize matrices
coeff_outbound     = zeros(num_node,num_node-1,2); coeff_inbound = zeros(num_node,num_node-1,2);
score_outbound     = zeros(num_node,num_node-1,2); score_inbound = zeros(num_node,num_node-1,2);
explained_outbound = zeros(num_node-1,2); explained_inbound = zeros(num_node-1,2);

% Run PCA
for day_id = 1:num_days
    [coeff_outbound(:,:,day_id),score_outbound(:,:,day_id),~,~,explained_outbound(:,day_id),~] = pca(squeeze(outbound_mat(:,:,day_id)));
    [coeff_inbound(:,:,day_id),score_inbound(:,:,day_id),~,~,explained_inbound(:,day_id),~]    = pca(squeeze(inbound_mat(:,:,day_id)));
end

% Make scree plots
subplot(2,1,1)
plot(explained_outbound(1:10,1),'-ko','MarkerIndices',1:10,'MarkerFaceColor','b','MarkerEdgeColor','k')
ylabel('Percentage of Variance Explained'); xlabel('Principal Component')
pbaspect([1 1 1])
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')
ylim([0 40])

subplot(2,1,2)
plot(explained_outbound(1:10,2),'-ko','MarkerIndices',1:10,'MarkerFaceColor','b','MarkerEdgeColor','k')
ylabel('Percentage of Variance Explained'); xlabel('Principal Component')
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')
pbaspect([1 1 1])
ylim([0 40])

% Make scalars
if create_scalars == 1
    for coeffnum = 1:6
        for day_id = 1:num_days
            
            % coefficients
            write_vec  = squeeze(coeff_outbound(:,coeffnum,day_id));
            group_name = strcat('hcp',string(day_id));
            type       = strcat(string(coeffnum),'coeff_outbound');
            file_name  = strcat(parc_name,'_',group_name,'_',type);
            hcp_vis(write_vec,scalar_template,file_name);
            
            write_vec  = squeeze(coeff_inbound(:,coeffnum,day_id));
            group_name = strcat('hcp',string(day_id));
            type       = strcat(string(coeffnum),'coeff_inbound');
            file_name  = strcat(parc_name,'_',group_name,'_',type);
            hcp_vis(write_vec,scalar_template,file_name);
            
            % scores
            write_vec  = squeeze(score_outbound(:,coeffnum,day_id));
            group_name = strcat('hcp',string(day_id));
            type       = strcat(string(coeffnum),'score_outbound');
            file_name  = strcat(parc_name,'_',group_name,'_',type);
            hcp_vis(write_vec,scalar_template,file_name);
            
            write_vec  = squeeze(score_inbound(:,coeffnum,day_id));
            group_name = strcat('hcp',string(day_id));
            type       = strcat(string(coeffnum),'score_inbound');
            file_name  = strcat(parc_name,'_',group_name,'_',type);
            hcp_vis(write_vec,scalar_template,file_name);
        end
    end
end

% Do PCA using full 360x360 RSFC matrix for inbound analysis
square_matrix_full_inbound               = (square_matrix_r(:,:,:,[1,4]));
num_edge                      = (num_node*(num_node-1))/2;
square_matrix_z_coleanticevic = zeros(num_node,num_edge,2);
for day_id = 1:num_days
    for node = 1:num_node
        X = squeeze(square_matrix_full_inbound(node,:,:,day_id));
        square_matrix_z_coleanticevic(node,:,day_id) = X(logical(tril(ones(size(X)),-1)));
    end
end

coeff_vector     = zeros(num_edge,num_node-1,2);
score_vector     = zeros(num_node,num_node-1,2);
explained_square = zeros(num_node-1,2);

for day_id = 1:num_days
    [coeff_vector(:,:,day_id),score_vector(:,:,day_id),latent,tsquared,explained_square(:,day_id),mu] = pca(squeeze(square_matrix_z_coleanticevic(:,:,day_id)));
end

% Reshape to square
num_coef     = 6;
coeff_square = zeros(num_coef,num_node,num_node,2);
score_square = zeros(num_coef,num_node,num_node,2);

[i,j] = find(tril(ones(num_node), -1)); 
for day_id = 1:num_days
    for co = 1:num_coef
        for index = 1:num_edge
            coeff_square(co,i(index),j(index),day_id) = coeff_vector(index,co,day_id);
            coeff_square(co,j(index),i(index),day_id) = coeff_vector(index,co,day_id);
        end
    end
end

% Make full inbound PCA pscalars
if create_scalars == 1
    for coeffnum = 1:num_coef
        for day_id = 1:num_days
            
            % scores
            write_vec  = squeeze(score_vector(:,coeffnum,day_id));
            group_name = strcat('hcp',string(day_id));
            type       = strcat(string(coeffnum),'square_score_inbound');
            file_name  = strcat(parc_name,'_',group_name,'_',type);
            hcp_vis(write_vec,scalar_template,file_name);
        end
    end
end

% Re-order to match Cole/Anticevic RSNs
coeff_square_coleanticevic = coeff_square;
for day_id = 1:num_days
    for co = 1:num_coef
        for node1 = 1:num_node
            for node2 = 1:num_node
                coeff_square_coleanticevic(co,node1,node2,day_id) = coeff_square(co,indices_reshape(node1,1),indices_reshape(node2,1),day_id);
            end
        end
    end
end

% Create heatmap to visualize PC loadings
for coeff_num = 1:2
    
    % Set the same scale for both heatmaps
    curr_img1 = squeeze(coeff_square_coleanticevic(coeff_num,:,:,1));
    curr_img2 = squeeze(coeff_square_coleanticevic(coeff_num,:,:,2));
    if max(max(cat(1,curr_img1,curr_img2))) > abs(min(min(cat(1,curr_img1,curr_img2))))
        use_scale = max(max(cat(1,curr_img1,curr_img2))) ;
    else
        use_scale = abs(min(min(cat(1,curr_img1,curr_img2)))) ;
    end
    
    day_id = 1;
    figure
    subplot(1,2,1)
    curr_img = squeeze(coeff_square_coleanticevic(coeff_num,:,:,day_id));
    
    imagesc(curr_img)
    colormap(fireice(1000))
    title(strcat(['PC ' string(coeff_num) ' Day ' string(day_id)]))
    ticks = [1,7,63,100,156,179,202,252,267,344,351,355];
    xticks(ticks); yticks(ticks)
    network_names = {'Visual1','Visual2','Somatomotor','Cingulo-Opercular','Dorsal-attention','Language','Frontoparietal','Auditory','Default','Posterior-Multimodal','Ventral-Multimodal','Orbito-Affective'};
    yticklabels(network_names); xticklabels(network_names); xtickangle(45)
    colorbar; caxis([-use_scale,use_scale])
    daspect([1 1 1 ])
    
    day_id = 2;
    subplot(1,2,2)
    curr_img = squeeze(coeff_square_coleanticevic(coeff_num,:,:,day_id));

    % Flip axes for consistency
    if coeff_num == 1
        curr_img = -curr_img;
    end
    
    imagesc(curr_img)
    colormap(fireice(256))
    title(strcat(['PC ' string(coeff_num) ' Day ' string(day_id)]))
    ticks = [1,7,63,100,156,179,202,252,267,344,351,355];
    xticks(ticks); yticks(ticks)
    network_names = {'Visual1','Visual2','Somatomotor','Cingulo-Opercular','Dorsal-attention','Language','Frontoparietal','Auditory','Default','Posterior-Multimodal','Ventral-Multimodal','Orbito-Affective'};
    yticklabels(network_names); xticklabels(network_names); xtickangle(45)
    colorbar; caxis([-use_scale,use_scale])
    counter = counter+1;
    daspect([1 1 1 ])
end

% Make scree plots
subplot(2,1,1)
plot(explained_square(1:10,1),'-ko','MarkerIndices',1:10,'MarkerFaceColor','b','MarkerEdgeColor','k')
ylabel('Percentage of Variance Explained'); xlabel('Principal Component')
pbaspect([1 1 1])
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')
ylim([0 15])
subplot(2,1,2)
plot(explained_square(1:10,2),'-ko','MarkerIndices',1:10,'MarkerFaceColor','b','MarkerEdgeColor','k')
ylabel('Percentage of Variance Explained'); xlabel('Principal Component')
pbaspect([1 1 1])
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')
ylim([0 10])

%% Topology

% First, do with all subjects to generate heatmaps
sim_mat_isc = NaN(length(all_used),length(all_used),2);
sim_mat_fc  = NaN(length(all_used),length(all_used),2);
for day_id=1:2
    tmp_isc = squeeze(hcp_isc(:,:,day_id))';
    tmp_fc  = reshape(hcp_fc{day_id,1},[num_node*num_node,num_subj]);
    
    sim_mat_isc(:,:,day_id) = corr(tmp_isc(:,all_used),'rows','complete','type','spearman');
    sim_mat_fc(:,:,day_id) = corr(tmp_fc(:,all_used),'rows','complete','type','spearman');
end

for day_id = 1:2
    for subj = 1:length(all_used)
        sim_mat_isc(subj,subj,day_id) = NaN;
        sim_mat_fc(subj,subj,day_id)  = NaN;
    end
end

detect_outliers = squeeze(nanmean(sim_mat_isc(:,:,:),2));
detect_outliers(:,3) = 1:176;
crit = nanmean(detect_outliers) - 2.5*std(detect_outliers);
detect_outliers(:,1:2)
colormin       = [50 50 200];
colormax       = [200 50 50];
colorneutral   = [255 255 255];

n    = 100;

Rlow  = linspace(colormin(1)/255,colorneutral(1)/255,n);
Blow  = linspace(colormin(3)/255,colorneutral(3)/255,n);
Glow  = linspace(colormin(2)/255,colorneutral(2)/255,n);
Rhigh = linspace(colorneutral(1)/255,colormax(1)/255,n);
Bhigh = linspace(colorneutral(3)/255,colormax(3)/255,n);
Ghigh = linspace(colorneutral(2)/255,colormax(2)/255,n);
colormap( [cat(1,Rlow(:),Rhigh(:)), cat(1,Glow(:),Ghigh(:)), cat(1,Blow(:),Bhigh(:))] );  %// create colormap
c     = colorbar; c.Location = 'southoutside'; set(c,'xtick',[]);

% Generate RSFC/ISC similarity heatmaps
for day_id = 1:2
    subplot(2,2,day_id)
    imagesc(sim_mat_fc(:,:,day_id))
    colormap(fireice)
    caxis([-1 1])
    
    subplot(2,2,day_id+2)
    imagesc(sim_mat_isc(:,:,day_id))
    colormap()
    caxis([-1 1])
end

% These subjects are marked outliers on the ISC similarity heatmaps, so
% exclude them from this analysis
exclude1 = [88, 106];
exclude2 = [34, 94];

all_used_topo(:,1) = setdiff(all_used,exclude1);
all_used_topo(:,2) = setdiff(all_used,exclude2);


% Are subjects with more similar FC patterns more similar in their ISC
% patterns as well?

% Initialize matrices
sim_mat_isc = NaN(length(all_used_topo),length(all_used_topo),2);
sim_mat_fc  = NaN(length(all_used_topo),length(all_used_topo),2);

% Calculate RSFC and ISC similarity between all pairs of subjects
for day_id=1:2
    tmp_isc = squeeze(hcp_isc(:,:,day_id))';
    tmp_fc  = reshape(hcp_fc{day_id,1},[num_node*num_node,num_subj]);
    
    sim_mat_isc(:,:,day_id) = corr(tmp_isc(:,all_used_topo(:,day_id)),'rows','complete','type','spearman');
    sim_mat_fc(:,:,day_id)  = corr(tmp_fc(:,all_used_topo(:,day_id)),'rows','complete','type','spearman');
end

% Calculate motion similarity between all pairs of subjects (higher = more
% similar)
sim_mat_motion = zeros(num_subj,num_subj,4);
for scan_id = 1:4
    for subj1 = 1:num_subj
        for subj2 = 1:num_subj
            if subj1==subj2
                sim_mat_motion(subj1,subj2,scan_id) = NaN;
            else
                sim_mat_motion(subj1,subj2,scan_id) = -abs(mean_FD_uncensored_cat(subj1,scan_id)-mean_FD_uncensored_cat(subj2,scan_id));
            end
        end
    end
end
sim_mat_motion = cat(3,sim_mat_motion(all_used_topo(:,1),all_used_topo(:,1),1),sim_mat_motion(all_used_topo(:,2),all_used_topo(:,2),2),sim_mat_motion(all_used_topo(:,1),all_used_topo(:,1),3),sim_mat_motion(all_used_topo(:,2),all_used_topo(:,2),4));

% Calculate mean ISC similarity between all pairs of subjects (higher = more
% similar)
sim_mat_mean_isc = zeros(num_subj,num_subj,2);
for scan_id = 1:2
    for subj1 = 1:num_subj
        for subj2 = 1:num_subj
            if subj1==subj2
                sim_mat_mean_isc(subj1,subj2,scan_id) = NaN;
            else
                sim_mat_mean_isc(subj1,subj2,scan_id) = -abs(conv_r2z(mean_isc(subj1,scan_id))-conv_r2z(mean_isc(subj2,scan_id)));
            end
        end
    end
end
sim_mat_mean_isc = cat(3,sim_mat_mean_isc(all_used_topo(:,1),all_used_topo(:,1),1),sim_mat_mean_isc(all_used_topo(:,2),all_used_topo(:,2),2));

% Get means/stds for similarity matrices
rsfc_sim_mean = conv_z2r(nanmean(nanmean(conv_r2z(sim_mat_fc))));
isc_sim_mean = conv_z2r(nanmean(nanmean(conv_r2z(sim_mat_isc))));

rsfc_sim_std = conv_z2r(nanstd(nanstd(conv_r2z(sim_mat_fc))));
isc_sim_std = conv_z2r(nanstd(nanstd(conv_r2z(sim_mat_isc))));

% Get coordinates of upper triangle
square     = ones(length(all_used_topo));
square_upp = tril(square,-1);
upp_id     = find(square_upp);

% Set diagonals to NaN
for day_id = 1:2
    for subj = 1:length(all_used_topo)
        sim_mat_motion(subj,subj,day_id)   = NaN;
        sim_mat_mean_isc(subj,subj,day_id) = NaN;
        sim_mat_isc(subj,subj,day_id)      = NaN;
        sim_mat_fc(subj,subj,day_id)       = NaN;
    end
end

topo_r = zeros(2,1);
for day_id = 1:2
    tmp_isc          = sim_mat_isc(:,:,day_id);
    tmp_mean_isc     = sim_mat_mean_isc(:,:,day_id);
    tmp_fc           = sim_mat_fc(:,:,day_id);
    tmp_motion_fc    = sim_mat_motion(:,:,day_id);
    tmp_motion_isc   = sim_mat_motion(:,:,day_id+2);
    topo_r(day_id,1) = partialcorr(tmp_isc(upp_id),tmp_fc(upp_id),cat(2,tmp_motion_fc(upp_id),tmp_motion_isc(upp_id),tmp_mean_isc(upp_id)),'rows','complete','type','spearman');
    
    mdl1  = fitlm(cat(2,tmp_motion_fc(upp_id),tmp_motion_isc(upp_id),tmp_mean_isc(upp_id)),tmp_isc(upp_id));
    mdl2  = fitlm(cat(2,tmp_motion_fc(upp_id),tmp_motion_isc(upp_id),tmp_mean_isc(upp_id)),tmp_fc(upp_id));
    
    %mdl1_ranked  = fitlm(cat(2,tiedrank(tmp_motion_fc(upp_id)),tiedrank(tmp_motion_isc(upp_id))),tiedrank(tmp_isc(upp_id)));
    %mdl2_ranked  = fitlm(cat(2,tiedrank(tmp_motion_fc(upp_id)),tiedrank(tmp_motion_isc(upp_id))),tiedrank(tmp_fc(upp_id)));
    
    % Now predict gISC similarity controlling for topology
    %topo_r(day_id,1) = partialcorr(tmp_mean_isc(upp_id),tmp_fc(upp_id),cat(2,tmp_motion_fc(upp_id),tmp_motion_isc(upp_id),tmp_isc(upp_id)),'rows','complete','type','spearman')
    
    %mdl1  = fitlm(cat(2,tmp_motion_fc(upp_id),tmp_motion_isc(upp_id),tmp_isc(upp_id)),tmp_mean_isc(upp_id));
    %mdl2  = fitlm(cat(2,tmp_motion_fc(upp_id),tmp_motion_isc(upp_id),tmp_mean_isc(upp_id)),tmp_fc(upp_id));
    
    
    % Create scatter of day1/day2 gISC
    figure
    % Create scatter with raw values
    scatter((tmp_isc(upp_id)),(tmp_fc(upp_id)),5,'filled'); lsline
    set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')
    xlabel('ISC topological similarity')
    ylabel('RSFC similarity')
    
    % Create scatter with residuals
    figure
    scatter(mdl1.Residuals.Raw,mdl2.Residuals.Raw,5,'filled'); lsline
    xlim([-.5 .3])
    set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')
    xlabel('ISC topological similarity (Residuals)')
    ylabel('RSFC similarity (Residuals)')
end

% Assess significance of RSFC/ISC similarity relationship with permutation
% testing

% Initialize matrices
isc_all_used                    = cat(3,squeeze(hcp_isc(all_used_topo(:,1),:,1))',squeeze(hcp_isc(all_used_topo(:,2),:,2))');
mean_FD_uncensored_cat_all_used = cat(3,mean_FD_uncensored_cat(all_used_topo(:,1),1),mean_FD_uncensored_cat(all_used_topo(:,2),2),mean_FD_uncensored_cat(all_used_topo(:,1),3),mean_FD_uncensored_cat(all_used_topo(:,2),4));
mean_isc_all_used               = cat(2,mean_isc(all_used_topo(:,1),1),mean_isc(all_used_topo(:,2),2));
num_perm                        = 10000;
reshape_out                     = zeros(num_perm,length(all_used_topo));

sss        = RandStream.create('mlfg6331_64','Seed',3,'NumStreams',12);
paroptions = statset('UseParallel',true, ...
    'Streams',sss,'UseSubstreams',true);

perm_mat = zeros(num_perm,2);
for perm = 1:num_perm
    perm
    sim_mat_isc_perm = NaN(length(all_used_topo),length(all_used_topo),2);
    reshape_vec      = randperm(sss,length(all_used_topo));
    
    for day_id = 1:2
        % Calculate ISC similarity across shuffled subjects
        tmp_isc                      = squeeze(isc_all_used(:,:,day_id))';
        tmp_isc                      = tmp_isc(:,reshape_vec);
        sim_mat_isc_perm(:,:,day_id) = corr(tmp_isc,'rows','complete','type','spearman');
        
        % Initialize matrices
        tmp_mean_isc   = zeros(length(all_used_topo));
        tmp_motion_isc = zeros(length(all_used_topo));
        
        % Re-calculate mean ISC and motion similarity across shuffled subjects
        for subj1 = 1:length(all_used_topo)
            for subj2 = 1:length(all_used_topo)
                if subj1==subj2
                    tmp_motion_isc(subj1,subj2) = NaN;
                    tmp_mean_isc(subj1,subj2)   = NaN;
                else
                    tmp_motion_isc(subj1,subj2) = -abs(mean_FD_uncensored_cat_all_used(subj1,day_id+2)-mean_FD_uncensored_cat_all_used(subj2,day_id+2));
                    tmp_mean_isc(subj1,subj2)   = -abs(conv_r2z(mean_isc_all_used(subj1,scan_id))-conv_r2z(mean_isc_all_used(subj2,scan_id)));
                end
            end
        end
        
        % Reshape matrices
        tmp_fc        = sim_mat_fc(:,:,day_id);
        tmp_motion_fc = sim_mat_motion(:,:,day_id);
        fc_vector     = cat(2,tmp_fc(upp_id),tmp_motion_fc(upp_id));
        isc_perm      = cat(2,sim_mat_isc_perm(upp_id),tmp_motion_isc(upp_id));
        
        % Correlate shuffled ISC matrix with RSFC matrices
        perm_mat(perm,day_id) = partialcorr(fc_vector(:,1),isc_perm(:,1),cat(2,fc_vector(:,2),isc_perm(:,2),tmp_mean_isc(upp_id)),'rows','complete','type','spearman');
    end
end

% Calculate permutation p-statistic
topo_p(1,1) = (1 + sum(perm_mat(:,1)>topo_r(1,1)))/(num_perm+1);
topo_p(1,2) = (1 + sum(perm_mat(:,2)>topo_r(2,1)))/(num_perm + 1);

% Compare null distribution from permutation tests with actual results
clear h

subplot(1,2,1)
X          = perm_mat(:,1);
color      = [.5 .5 .5];
[f, Xi, ~] = ksdensity(X, 'bandwidth', []);
h{1}       = area(Xi, f); hold on
set(h{1}, 'FaceColor', color); set(h{1}, 'EdgeColor', [0 0 0]); set(h{1}, 'LineWidth', 2); set(h{1}, 'FaceAlpha', 1);
xline(median(topo_r(1,1)),'LineWidth',2,'Color',[1,0,0])
yl = get(gca, 'YLim');
set(gca, 'YLim', [0 yl(2)]); set(gca, 'XLim', [-.3 .3]);
rtxt = strcat(['Correlation between RSFC and ISC Similarity = ' num2str(topo_r(1,1)) ' permutation p = ' num2str(topo_p(1,1))]);
t1 = TextLocation(rtxt,'Location','Best'); t1.FontSize = 14;
set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')
ylim([0 9])
xlim([-.3 .3])


subplot(1,2,2)
X          = perm_mat(:,2);
color      = [.5 .5 .5];
[f, Xi, u] = ksdensity(X, 'bandwidth', []);
h{1}       = area(Xi, f); hold on
set(h{1}, 'FaceColor', color); set(h{1}, 'EdgeColor', [0 0 0]); set(h{1}, 'LineWidth', 2); set(h{1}, 'FaceAlpha', 1);
xline(median(topo_r(2,1)),'LineWidth',2,'Color',[1,0,0])
yl = get(gca, 'YLim');
set(gca, 'YLim', [0 yl(2)]);
set(gca, 'XLim', [-.3 .3]);
rtxt = strcat(['Correlation between RSFC and ISC Similarity = ' num2str(topo_r(2,1)) ' permutation p = ' num2str(topo_p(1,2))]);
t1   = TextLocation(rtxt,'Location','Best'); t1.FontSize = 14;
ylim([0 9])
xlim([-.3 .3])

set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')

%% Pipeline figures

% RSFC timecourses
plot(smoothdata(sub_data_cat{1,1}(2,200:300,1),'movmedian',2))
plot(smoothdata(sub_data_cat{1,1}(3,200:300,1),'movmedian',2))

% ISC timecourses
plot(smoothdata(sub_data_cat{3,1}(2,200:300,1),'movmedian',2))
plot(smoothdata(nanmean(squeeze(sub_data_cat{3,1}(2,200:300,2:4)),2),'movmedian',2))

% FC square
demo_fc = hcp_fc{1,1}(:,:,1);
demo_fc_reshaped = zeros(num_node,num_node);
for node1 = 1:num_node
    for node2 = 1:num_node
        if node1 == node2
            demo_fc_reshaped(node1,node2) = 0;
        else
            demo_fc_reshaped((node1),(node2)) = demo_fc(indices_reshape(node1),indices_reshape(node2));
        end
    end
end

imagesc(demo_fc_reshaped);colormap(fireice(1000))

if max(max(demo_fc_reshaped)) > abs(min(min(demo_fc_reshaped)))
    use_scale = max(max(demo_fc_reshaped)) ;
else
    use_scale = abs(min(min(demo_fc_reshaped))) ;
end

set(gca,'xtick',[]); set(gca,'xticklabel',[])
set(gca,'ytick',[]); set(gca,'yticklabel',[])
caxis([-use_scale,use_scale]); axis square

% ISC square
demo_isc = squeeze(hcp_isc(all_used,:,1));
demo_isc_reshaped = zeros(length(all_used),num_node);
for node1 = 1:num_node
    demo_isc_reshaped(:,node1) = demo_isc(:,indices_reshape(node1));
end

if max(max(demo_isc_reshaped)) > abs(min(min(demo_isc_reshaped)))
    use_scale = max(max(demo_isc_reshaped));
else
    use_scale = abs(min(min(demo_isc_reshaped)));
end

imagesc(demo_isc_reshaped);colormap(fireice(1000))
set(gca,'xtick',[]); set(gca,'xticklabel',[])
set(gca,'ytick',[]); set(gca,'yticklabel',[])
caxis([-use_scale,use_scale]); axis square

% GISC square
demo_gisc_reshaped = nanmean(demo_isc_reshaped,2);

imagesc(demo_gisc_reshaped);colormap(fireice(1000))
set(gca,'xtick',[]); set(gca,'xticklabel',[])
set(gca,'ytick',[]); set(gca,'yticklabel',[])
caxis([-max(demo_gisc_reshaped),max(demo_gisc_reshaped)])

% Ridge predictions
condensed_isc = conv_z2r(ridgeCPMoutput.observed_ISC(:,1));
condensed_isc(condensed_isc<.2) = condensed_isc(condensed_isc<.2) + .2;
scatter(conv_z2r(ridgeCPMoutput.predicted_ISC'),condensed_isc,'filled'); lsline
r_rank = corr(condensed_isc,ridgeCPMoutput.predicted_ISC','rows','complete','type','spearman');
title('Ridge CPM Performance: Actual gISC vs. Predicted gISC')
xlabel('Predicted gISC'); ylabel('Actual gISC');
rtxt = strcat(['spearman rho = ' num2str(r_rank) ' p-value = ' num2str(p_rank)]);
t1   = TextLocation(rtxt,'Location','Best'); t1.FontSize = 14;

set(gca,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}); set(gca,'TickDir','out'); set(gca, 'DefaultTextInterpreter', 'none')

%% Save results
A = datestr(now);
curr_date = A(~isspace(A));
save(strcat('/data/data7/restsync/restsync_results_',curr_date,'.mat'),'-v7.3')