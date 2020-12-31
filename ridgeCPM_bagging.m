%%%%%%%%
% Author: David Gruskin
% Contact: dcg2153@cumc.columbia.edu
% Last updated: 12/2020
% Project: Bridging Resting State Functional Connectivity and Movie-Evoked
% Brain Activity
% Description: This code is called by restsync_main.m to run the rCPM analysis 

% Originally by: Siyuan Gao, Yale University, 2018-2019, Adapted by Abby Greene, 2019
% Adapted by DG, 2020
%%%%%%%%

%%
function   [mdl_bagged, ridgeCPMoutput ] = ridgeCPM_bagging(hcp_fc, hcp_isc, all_used, sub_data_cat, mean_FD_uncensored_cat, family_ids, num_iter, randomize)
%
%   Input:      hcp_fc,                FC matrices for all subjects
%                                      [scans x 1 cell] [parcels x parcels x
%                                      subjects]
%                                  
%
%               hcp_isc,               ISC matrices for all subjects
%                                      [subjects x parcels x days]
%
%               all_used,              Vector of subject IDs to be used 
%                                      [subjects x 1]
%
%               sub_data_cat,          Matrix containing BOLD timecourse data
%                                      [scans x 1 cell] [parcels x TRs x subjects]
%
%               mean_FD_uncensored_cat Matrix containing mean FD values 
%                                      [subjects x scans]
%
%               family_ids,            Matrix containing family IDs for
%                                      each subject
%                                      [subjects x 2]
%
%               num_iter,              Scalar reflecting number of CV iterations to perform
%
%               randomize,             Binary reflecting whether to run
%                                      actual model or null permutations
%
%   Output:     ridgeCPMoutput,        Struct containing the bagged model's
%                                      predictions as well as the CV
%                                      outputs
%                                   
%
%               mdl_bagged,            Struct containing the bagged model 
%                                      parameters        
%
%

%% Initialization

% Set hyperparameters
feature_thresh = .1; % This value determines the percentage of CV models in which a given feature must be present for it to be included in the bagged model
k              = 2; % This value determines the number of CV folds
thresh         = .01; % This is the p-value threshold for the initial feature selection step
v_alpha        = 1e-9; % This is the alpha value (0 = Ridge, 1 = Lasso)
num_lambda     = 5; % This is the number of lambdas to search 
num_loops      = 1; % This can be changed to loop through the full analysis multiple times

% Define training matrices
all_mats                   = hcp_fc{1,1}(:,:,all_used);
family_ids_tmp             = family_ids(all_used,:);
movie_timecourses          = sub_data_cat{3,:}(:,:,all_used);
mean_FD_uncensored_cat_tmp = mean_FD_uncensored_cat(all_used,:);

% Parse dimensions of the data
num_sub_total = size(all_mats, 3);
num_node      = size(all_mats, 1);

% Convert connectivity to edge matrix 
num_edge  = num_node * (num_node - 1) / 2;
all_edges = zeros(num_edge, num_sub_total);
for i_sub = 1 : num_sub_total
    all_edges(:, i_sub) = squareform(tril(all_mats(:, :, i_sub), -1));
end

all_edges = permute(all_edges, [1, 3, 2]);
all_edges = reshape(all_edges, [], num_sub_total);

% Initialize matrices
predicted_isc = zeros(num_loops,num_sub_total);
mse_output    = zeros(num_loops,3);
qs_output     = zeros(num_loops,1);
rank_p_output = zeros(num_loops,1);
rank_r_output = zeros(num_loops,1);
spearman_out  = zeros(num_loops,num_iter);

% Set randomization seed for replicability (same random numbers each time
% this script is run)
s       = RandStream('mlfg6331_64','Seed',1);
options = statset('UseParallel',true, ...
    'Streams',s,'UseSubstreams',true);


for loop = 1:num_loops
    disp('Loop #' + string(loop))
    
    % Initialize model outputs
    intercept_out = zeros(num_iter,k);
    lambda_out    = zeros(num_iter,k);
    edge_idx_out  = cell(num_iter,k);
    betas_out     = cell(num_iter,k);
   
    parfor iter = 1:num_iter
        disp('Iteration #' + string(iter))
        
        % Split the sample into two groups (accounting for family
        % structure)
        randomsplit      = zeros(176,2);
        randomsplit(:,1) = 1:176;
        indices          = kfold_family(family_ids_tmp,2,(iter+(num_loops*(loop-1))));
        randomsplit(:,2) = indices;
        
        % Calculate LOO ISC in the training and test groups separately
        hcp_isc1 = NaN(176,num_node);
        for subj = 1:176
            tmp = movie_timecourses;
            group_id = randomsplit(subj,2);
            if ~isnan(group_id)
                targ             = squeeze(movie_timecourses(:,:,subj))';
                tmp(:,:,subj)    = NaN;
                tmp(:,:,randomsplit(:,2)~=group_id) = NaN;
                tmp_mean         = nanmean(tmp,3)';
                hcp_isc1(subj,:) = diag(corr(targ(:,:),tmp_mean(:,:),'type','pearson','rows','complete'));
            end
        end
        
        % Create concatenated behavioral matrix (ISC and motion covariates)
        % for later use
        all_behav1 = cat(2,nanmean(conv_r2z(hcp_isc1(:,:,1)),2),mean_FD_uncensored_cat_tmp(:,[1,3])); 
        if randomize == 1
            all_behav_tmp = all_behav1(randperm(s,size(all_behav1,1)),:);
            
            % Make sure the mean FD values from the rest scans are kept unshuffled
            all_behav_tmp(:,2) = mean_FD_uncensored_cat_tmp(:,1);

        else
            all_behav_tmp = all_behav1;
        end
        
        
        %% Ridge CPM
        y = zeros(num_sub_total, 1);
        for i_fold = 1 : k
            
            test_idx  = (indices==i_fold);
            train_idx = (indices~=i_fold);
            
            train_mats = all_edges(:, train_idx);
            test_mats  = all_edges(:, test_idx);
            
            train_behav             = all_behav_tmp;
            train_behav(test_idx,:) = [];
            
            % Feature selection step
            [~, edge_p] = partialcorr(train_mats', train_behav(:,1),train_behav(:,2:end),'rows','complete','type','spearman');
            edge_idx    = find(edge_p < thresh);
            
            % Find optimal value for lambda hyperparameter
            lmax = computeLambdaMax(cat(2,train_mats(edge_idx, :)',train_behav(:,2:end)), train_behav(:,1), [],...
                0.01, true);
            lmin           = lmax * 0.001;
            loghi          = log(lmax);
            loglo          = log(lmin);
            lambda_squence = exp(linspace(loghi,loglo,num_lambda));
            
            opt_lambda     = find_lambda(train_mats, train_behav,...
                lambda_squence, thresh); % TIME-CONSUMING            
            
            % Fit the ridge model
            [beta, fit_info] = lasso(cat(2,train_mats(edge_idx, :)',train_behav(:,2:end)),...
                train_behav(:,1), 'Alpha',v_alpha, 'Lambda', opt_lambda);
           
            intercept = fit_info.Intercept;
            beta      = beta(1:size(beta,1)-2,:);
            
            % Run the model on test subjects
            y(test_idx) = test_mats(edge_idx,:)'*beta+intercept;
            
            % Store model outputs
            intercept_out(iter,i_fold) = intercept;
            lambda_out(iter,i_fold)    = opt_lambda;
            edge_idx_out{iter,i_fold}  = edge_idx;
            betas_out{iter,i_fold}     = beta;
        end
        
        % Calculate model fit for this train/test split and save output
        [r_spearman, ~]         = partialcorr(y, all_behav_tmp(:,1),all_behav_tmp(:,2:end),'type','spearman');
        spearman_out(loop,iter) = r_spearman;
    end
    
    %% Building the bagged model
    
    % Incorporate all features that appear in 10% of the models trained
    % above into the bagged model
    betas_perc = NaN(num_edge,num_iter*k);
    
    counter = 1;
    for iter = 1:num_iter
        for fold = 1:k
            for feat = 1:length(edge_idx_out{iter,fold})
                betas_perc(edge_idx_out{iter,fold}(feat,1),counter) = betas_out{iter,fold}(feat,1);
            end
            counter = counter +1;
        end
    end
    
    % Get the average ridge coefficient value for the previously identified
    % features
    betas_mean    = nanmean(betas_perc,2);
    betas_counter = zeros(size(betas_perc,1),1);
    for row = 1:size(betas_perc,1)
        betas_counter(row,1) = length(find(betas_perc(row,:)>0));
    end
    
    betas_counter = betas_counter / (num_iter*k);
    
    % Store the bagged model parameters
    [edge_idx_bagged,~] = find(betas_counter > feature_thresh);
    betas_bagged        = betas_mean(edge_idx_bagged);
    intercept_bagged    = nanmean(nanmean(intercept_out(1:iter,:)));
    lambda_bagged       = nanmean(nanmean(lambda_out(1:iter,:)));
    
    mdl_bagged           = struct;
    mdl_bagged.edge_idx  = edge_idx_bagged;
    mdl_bagged.betas     = betas_bagged;
    mdl_bagged.intercept = intercept_bagged;
    mdl_bagged.lambda    = lambda_bagged;
    
    % Isolate the (Day 2) test data
    all_mats2  = hcp_fc{2,1};
    all_behav2 = cat(2,nanmean(conv_r2z(hcp_isc(:,:,2)),2),mean_FD_uncensored_cat(:,[2,4]));
    
    nan_idx2                = find(isnan(all_behav2(:,1)));
    all_behav2(nan_idx2,:)  = [];
    all_mats2(:,:,nan_idx2) = [];
    
    ridgeCPMoutput              = struct;
    ridgeCPMoutput.observed_ISC = all_behav2;
    
    all_edges2 = zeros(num_edge, num_sub_total);
    for i_sub = 1 : num_sub_total
        all_edges2(:, i_sub) = squareform(tril(all_mats2(:, :, i_sub), -1));
    end
    
    test_mats2 = reshape(permute(all_edges2, [1, 3, 2]), [], num_sub_total);
    
    % Predict Day 2 gISC by feeding Day 2 RSFC into the bagged model
    y = test_mats2(mdl_bagged.edge_idx, :)'*mdl_bagged.betas+mdl_bagged.intercept;
    
    % Evaluate predictive performance of the bagged model
    [r_rank, p_rank] = partialcorr(y, all_behav2(:,1),all_behav2(:,2:end), 'type', 'spearman','rows','complete');
    
    mse = sum((y - all_behav2).^2) / num_sub_total;
    q_s = 1 - mse / var(all_behav2, 1);
    
    % Store outputs
    rank_p_output(loop,1) = p_rank;
    rank_r_output(loop,1) = r_rank;
    qs_output(loop,1)     = q_s;
    mse_output(loop,:)    = mse;
    predicted_isc(loop,:) = y;
    
end

% Store outputs
ridgeCPMoutput.rank_r        = rank_r_output;
ridgeCPMoutput.rank_p        = rank_p_output;
ridgeCPMoutput.mse           = mse_output;
ridgeCPMoutput.qs            = qs_output;
ridgeCPMoutput.predicted_ISC = predicted_isc;
ridgeCPMoutput.spearman_out  = spearman_out;

end
