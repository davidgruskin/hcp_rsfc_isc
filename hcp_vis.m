%%%%%%%%
% Author: David Gruskin
% Contact: dcg2153@cumc.columbia.edu
% Last updated: 03/2022
% Project: Brain connectivity at rest predicts individual differences in normative activity during movie watching 
% Description: This code is called by restsync_main.m to create cifti files for visualization in HCP Workbench

% Dependencies: gifti toolbox from https://github.com/gllmflndn/gifti
%%%%%%%%
%% 
function hcp_vis(wb_path,save_stem,write_vec,scalar_template,fname)

% Set paths

suffix = scalar_template(end-11:end);

% Save dscalars
if strcmp(suffix,'.dscalar.nii')
    cii_scalar = ciftiopen(scalar_template, wb_path);
    % map
    for row = 1:length(write_vec)
        cii_scalar.cdata(cii_scalar.cdata == row) = write_vec(row);
    end
    
    save_path = sprintf('%s/%s%s',save_stem,fname,suffix);
    ciftisave(cii_scalar, save_path, wb_path)
end

% Save pscalars
if strcmp(suffix,'.pscalar.nii')
    cii_scalar = ciftiopen(scalar_template, wb_path);
    
    % map
    for row = 1:length(write_vec)
        cii_scalar.cdata(row,:) = write_vec(row);
    end
    
    save_path = sprintf('%s/%s%s',save_stem,fname,suffix);
    ciftisave(cii_scalar, save_path, wb_path)
end
end
