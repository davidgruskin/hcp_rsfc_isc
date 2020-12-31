%%%%%%%%
% Author: David Gruskin
% Contact: dcg2153@cumc.columbia.edu
% Last updated: 12/2020
% Project: Brain connectivity at rest predicts individual differences in normative activity during movie watching 
% Description: This script outputs a list of volumes that will be censored based on framewise displacement. 

%%%%%%%%

function [vols_to_censor] = fd_censoring(fd_trace,fd_threshold,num_contiguous)
 
fd_trace(fd_trace>fd_threshold,1) = NaN;
for vol = 1:length(fd_trace)
    if isnan(fd_trace(vol))
        if vol + num_contiguous > length(fd_trace)
            fd_trace(vol+1:end,1) = NaN;
        else
        if ~isnan(fd_trace(vol+1))
            for next_vol = 1:num_contiguous - 1
                if ~isnan(fd_trace(vol + 1 + next_vol))
                    island_status = 0;
                else
                    island_status = 1;
                    break
                end
            end
            if island_status == 1
                fd_trace(vol+1:vol+1+next_vol,1) = NaN;
            end
        end
        end
    end
end

vols_to_censor = find(isnan(fd_trace));
vols_to_censor = vols_to_censor + 1;
end