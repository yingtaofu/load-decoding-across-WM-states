function data_mirrored = mirror_trials_lr_WML(data)

    left_hemifield_chanidx = [1:27];
    right_hemifield_chanidx = [34:36 39:46 49:64];
    
    data_mirrored = data;
    
    % select trials in which the attended state appear on the right
    trial_idx_to_mirror = find(ismember(data.trialinfo{:,2},[2 4 6]));
    
    %mirror
    for trl = trial_idx_to_mirror'
        left_hemifield_data = data_mirrored.trial{trl}(left_hemifield_chanidx,:);
        right_hemifield_data = data_mirrored.trial{trl}(right_hemifield_chanidx,:);
        data_mirrored.trial{trl}(left_hemifield_chanidx,:) = right_hemifield_data;
        data_mirrored.trial{trl}(right_hemifield_chanidx,:) = left_hemifield_data;
    end
   
end
