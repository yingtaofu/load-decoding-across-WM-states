%% searchlight anlysis, remap for different memory configurations

% decoding load: remap for different positions of states (cue left or cue right)


% to see the distribution of channels
clear;

dataAll = dir([cd,'\cleaned_data\*.set']);
save_path = [cd,'\searchlight_chantime_remap\data\'];

if ~exist(save_path)
    mkdir(save_path);
end

%%
for sub = 1:length(dataAll)
        
    % import data
    EEG = pop_loadset('filename',dataAll(sub).name,'filepath',dataAll(sub).folder);
    CurrData = eeglab2fieldtrip(EEG, 'raw');
    
    % latency window
    cfg = [];
    cfg.latency = [-0.2 2];
    CurrDataTrimmed = ft_selectdata(cfg, CurrData);
    
    % resample
    cfg = [];
    cfg.resamplefs     = 200;
    CurrData_resampled     = ft_resampledata(cfg,CurrDataTrimmed);
    
    % remove useless channels
    cfg = [];
    if isequal(dataAll(sub).name,'12bxx.set')% sub 12: the right HEOG channel was excluded previously due to missing data.
        cfg.channel = {'all','-HEOG-L','-VEOG-U','-VEOG-D'};
    elseif isequal(dataAll(sub).name,'28SF.set')% sub 28: HEOG and VEOG was excluded previously manually
        cfg.channel = {'all'};
    else
        cfg.channel = {'all','-HEOG-L','-HEOG-R','-VEOG-U','-VEOG-D'};
    end

    CurrData_cleaned = ft_preprocessing(cfg,CurrData_resampled);
    
    % mirror trials for different configurations
    
    CurrData_mirrored_4loaddecoding = mirror_trials_lr_WML(CurrData_cleaned);
    
    % split data in three conditions for load decoding
    
    cfg = [];
    cfg.trials = find(ismember(CurrData_mirrored_4loaddecoding.trialinfo{:,2},[1,2]));   
    CurrData_A1U1 = ft_selectdata(cfg,CurrData_mirrored_4loaddecoding);
    n_A1U1 = numel(CurrData_A1U1.trial);
    
    cfg = [];
    cfg.trials = find(ismember(CurrData_mirrored_4loaddecoding.trialinfo{:,2},[3,6]));   
    CurrData_A1U3 = ft_selectdata(cfg,CurrData_mirrored_4loaddecoding);
    n_A1U3 = numel(CurrData_A1U3.trial);
    
    cfg = [];
    cfg.trials = find(ismember(CurrData_mirrored_4loaddecoding.trialinfo{:,2},[4,5]));   
    CurrData_A3U1 = ft_selectdata(cfg,CurrData_mirrored_4loaddecoding);
    n_A3U1 = numel(CurrData_A3U1.trial);
    
    
    %% search across time and channels
    
    latency_window = [-0.2:0.2:2];
    
    time = CurrData_cleaned.time{1};
    time_win = [];
    
    % prepare neighbourhood structure
    
    cfg = [];
    cfg.method = 'triangulation';
    cfg.layout = 'biosemi64.lay';
    cfg.channel = CurrData_cleaned.label;
    neighbours = ft_prepare_neighbours(cfg);
    
    for t = 1:length(latency_window)-1
        
        % exact time window
        time_win(1) = min(time(time >= latency_window(t)));
        time_win(2) = max(time(time <= latency_window(t+1)));

        % A1U1 vs. A1U3

        cfg = [];
        cfg.method = 'mvpa';
        cfg.latency= time_win;
        cfg.avgovertime = 'no';
        cfg.features = 'time';
        cfg.neighbours = neighbours;
        cfg.mvpa.preprocess = {'zscore' 'average_samples'};
        cfg.mvpa.preprocess_param{2} = {'group_size' 5};
        cfg.mvpa.classifier = 'lda';
        cfg.mvpa.metric = 'accuracy';
        cfg.mvpa.k = 5;
        cfg.mvpa.repeat = 5;
        cfg.design = [ones(n_A1U1,1);2*ones(n_A1U3,1)];

        stat_mvpa_unattend_load{t} = ft_timelockstatistics(cfg,CurrData_A1U1,CurrData_A1U3);

        % A1U1 vs. A3U1

        cfg.design = [ones(n_A1U1,1);2*ones(n_A3U1,1)];
        stat_mvpa_attend_load{t} = ft_timelockstatistics(cfg,CurrData_A1U1,CurrData_A3U1);

    end
    
    save([save_path, dataAll(sub).name(1:end-4),'.mat'],'stat_mvpa_unattend_load','stat_mvpa_attend_load')

end


%% statistics

subAll = dir([save_path,'\*.mat']);
N = length(subAll);
latency_window = [-0.2:0.2:2];
n_channels = 64;

result_avgsub_u = {};
result_avgsub_a = {};

results_subAll_u = {};
results_subAll_a = {};

stat_group_u = {};
stat_group_a = {};

% reorg data
for sub = 1:N % for every sub
    
    load([subAll(sub).folder,'\',subAll(sub).name])
    
    for t = 1:length(latency_window)-1 % for every time window
        results_subAll_u{t}{sub}= stat_mvpa_unattend_load{t};
        results_subAll_u{t}{sub}.mvpa.perf = results_subAll_u{t}{sub}.mvpa.perf{1};
        results_subAll_u{t}{sub}.mvpa.perf_std = results_subAll_u{t}{sub}.mvpa.perf_std{1};
        results_subAll_a{t}{sub}= stat_mvpa_attend_load{t}; 
        results_subAll_a{t}{sub}.mvpa.perf = results_subAll_a{t}{sub}.mvpa.perf{1};
        results_subAll_a{t}{sub}.mvpa.perf_std = results_subAll_a{t}{sub}.mvpa.perf_std{1};
    end

end


% prepare neigbours
num_channels = length(neighbours); % 64 channels
channel_labels = {neighbours.label}; 
binary_neighbors = cell(1, num_channels);

for i = 1:num_channels
    
    binary_vector = zeros(1, num_channels);
    
    current_neighbors = neighbours(i).neighblabel;
    
    for j = 1:length(current_neighbors)
        idx = find(strcmp(channel_labels, current_neighbors{j}));
        if ~isempty(idx)
            binary_vector(idx) = 1;
        end
    end

    binary_neighbors{i} = binary_vector;
end


for t = 1:length(latency_window)-1 

    % average across subjects
    cfg = [];
    cfg.parameter = 'avg';
    for sub_idx = 1:numel(results_subAll_u{t})
        results_subAll_u{t}{sub_idx}.avg = results_subAll_u{t}{sub_idx}.accuracy(:);
        results_subAll_u{t}{sub_idx}.dimord = 'chan_time';
        results_subAll_a{t}{sub_idx}.avg = results_subAll_a{t}{sub_idx}.accuracy(:);
        results_subAll_a{t}{sub_idx}.dimord = 'chan_time';
    end
    result_avgsub_u{t} = ft_timelockgrandaverage(cfg,results_subAll_u{t}{:});
    result_avgsub_a{t} = ft_timelockgrandaverage(cfg,results_subAll_a{t}{:});

    % significant statistics on every time window
    cfg_stat = [];
    cfg_stat.metric          = 'accuracy';
    cfg_stat.test            = 'permutation';
    cfg_stat.correctm        = 'cluster';  % correction method is cluster
    cfg_stat.n_permutations  = 1000;
    cfg_stat.clusterstatistic = 'maxsum';
    cfg_stat.alpha           = 0.05; % use standard significance threshold of 5%
    cfg_stat.design          = 'within';
    cfg_stat.statistic       = 'ttest';
    cfg_stat.null            = 0.5;
    cfg_stat.clustercritval  = tinv(0.95,N-1); 
    cfg_stat.neighbours      = binary_neighbors;
    
    try
        stat_group_u{t} = mv_statistics(cfg_stat,cellfun(@(x) x.mvpa, results_subAll_u{t}, 'UniformOutput', false));
    catch
        stat_group_u{t}.mask = zeros(1,n_channels);
    end
    try
        stat_group_a{t} = mv_statistics(cfg_stat,cellfun(@(x) x.mvpa, results_subAll_a{t}, 'UniformOutput', false) );
    catch
        stat_group_a{t}.mask = zeros(1,n_channels);
    end

        
end

%% visualization

load biosemi64_layout.mat
cmap = jet;

picsave_path = [cd,'\searchlight_chantime_remap\results\'];

if ~exist(picsave_path)
    mkdir(picsave_path);
end

for t = 1:length(latency_window)-1
    
    cfg = [];
    cfg.layout = layout;
    cfg.parameter = 'avg';
    cfg.colorbar = 'yes';
    cfg.style = 'straight';
    cfg.zlim = [0.43 0.57];
    cfg.marker = 'off';
    cfg.highlight = 'on';
    cfg.highlightsymbol = 'o';
    cfg.highlightcolor     = 'k';
    cfg.comment = 'no';
    
    if ~nnz(stat_group_u{t}.mask')==0
        cfg.highlightchannel = result_avgsub_u{t}.label(stat_group_u{t}.mask');
    else
        cfg.highlight = 'off';
    end
    ft_topoplotER(cfg,result_avgsub_u{t});
    colormap(cmap);    
    
    if t == length(latency_window)-1
        colorbar;
        cb = colorbar;
        cb.Ticks = [0.43 0.5 0.57];
        cb.FontName = 'Calibri';
        cb.FontSize = 20;
        cb.Label.String = 'Classification Accuracy';
        cb.Label.Rotation = 270;
        cb.Label.Position = cb.Label.Position + [1 0 0];
    else
        colorbar off
    end
    saveas(gcf, [picsave_path,'\u',num2str(latency_window(t)),'s - ',num2str(latency_window(t+1)),'s','.svg']);
    close Figure 1
  
    if ~nnz(stat_group_a{t}.mask')==0
        cfg.highlight = 'on';
        cfg.highlightchannel = result_avgsub_a{t}.label(stat_group_a{t}.mask');
    else
        cfg.highlight = 'off';
    end
    ft_topoplotER(cfg,result_avgsub_a{t})
    colormap(cmap);   
    cb = colorbar;
    cb.Ticks = [0.43:0.02:0.57];
    if t == length(latency_window)-1
        colorbar;
        cb = colorbar;
        cb.Ticks = [0.43 0.5 0.57];
        cb.FontName = 'Calibri';
        cb.FontSize = 20;
        cb.Label.String = 'Classification Accuracy';
        cb.Label.Rotation = 270;
        cb.Label.Position = cb.Label.Position + [1 0 0];
    else
        colorbar off
    end
    saveas(gcf, [picsave_path,'\a',num2str(latency_window(t)),'s - ',num2str(latency_window(t+1)),'s','.svg']);
    close Figure 1

end
