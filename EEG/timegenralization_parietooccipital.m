%% mvpa, time generalization, parietooccipital

clear;

dataAll = dir([cd,'\cleaned_data\*.set']);

resultsAll_unattend_load = {};
resultsAll_attend_load = {};

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
    cfg.resamplefs     = 50;
    CurrData_resampled     = ft_resampledata(cfg,CurrDataTrimmed);
    
    % select channels
    cfg = [];
    cfg.channel = {'P1';'P3';'P5';'P7';'P9';'PO7';'PO3';'O1';'Iz';'Oz';'POz';'Pz';'P2';'P4';'P6';'P8';'P10';'PO8';'PO4';'O2'};

    CurrData_cleaned = ft_preprocessing(cfg,CurrData_resampled);
    
    
    % split data in three conditions for load decoding
    
    cfg = [];
    cfg.trials = find(ismember(CurrData_cleaned.trialinfo{:,2},[1,2]));   
    CurrData_A1U1 = ft_selectdata(cfg,CurrData_cleaned);
    dat_A1U1 = permute(cat(3,CurrData_A1U1.trial{:}),[3 1 2]);
    
    cfg = [];
    cfg.trials = find(ismember(CurrData_cleaned.trialinfo{:,2},[3,6]));  
    CurrData_A1U3 = ft_selectdata(cfg,CurrData_cleaned);
    dat_A1U3 = permute(cat(3,CurrData_A1U3.trial{:}),[3 1 2]);
    
    cfg = [];
    cfg.trials = find(ismember(CurrData_cleaned.trialinfo{:,2},[4,5]));   
    CurrData_A3U1 = ft_selectdata(cfg,CurrData_cleaned);
    dat_A3U1 = permute(cat(3,CurrData_A3U1.trial{:}),[3 1 2]);
    
    
    
    %% mvpa
    
    cfg_mvpa = [];
    cfg_mvpa.classifier = 'lda';
    cfg_mvpa.metric = 'accuracy';
    cfg_mvpa.cv = 'kfold';
    cfg_mvpa.k = 5;
    cfg_mvpa.repeat = 5;
    
    cfg_mvpa.preprocess = {'zscore' 'average_samples'};
    cfg_mvpa.preprocess_param = {};
    cfg_mvpa.preprocess_param{2} = [];
    cfg_mvpa.preprocess_param{2}.group_size = 5;

    cfg_mvpa.sampledimension = 1;   % trials/samples
    cfg_mvpa.featuredimension = 2;  % channels/features

    dat_unatload = vertcat(dat_A1U1, dat_A1U3);
    dat_unatload_label = [ones(size(dat_A1U1,1),1); 2*ones(size(dat_A1U3,1),1)];
    [~, stat_mvpa_unattend_load] = mv_classify_timextime(cfg_mvpa, dat_unatload, dat_unatload_label);

    dat_atload = vertcat(dat_A1U1, dat_A3U1);
    dat_atload_label = [ones(size(dat_A1U1,1),1); 2*ones(size(dat_A3U1,1),1)];
    [~, stat_mvpa_attend_load] = mv_classify_timextime(cfg_mvpa, dat_atload, dat_atload_label);
    
    resultsAll_unattend_load{sub} = stat_mvpa_unattend_load;
    resultsAll_attend_load{sub} = stat_mvpa_attend_load;    
    
    
end

%% statistics

% average
N = length(dataAll);
result_average_unattend_load = mv_combine_results(resultsAll_unattend_load, 'average');
result_average_attend_load = mv_combine_results(resultsAll_attend_load, 'average');

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
cfg_stat.clustercritval  = tinv(0.95,N-1); % p<0.05

stat_group_unattend_load = mv_statistics(cfg_stat, resultsAll_unattend_load);
stat_group_attend_load = mv_statistics(cfg_stat, resultsAll_attend_load);


%% visualization


cmap = flipud(brewermap(256,'RdBu'));
ticks = compose('%.1f',0:0.2:2);ticks{1} = '0';
time = [-0.1910:0.02:1.9890];


figure

imagesc(time,time,result_average_unattend_load.perf{:});
axis xy;
caxis([0.44 0.56]);
xlim([-0.2 2]);xticks(0:0.2:2);xticklabels(0:0.2:2);
ylim([-0.2 2]);yticks(-0.2:0.2:2);yticklabels(-0.2:0.2:2);
xline(0,'k');xline(0.2,'--k');xline(1,'--k');xline(1.2,'--k');
yline(0,'k');yline(0.2,'--k');yline(1,'--k');yline(1.2,'--k');
set(gca,'FontSize',12, 'FontName', 'Calibri');
colormap(cmap)
colorbar;
xlabel('Testing time (s)','FontSize',18);ylabel('Training time (s)','FontSize',18)
hold on 
contour(time, time, stat_group_unattend_load.mask, 1, 'LineColor', [0.3 0.3 0.3], 'LineWidth', 1.5);


figure
imagesc(time,time,result_average_attend_load.perf{:});
axis xy;
caxis([0.44 0.56]);
xlim([-0.2 2]);xticks(0:0.2:2);xticklabels(0:0.2:2);
ylim([-0.2 2]);yticks(-0.2:0.2:2);yticklabels(-0.2:0.2:2);
xline(0,'k');xline(0.2,'--k');xline(1,'--k');xline(1.2,'--k');
yline(0,'k');yline(0.2,'--k');yline(1,'--k');yline(1.2,'--k');
set(gca,'FontSize',12, 'FontName', 'Calibri');
colormap(cmap)
colorbar;
xlabel('Testing time (s)','FontSize',18);ylabel('Training time (s)','FontSize',18)
hold on 
contour(time, time, stat_group_attend_load.mask, 1, 'LineColor', [0.3 0.3 0.3], 'LineWidth', 1.5);

