%% mvpa, search across time

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
    
    % A1U1 vs. A1U3
   
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
    [~, stat_mvpa_unattend_load] = mv_classify_across_time(cfg_mvpa, dat_unatload, dat_unatload_label);

    dat_atload = vertcat(dat_A1U1, dat_A3U1);
    dat_atload_label = [ones(size(dat_A1U1,1),1); 2*ones(size(dat_A3U1,1),1)];
    [~, stat_mvpa_attend_load] = mv_classify_across_time(cfg_mvpa, dat_atload, dat_atload_label);

    
    resultsAll_unattend_load{sub} = stat_mvpa_unattend_load;
    resultsAll_attend_load{sub} = stat_mvpa_attend_load;
    
    
end

%% statistics

N = length(dataAll);
% average
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
cfg_stat.clustercritval  = tinv(0.95,N-1); 


stat_group_unattend_load = mv_statistics(cfg_stat, resultsAll_unattend_load);
stat_group_attend_load = mv_statistics(cfg_stat, resultsAll_attend_load);

%% visualization

attend_color = [0 0.4470 0.7410];
unattend_color = [0.9290 0.6940 0.1250];


time = [-0.1985:0.005:1.9965];

figure('Position',[100 100 1000 500])
avg = result_average_attend_load.perf{:};
upperbound = avg + result_average_attend_load.perf_std{:}/sqrt(N);
lowerbound = avg - result_average_attend_load.perf_std{:}/sqrt(N);
fill2 = fill([time,fliplr(time)],[upperbound',fliplr(lowerbound')],attend_color,'FaceAlpha',0.3,'EdgeColor','none');
hold on
plot(time,avg,'Color',attend_color,'LineWidth',2)
hold on
clusters = bwlabel(stat_group_attend_load.mask);
for c = 1:max(clusters)
    cluster{c} = time(find(clusters == c));
    sig_cluster_attend{c} = plot([min(cluster{c}):0.001:max(cluster{c})],repmat(0.47,size([min(cluster{c}):0.001:max(cluster{c})])),'o',...
        'Color',attend_color,'MarkerSize',5,'MarkerFaceColor',attend_color);
    hold on
end


avg = result_average_unattend_load.perf{:};
upperbound = avg + result_average_unattend_load.perf_std{:}/sqrt(N);
lowerbound = avg - result_average_unattend_load.perf_std{:}/sqrt(N);
fill1 = fill([time,fliplr(time)],[upperbound',fliplr(lowerbound')],unattend_color,'FaceAlpha',0.3,'EdgeColor','none');
hold on
plot(time,avg,'Color',unattend_color,'LineWidth',2)
hold on
clusters = bwlabel(stat_group_unattend_load.mask);
for c = 1:max(clusters)
    cluster{c} = time(find(clusters == c));
    sig_cluster_unattend{c} = plot([min(cluster{c}):0.001:max(cluster{c})],repmat(0.48,size([min(cluster{c}):0.001:max(cluster{c})])),'o',...
        'Color',unattend_color,'MarkerSize',5,'MarkerFaceColor',unattend_color);
    hold on
end

xticks(-0.2:0.2:2);xticklabels(-0.2:0.2:2);xlim([-0.2,2]);xl = xlim;
yline(0.5,'k','LineWidth',1.5);ylim([0.44 0.64]);yl = ylim;yticks(0.44:0.02:0.64);
set(gca,'FontSize',18, 'FontName', 'Calibri');
yt = ylabel('Classification Accuracy'); yt.FontSize = 20;
xt = xlabel('Times(s)'); xt.FontSize = 20;
h1 = patch([0 0.2 0.2 0], [yl(1) yl(1) yl(2) yl(2)], [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
h2 = patch([1 1.2 1.2 1], [yl(1) yl(1) yl(2) yl(2)], [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
uistack(h1, 'bottom'); uistack(h2, 'bottom');
ax = gca;ax.TickDir = 'out';ax.LineWidth = 1.2; box off;
                                                                                                                                                                                                                                                                         
set([fill1,fill2,sig_cluster_unattend{:},sig_cluster_attend{:},h1,h2],'HandleVisibility','off')
legend('Attended-load classification (A1U1 vs A3U1)','Unattended-load classification (A1U1 vs A1U3)')
