%% mvpa, cross decoding, parietooccipital

clear;

dataAll = dir([cd,'\cleaned_data\*.set']);
resultsAll= {};

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
    cfg.resamplefs  = 200;
    CurrData_resampled  = ft_resampledata(cfg,CurrDataTrimmed);

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
    
    
    %% mvpa, attend->unattend
    time = CurrData_A1U1.time{1};
    timepoints = length(time);

    Niter = 50; 
    ntrial_A1U1 = length(CurrData_A1U1.trial);ntrial_A3U1 = length(CurrData_A3U1.trial);ntrial_A1U3 = length(CurrData_A1U3.trial);
    perf_cross_iter = [];
    stat_cross_iter = cell(Niter,1);


    for iter = 1:Niter
        
        % prepare train data and test data
            idx_A1U1 = randperm(ntrial_A1U1);n_div = floor(ntrial_A1U1/1.5);
            idx_A1U1_train = idx_A1U1(1:n_div); idx_A1U1_test = idx_A1U1(n_div+1:end);
            dat_A1U1_train = dat_A1U1(idx_A1U1_train,:,:);dat_A1U1_test = dat_A1U1(idx_A1U1_test,:,:);

            idx_A3U1 = randperm(ntrial_A3U1);
            idx_A3U1_train = idx_A3U1(1:n_div); dat_A3U1_train = dat_A3U1(idx_A3U1_train,:,:);

            idx_A1U3 = randperm(ntrial_A1U3);
            idx_A1U3_test = idx_A1U3(1:ntrial_A1U1 - n_div); dat_A1U3_test = dat_A1U3(idx_A1U3_test,:,:);
            
            dat_train = vertcat(dat_A1U1_train, dat_A3U1_train);
            dat_train_label = [ones(size(dat_A1U1_train,1),1); 2*ones(size(dat_A3U1_train,1),1)];

            dat_test = vertcat(dat_A1U1_test, dat_A1U3_test);
            dat_test_label = [ones(size(dat_A1U1_test,1),1); 2*ones(size(dat_A1U3_test,1),1)];

            cfg_mvpa = [];
            cfg_mvpa.classifier = 'lda';
            cfg_mvpa.metric = 'accuracy';
            cfg_mvpa.cv = 'none';
            
            cfg_mvpa.preprocess = {'zscore' 'average_samples'};
            cfg_mvpa.preprocess_param = {};
            cfg_mvpa.preprocess_param{2} = [];
            cfg_mvpa.preprocess_param{2}.group_size = 5;

            cfg_mvpa.sampledimension = 1;   % trials/samples
            cfg_mvpa.featuredimension = 2;  % channels/features
            
            [~, stat_cross_iter{iter}] = mv_classify_across_time( ...
                cfg_mvpa, ...
                dat_train, dat_train_label, ...
                dat_test, dat_test_label);
            
            perf_cross_iter(iter,:) = stat_cross_iter{iter}.perf(:)';
    end

    stat_mvpa_cross = stat_cross_iter{1};
    stat_mvpa_cross.perf = mean(perf_cross_iter, 1)';
    stat_mvpa_cross.perf_std = std(perf_cross_iter, 0, 1)';
    stat_mvpa_cross.perf_iter = perf_cross_iter;
    stat_mvpa_cross.iter = stat_cross_iter;
    stat_mvpa_cross.time = time;
    
    resultsAll{sub} = stat_mvpa_cross;

end

%% statistics

N = length(resultsAll);
result_average = mv_combine_results(resultsAll, 'average');

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
stat_group = mv_statistics(cfg_stat, resultsAll);


%% visualization

cross_color = [0.8500 0.3250 0.0980];

if iscell(result_average.perf)
    avg = result_average.perf{:};
else
    avg = result_average.perf;
end
if iscell(result_average.perf_std)
    avg_std = result_average.perf_std{:};
else
    avg_std = result_average.perf_std;
end

sig_mask = stat_group.mask;
avg = avg(:);
avg_std = avg_std(:);
sig_mask = sig_mask(:);

if ~exist('time','var') || numel(time) ~= numel(avg)
    time = linspace(-0.2, 2, numel(avg));
end
time = time(:)';

figure('Position',[100 100 700 500])
upperbound = avg + avg_std/sqrt(N);
lowerbound = avg - avg_std/sqrt(N);
fill1 = fill([time,fliplr(time)],[upperbound',fliplr(lowerbound')],cross_color,'FaceAlpha',0.3,'EdgeColor','none');
hold on
plot(time,avg,'Color',cross_color,'LineWidth',2)
hold on
sig_cluster_cross = {};
clusters = bwlabel(sig_mask);
for c = 1:max(clusters)
    cluster{c} = time(find(clusters == c));
    sig_cluster_cross{c} = plot([min(cluster{c}):0.001:max(cluster{c})],repmat(0.48,size([min(cluster{c}):0.001:max(cluster{c})])),'o',...
        'Color',cross_color,'MarkerSize',5,'MarkerFaceColor',cross_color);
    hold on
end

xlim([-0.2 2]);xticks(-0.2:0.2:2);
yline(0.5,'k');ylim([0.44 0.64]);yl = ylim;yticks(0.44:0.02:0.64);
set(gca,'FontSize',13, 'FontName', 'Calibri');
yt = ylabel('Classification Accuracy'); yt.FontSize = 18;
xt = xlabel('Times(s)'); xt.FontSize = 18;
h1 = patch([0 0.2 0.2 0], [yl(1) yl(1) yl(2) yl(2)], [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
h2 = patch([1 1.2 1.2 1], [yl(1) yl(1) yl(2) yl(2)], [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
uistack(h1, 'bottom'); uistack(h2, 'bottom');
ax = gca;ax.TickDir = 'out';ax.LineWidth = 1.2; box off;
if isempty(sig_cluster_cross)
    set([fill1,h1,h2],'HandleVisibility','off')
else
    set([fill1,sig_cluster_cross{:},h1,h2],'HandleVisibility','off')
end
title('Load Classification (At. to Unat.)','FontSize', 20);