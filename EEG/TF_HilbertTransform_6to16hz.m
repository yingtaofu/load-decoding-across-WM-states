%% Time frequency analysis, Hilbert transform, 6 - 16 Hz
clear; clc;


dataAll = dir([cd,'\cleaned_data\*.set']);
save_path = [cd,'\6-16hz_frequency_result_Hilbert\data\'];

if ~exist(save_path)
    mkdir(save_path);
end

%% TF
% parameters
min_freq =  6;
max_freq = 16;
frex = 6:0.5:16;
num_freq = length(frex);
freq_bandwidth = 1;


chan_ipsi = {'P1';'P3';'P5';'P7';'P9';'PO7';'PO3';'O1'};
chan_contra = {'P2';'P4';'P6';'P8';'P10';'PO8';'PO4';'O2'};

for sub = 1:length(dataAll)

    eegpower = {}; 

    for fi = 1:num_freq
        
        temp_band = [frex(fi)-freq_bandwidth/2 frex(fi)+freq_bandwidth/2];
        hilbert_tempdata = {}; tempdata = {};
    
        % import data
        EEG = pop_loadset('filename',dataAll(sub).name,'filepath',dataAll(sub).folder);
        baseidx = dsearchn(EEG.times',[-500 -200]');

        % filter
        EEG = pop_eegfiltnew(EEG, 'locutoff',temp_band(1));
        EEG = pop_eegfiltnew(EEG, 'hicutoff',temp_band(2));

        CurrData = eeglab2fieldtrip(EEG, 'raw');
        CurrData_mirrored = mirror_trials_lr_WML(CurrData);
        [~,chan_ipsi_idx] = ismember(chan_ipsi ,CurrData_mirrored.label);
        [~,chan_contra_idx] = ismember(chan_contra ,CurrData_mirrored.label);
    
        % split data
        cfg = [];
        cfg.trials = find(ismember(CurrData_mirrored.trialinfo{:,2},[1,2]));   
        CurrData_A1U1 = ft_selectdata(cfg,CurrData_mirrored);
        tempdata{1} = permute(cat(3,CurrData_A1U1.trial{:}),[2,3,1]);
        
        cfg = [];
        cfg.trials = find(ismember(CurrData_mirrored.trialinfo{:,2},[3,6]));  
        CurrData_A1U3 = ft_selectdata(cfg,CurrData_mirrored);
        tempdata{2}= permute(cat(3,CurrData_A1U3.trial{:}),[2,3,1]);
        
        cfg = [];
        cfg.trials = find(ismember(CurrData_mirrored.trialinfo{:,2},[4,5]));   
        CurrData_A3U1 = ft_selectdata(cfg,CurrData_mirrored);
        tempdata{3} = permute(cat(3,CurrData_A3U1.trial{:}),[2,3,1]);

        for con = 1:3
            % apply hilbert tranform
            hilbert_tempdata{con} = hilbert(tempdata{con});
            % average across trials
            temppower = squeeze(mean(abs(hilbert_tempdata{con}).^2,2))';
            eegpower{con}(fi,:,:) = 10*log10(temppower./mean(temppower(:, baseidx(1):baseidx(2)), 2));
        end
    end
    eegpower_A1U1_ipsi = eegpower{1}(:,chan_ipsi_idx,:);eegpower_A1U1_contra = eegpower{1}(:,chan_contra_idx,:);
    eegpower_A1U3_ipsi = eegpower{2}(:,chan_ipsi_idx,:);eegpower_A1U3_contra = eegpower{2}(:,chan_contra_idx,:);
    eegpower_A3U1_ipsi = eegpower{3}(:,chan_ipsi_idx,:);eegpower_A3U1_contra = eegpower{3}(:,chan_contra_idx,:);
    
    save([save_path, dataAll(sub).name(1:end-4),'.mat'],"eegpower_A1U1_ipsi","eegpower_A1U1_contra",...
        "eegpower_A1U3_ipsi","eegpower_A1U3_contra","eegpower_A3U1_ipsi","eegpower_A3U1_contra");
end

%% average across subjects, calculate lateralized power

subAll = dir([cd,'\6-16hz_frequency_result_Hilbert\data\*.mat']);

% freq x time x sub
A1U1_ipsi = []; A1U1_contra = []; A1U1_lateralized = [];
A1U3_ipsi = []; A1U3_contra = []; A1U3_lateralized = [];
A3U1_ipsi = []; A3U1_contra = []; A3U1_lateralized = [];

for sub = 1:length(subAll)

    load([subAll(sub).folder,'\',subAll(sub).name])

    A1U1_ipsi(:,:,sub) = mean(eegpower_A1U1_ipsi,2);
    A1U1_contra(:,:,sub) = mean(eegpower_A1U1_contra,2);
    A1U1_lateralized(:,:,sub) = mean(eegpower_A1U1_contra - eegpower_A1U1_ipsi,2);

    A1U3_ipsi(:,:,sub) = mean(eegpower_A1U3_ipsi,2);
    A1U3_contra(:,:,sub) = mean(eegpower_A1U3_contra,2);
    A1U3_lateralized(:,:,sub) = mean(eegpower_A1U3_contra - eegpower_A1U3_ipsi,2);

    A3U1_ipsi(:,:,sub) = mean(eegpower_A3U1_ipsi,2);
    A3U1_contra(:,:,sub) = mean(eegpower_A3U1_contra,2);
    A3U1_lateralized(:,:,sub) = mean(eegpower_A3U1_contra - eegpower_A3U1_ipsi,2);
    
end

avg_A1U1_ipsi = mean(A1U1_ipsi,3); avg_A1U1_contra = mean(A1U1_contra,3); avg_A1U1_lateralized = mean(A1U1_lateralized,3);

avg_A1U3_ipsi = mean(A1U3_ipsi,3); avg_A1U3_contra = mean(A1U3_contra,3); avg_A1U3_lateralized = mean(A1U3_lateralized,3);

avg_A3U1_ipsi = mean(A3U1_ipsi,3); avg_A3U1_contra = mean(A3U1_contra,3); avg_A3U1_lateralized = mean(A3U1_lateralized,3);



%% statistics

% cluster permutation
stat_group_A1U1 = [];
stat_group_A1U3 = [];
stat_group_A3U1 = [];


clustercritval = 0.05;
alpha = 0.05;
n_permutation = 1000;
nullvalue = 0;
tail = -1; % whether lower than nullvalue

stat_group_A1U1.mask1 = cluster_based_permutation_test(A1U1_lateralized,clustercritval,alpha,n_permutation,nullvalue,tail);
stat_group_A1U3.mask1 = cluster_based_permutation_test(A1U3_lateralized,clustercritval,alpha,n_permutation,nullvalue,tail);
stat_group_A3U1.mask1 = cluster_based_permutation_test(A3U1_lateralized,clustercritval,alpha,n_permutation,nullvalue,tail);

tail = 1;
stat_group_A1U1.mask2 = cluster_based_permutation_test(A1U1_lateralized,clustercritval,alpha,n_permutation,nullvalue,tail);
stat_group_A1U3.mask2 = cluster_based_permutation_test(A1U3_lateralized,clustercritval,alpha,n_permutation,nullvalue,tail);
stat_group_A3U1.mask2 = cluster_based_permutation_test(A3U1_lateralized,clustercritval,alpha,n_permutation,nullvalue,tail);

stat_group_A1U1.mask = stat_group_A1U1.mask1 + stat_group_A1U1.mask2;
stat_group_A1U3.mask = stat_group_A1U3.mask1 + stat_group_A1U3.mask2;
stat_group_A3U1.mask = stat_group_A3U1.mask1 + stat_group_A3U1.mask2;
%%

cmap = flipud(brewermap(256,'RdBu'));
ticks = compose('%.1f',0:0.2:2);ticks{1} = '0';

figure('Position',[100 100 900 250])
contourf(CurrData_mirrored.time{1},frex,avg_A1U1_lateralized,30,'linecolor','none')
set(gca,'clim',[-0.4 0.4],'xlim',[-0.2 2],'xtick',-0.2:0.2:2,'ytick',6:2:16)
title('A1U1 lateralized')
xlabel('Times(s)');ylabel('Frequency(Hz)');
set(gca, 'TickLength',[0.005 0.005], 'LineWidth',0.8 , 'FontSize',13, 'FontName', 'Calibri')
colormap(cmap); cb = colorbar; cb.Label.String = 'dB power'; cb.Label.FontSize = 13;cb.Label.FontName = 'Calibri';
hold on
contour(CurrData_mirrored.time{1}, frex, stat_group_A1U1.mask, 1, 'LineColor', [0.3 0.3 0.3], 'LineWidth', 0.8);
fig = gcf;
set(fig,'Renderer','painters')
print(fig,'6-16hz_frequency_result_Hilbert\A1U1lateralized_6-16hz.svg','-dsvg','-vector')

figure('Position',[100 100 900 250])
contourf(CurrData_mirrored.time{1},frex,avg_A1U3_lateralized,30,'linecolor','none')
set(gca,'clim',[-0.4 0.4],'xlim',[-0.2 2],'xtick',-0.2:0.2:2,'ytick',6:2:16)
title('A1U3 lateralized')
xlabel('Times(s)');ylabel('Frequency(Hz)');
set(gca, 'TickLength',[0.005 0.005], 'LineWidth',0.8, 'FontSize',13, 'FontName', 'Calibri')
colormap(cmap); cb = colorbar; cb.Label.String = 'dB power'; cb.Label.FontSize = 13;cb.Label.FontName = 'Calibri';
hold on
contour(CurrData_mirrored.time{1}, frex, stat_group_A1U3.mask, 1, 'LineColor', [0.3 0.3 0.3], 'LineWidth', 0.8);
fig = gcf;
set(fig,'Renderer','painters')
print(fig,'6-16hz_frequency_result_Hilbert\A1U3lateralized_6-16hz.svg','-dsvg','-vector')

figure('Position',[100 100 900 250])
contourf(CurrData_mirrored.time{1},frex,avg_A3U1_lateralized,30,'linecolor','none')
set(gca,'clim',[-0.4 0.4],'xlim',[-0.2 2],'xtick',-0.2:0.2:2,'ytick',6:2:16)
title('A3U1 lateralized')
xlabel('Times(s)');ylabel('Frequency(Hz)');
set(gca, 'TickLength',[0.005 0.005], 'LineWidth',0.8, 'FontSize',13, 'FontName', 'Calibri')
colormap(cmap); cb = colorbar; cb.Label.String = 'dB power'; cb.Label.FontSize = 13;cb.Label.FontName = 'Calibri';
hold on
contour(CurrData_mirrored.time{1}, frex, stat_group_A3U1.mask, 1, 'LineColor', [0.3 0.3 0.3], 'LineWidth', 0.8);
fig = gcf;
set(fig,'Renderer','painters')
print(fig,'6-16hz_frequency_result_Hilbert\A3U1lateralized_6-16hz.svg','-dsvg','-vector')