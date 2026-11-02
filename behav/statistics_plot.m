clear; clc;
data_all = readtable('DataAll.xlsx');

%% data reorganize
data_all.load = strcat("A",string(data_all.AtLoad), "U", string(data_all.UnatLoad));
data_all.Probe1ACC = double(data_all.Probe1ACC);
data_all.Probe2ACC = double(data_all.Probe2ACC);

summary_attend = groupsummary( ...
    data_all, ...
    {'Name','load'}, ...
    'mean','Probe1ACC');

summary_unattend = groupsummary( ...
    data_all, ...
    {'Name','load'}, ...
    'mean','Probe2ACC');

S = stack( ...
    [summary_attend, summary_unattend(:,{'mean_Probe2ACC'})], ...
    {'mean_Probe1ACC','mean_Probe2ACC'}, ...
    'NewDataVariableName','acc', ...
    'IndexVariableName','state');

summary_data = table;
summary_data.Id = repelem([1:1:33]',3*2);
summary_data.load = S.load;
summary_data.state = repmat({'attend','unattend'},1,33*3)';
summary_data.ACC = S.acc;


%% statistics

% data reorg

summary_data.load  = categorical(summary_data.load);
summary_data.state = categorical(summary_data.state);
loads = categories(summary_data.load);
states = categories(summary_data.state);
subs = unique(summary_data.Id);
tab_data = nan(numel(subs), numel(loads)*numel(states));

colNames = {};
c = 1;

for s = 1:numel(states)
    for l = 1:numel(loads)
        colNames{c} = sprintf('%s_%s', states{s}, loads{l});
        c = c + 1;
    end
end

for i = 1:numel(subs)
    for s = 1:numel(states)
        for l = 1:numel(loads)
            idx = summary_data.Id == subs(i) & ...
                  summary_data.state == states{s} & ...
                  summary_data.load == loads{l};
            tab_data(i, (s-1)*numel(loads)+l) = summary_data.ACC(idx);
        end
    end
end

tab_data = array2table(tab_data, 'VariableNames', colNames);
tab_data.subs = subs;
tab_data = movevars(tab_data,'subs','Before',1);

% r-anova
WithinDesign = table(repelem(states, numel(loads)),repmat(loads', 1, numel(states))','VariableNames', {'state','load'});
rm = fitrm(tab_data, 'attend_A1U1-unattend_A3U1 ~ 1', 'WithinDesign',WithinDesign);
ranovatbl = ranova(rm, 'WithinModel', 'state*load');
disp(ranovatbl)
F = ranovatbl.F;

effects = {'state','load','state:load'};
eta_p2  = struct();

% state
SS_state = ranovatbl.SumSq( ...
    strcmp(ranovatbl.Properties.RowNames,'(Intercept):state'));
SS_state_err = ranovatbl.SumSq( ...
    strcmp(ranovatbl.Properties.RowNames,'Error(state)'));
eta_p2.state = SS_state / (SS_state + SS_state_err);

% load
SS_load = ranovatbl.SumSq( ...
    strcmp(ranovatbl.Properties.RowNames,'(Intercept):load'));
SS_load_err = ranovatbl.SumSq( ...
    strcmp(ranovatbl.Properties.RowNames,'Error(load)'));
eta_p2.load = SS_load / (SS_load + SS_load_err);

% interaction
SS_int = ranovatbl.SumSq( ...
    strcmp(ranovatbl.Properties.RowNames,'(Intercept):state:load'));
SS_int_err = ranovatbl.SumSq( ...
    strcmp(ranovatbl.Properties.RowNames,'Error(state:load)'));
eta_p2.state_load = SS_int / (SS_int + SS_int_err);

% effect size
disp(eta_p2)

% simple comparisons

[~, p11, ~, stats11] = ttest(tab_data.attend_A1U1, tab_data.unattend_A1U1);
dz11 = stats11.tstat / sqrt(stats11.df + 1);

[~, p13, ~, stats13] = ttest(tab_data.attend_A1U3, tab_data.unattend_A1U3);
dz13 = stats13.tstat / sqrt(stats13.df + 1);

[~, p31, ~, stats31] = ttest(tab_data.attend_A3U1, tab_data.unattend_A3U1);
dz31 = stats31.tstat / sqrt(stats31.df + 1);

%% plot

figure; hold on;

colors_a = [0 112 192;
          112 173 71;
          255 192 0]/255;

colors_u = [143 207 255;
          183,216,160;
          255 219 105]/255;

offset = 0.3; 
violin_jitter = 0.08;
scale = 2.4;


for i = 1:numel(loads)

    for k = 1:numel(subs)

        sid = subs(k);

        % attend
        idx_a = summary_data.Id == sid & ...
                summary_data.load == loads{i} & ...
                summary_data.state == "attend";

        % unattend
        idx_u = summary_data.Id == sid & ...
                summary_data.load == loads{i} & ...
                summary_data.state == "unattend";

        y_a = summary_data.ACC(idx_a);
        y_u = summary_data.ACC(idx_u);

        x_a = i*scale - offset + violin_jitter;
        x_u = i*scale + offset - violin_jitter;

        % --------
        % lines
        % --------
        plot([x_a, x_u], [y_a, y_u], ...
            '-', ...
            'Color', [0.9 0.9 0.9], ...
            'LineWidth', 0.2);

        % --------
        % individual data point
        % --------
        scatter(x_a, y_a, ...
            30, colors_a(i,:), ...
            'filled', ...
            'MarkerFaceAlpha', 0.8,...
            'MarkerEdgeColor', [1 1 1],...
            'MarkerEdgeAlpha', 1);

        scatter(x_u, y_u, ...
            30, colors_u(i,:), ...
            'filled', ...
            'MarkerFaceAlpha', 0.8, ...
            'MarkerEdgeColor', [1 1 1],...
            'MarkerEdgeAlpha', 1);
    end

    for s = 1:numel(states)

        % get data
        idx = summary_data.load  == loads{i} & ...
              summary_data.state == states{s};

        data = summary_data.ACC(idx);

        if states{s} == "attend"
            x0 = i*scale - offset;
            side = -1;
            c = colors_a(i,:);
            x_mean = x0 + violin_jitter;
            x_a = x0 + violin_jitter; 
            y_a = mean(data);
        else
            x0 = i*scale + offset;
            side = +1;
            c = colors_u(i,:);
            x_mean = x0 - violin_jitter;
            x_u = x0 - violin_jitter;
            y_u = mean(data);
        end

        [f, y] = ksdensity(data);
        f = f ./ 18;
        
        x_violin = x0 + side * f;
        
        X = [x0*ones(size(y)), x_violin];
        Y = [fliplr(y), y];
        
        fill(X, Y, c, 'FaceAlpha', 0.8, 'EdgeColor', 'none');

        % error bar
        sem = std(data) / sqrt(numel(data));
        errorbar(x_mean, mean(data), sem, 'k', 'LineWidth', 1, 'CapSize', 0);
        
        % mean marker
        plot(x_mean, mean(data), 'k.', 'MarkerSize',8);
    end
    
    plot([x_a x_u], [y_a, y_u], '-','Color', [0.3 0.3 0.3], 'LineWidth', 0.8);
    
    % significance
    if i == 1
        y_temp = 0.73; y_temp_length = y_temp + 0.01;

    else
        y_temp = 0.45; y_temp_length = y_temp + 0.01;
    end
    plot ([x_a x_u], [y_temp y_temp], 'k-','LineWidth', 1);
    plot ([x_a x_a], [y_temp y_temp_length], 'k-','LineWidth', 1);
    plot ([x_u x_u], [y_temp y_temp_length], 'k-','LineWidth', 1);
    text(mean([x_a x_u]), y_temp - 0.04, '***', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 11);

end

% legend
x0 = 0.18;
y0 = 0.18;

boxW = 0.025;   
boxH = 0.025;   
gapX = 0.004;   
gapY = 0.04;    

for i = 1:3
    annotation('rectangle', [x0 + (i-1)*(boxW+gapX), y0+gapY, boxW, boxH], 'FaceColor', colors_a(i,:),'FaceAlpha',0.8, 'EdgeColor', 'none');
end

annotation('textbox', [x0 + 3*(boxW+gapX) , y0+gapY, 0.2, boxH],'String', 'Attended', 'LineStyle', 'none', ...
    'FontName', 'Calibri', 'FontSize', 12, 'VerticalAlignment', 'middle');

for i = 1:3
    annotation('rectangle', [x0 + (i-1)*(boxW+gapX), y0, boxW, boxH], 'FaceColor', colors_u(i,:),'FaceAlpha',0.8,  'EdgeColor', 'none');
end

annotation('textbox', [x0 + 3*(boxW+gapX), y0, 0.2, boxH], 'String', 'Unattended', 'LineStyle', 'none', ...
    'FontName', 'Calibri', 'FontSize', 12,  'VerticalAlignment', 'middle');



xlim([0.6 numel(loads)*scale+1.5]);ylim([0.3 1.1])
xticks(scale:scale:numel(loads)*scale);yticks(0.3:0.1:1)
xticklabels({'A1U1','A1U3','A3U1'}); yticklabels(0.3:0.1:1)
set(gca,'FontSize',14, 'FontName', 'Calibri')
ax = gca;ax.TickDir = 'out';ax.LineWidth = 1; 
yl = ylabel('Accuracy');yl.FontSize = 18;
xl = ax.XLim;
yl = ax.YLim;
hold on
rectangle('Position', [xl(1), yl(1), diff(xl), diff(yl)], ...
          'EdgeColor', 'k', ...
          'LineWidth', 1, ...
          'FaceColor', 'none');

