function [significant_point] = cluster_based_permutation_test(data,clustercritval,alpha,n_permutation,nullvalue,tail)

    % data: freq x timepoint x sample matrix
    % clustercritval: the significance threshold (p value) for identifying
    % clusters, one-sided
    % alpha: the significance threshold (p value) for identifying
    % significant clusters (one-sided)
    % n_permutation: permutation numbers
    % nullvalue: null value for one-sample t-test(only one value)
    % tail: -1 or 1,  specifies whether the lower tail (-1), or the upper tail (+1) is computed Typically,
    

    n_freq = size(data,1);
    n_points = size(data,2);
    n_subjects = size(data,3);

    t_values = zeros(n_freq,n_points);
    
    significant_point = zeros(n_freq,n_points);
    
    if tail == 1
        cluster_threshold = tinv(1-clustercritval,n_subjects-1);
    else
        cluster_threshold = tinv(clustercritval,n_subjects-1);
    end

    for fi = 1:n_freq
        for t = 1:n_points
            [~,~,~,stats] = ttest(squeeze(data(fi,t,:)),nullvalue);
            t_values(fi,t) = stats.tstat;
        end
    end
    
    if tail == 1
        significant_mask = t_values > cluster_threshold;% one sided
    else
        significant_mask = t_values < cluster_threshold;% one sided
    end

    clusters = bwlabel(significant_mask,4);
    
    if max(clusters,[],'all') ~=0

        for c = 1:max(clusters,[],'all')
            real_clusters_stats(c) = sum(t_values(clusters==c));
        end

        % randomly shuffled and test
        
        max_perm_clusters_stats = zeros(1, n_permutation);

        parfor perm = 1:n_permutation
            
            disp(perm)

            % randomly shuffle the sign
            flipsign = 2*(rand(n_subjects,1)>0.5)-1;
            flipsign = reshape(flipsign,1,1,n_subjects);
            perm_data = (data-nullvalue).*flipsign;

            % t test
            perm_t_values = zeros(n_freq,n_points);

            for fi = 1:n_freq
                for t = 1:n_points
                    [~,~,~,stats] = ttest(squeeze(perm_data(fi,t,:)),0);
                    perm_t_values(fi,t) = stats.tstat;
                end
            end

            if tail == 1
                perm_significant_mask = perm_t_values > cluster_threshold;
            else
                perm_significant_mask = perm_t_values < cluster_threshold;
            end
            
            perm_clusters = bwlabel(perm_significant_mask,4);
            
            perm_clusters_stats = [];
            
            if max(perm_clusters,[],'all')~=0
                for c = 1:max(perm_clusters,[],'all')
                    perm_clusters_stats(c) = sum(perm_t_values(perm_clusters==c));
                end
                if tail == 1
                    max_perm_clusters_stats(perm) = max(perm_clusters_stats);
                else
                    max_perm_clusters_stats(perm) = min(perm_clusters_stats);
                end
            else
                max_perm_clusters_stats(perm) = 0;
            end

        end
        
        if tail == 1
            threshold = prctile(max_perm_clusters_stats,(1-alpha)*100);
            for c = 1:max(clusters,[],'all')
                if real_clusters_stats(c) > threshold
                    significant_point(clusters == c) = 1;
                end
            end
        else
            threshold = prctile(max_perm_clusters_stats,alpha*100);
            for c = 1:max(clusters,[],"all")
                if real_clusters_stats(c) < threshold
                    significant_point(clusters == c) = 1;
                end
            end
        end
    end
    
end