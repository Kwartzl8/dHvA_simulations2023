function groups = groupEmUp(data, k_dist_perc)
% scale data - equivalent to defining the distance function for DBSCAN
    data_scaled = data;
    data_scaled(:,1) = data(:,1) ./ mean(data(:,1));
    data_scaled(:,2) = data(:,2) ./ mean(data(:,2));
    
    minpts = 3;
    kD = pdist2(data_scaled, data_scaled, 'euc', 'Smallest', minpts);
    sorted_kD = sort(kD(end,:));

    epsilon = sorted_kD(floor(length(sorted_kD) .* k_dist_perc));
    labels = dbscan(data_scaled, epsilon, minpts);
    
    group_no = length(unique(labels));
    % check if the grouping has outliers and remove them from the groups
    if any(labels==-1)
        group_no = group_no - 1;
    end
    groups = cell(1, group_no);
    for j = 1:(group_no)
        groups{j} = find(labels==j);
    end
end