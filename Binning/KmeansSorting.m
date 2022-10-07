function [index,phase]=KmeansSorting( inhale_idx, exhale_idx , par)



% split inhale and exhale?
if ~par.split_in_exhale
    inhale_idx = [inhale_idx(:);exhale_idx(:)];
    exhale_idx = [];
end


nclusters=par.nclusters;
[cluster_idx_inhale,b_in,ptcd_in,D_inhale]=kmeans(par.surrogate_signal([inhale_idx(:)]),nclusters,'MaxIter',300);


% sort the cluster centroids
[~,s_ind_in]=sort(b_in,'descend');
cluster_ordering = [s_ind_in(:)];


if par.split_in_exhale
    [cluster_idx_exhale,b_ex,ptcd_ex,D_exhale]=kmeans(par.surrogate_signal([exhale_idx(:)]),nclusters,'MaxIter',300);
    [~,s_ind_ex]=sort(b_ex,'ascend');
    cluster_ordering = [s_ind_in(:);s_ind_ex(:)+nclusters];
    nclusters = nclusters*2;
end

k=1;
for i=cluster_ordering.'


    disp(['Reconstructing cluster ',num2str(k),'/',num2str(nclusters)]);


    if k>nclusters/2 && par.split_in_exhale % exhale
        % select N_per_bin spoke sets closest to the cluster centroid
        [~,top_X_points] = sort(D_exhale(:,i-nclusters/2),'ascend');
        N_per_bin = ceil(length(par.surrogate_signal)/nclusters);
        top_X_points=top_X_points(1:N_per_bin);

        respiratory_bin_idx = exhale_idx(top_X_points);

    else %inhale
        % select N_per_bin spoke sets closest to the cluster centroid
        [~,top_X_points] = sort(D_inhale(:,i),'ascend');
        N_per_bin = ceil(length(par.surrogate_signal)/nclusters);
        top_X_points=top_X_points(1:N_per_bin);
        respiratory_bin_idx = inhale_idx(top_X_points);
    end

    index(:,k)=respiratory_bin_idx;


%     respiratory_bin_idx = inhale_exhale_sorted(find(cluster_idx(:)==i));

    cmap = colormap(jet(nclusters));
    cmap = cmap(randperm(length(cmap)),:);

%     if k==1
    figure(99);
%     end
    plot(par.surrogate_signal([respiratory_bin_idx]),'Color',cmap(i,:));
    drawnow;
    hold on;


    k=k+1;
end
phase=[];

end
