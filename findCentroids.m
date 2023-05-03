% function that taks in 3D data points, cluster labels and the number of clusters
% and outputs the centroids of the clusters
function centroids = getCentroids(data, labels, k)
    centroids = zeros(k, 3);
    for i = 1:k
        cluster = data(labels == i, :);
        centroids(i, :) = mean(cluster);
    end
end
