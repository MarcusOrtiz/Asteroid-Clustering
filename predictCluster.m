% Function that takes in 3D centroids and a data sample and returns the
% closest centroid to each data point

function [closestCentroid] = findClosestCentroids(centroids, X)
    % Iterate through each data point
    for i = 1:size(X,1)
        % Initialize the minimum distance to be the distance between the
        % first centroid and the data point
        minDistance = sum((X(i,:) - centroids(1,:)).^2);
        % Initialize the closest centroid to be the first centroid
        closestCentroid = 1;
        % Iterate through each centroid
        for j = 2:size(centroids,1)
            % Calculate the distance between the current centroid and the
            % data point
            distance = sum((X(i,:) - centroids(j,:)).^2);
            % If the distance is less than the minimum distance, update the
            % minimum distance and the closest centroid
            if distance < minDistance
                minDistance = distance;
                closestCentroid = j;
            end
        end
        % Set the closest centroid for the current data point
        closestCentroid(i) = closestCentroid;
    end