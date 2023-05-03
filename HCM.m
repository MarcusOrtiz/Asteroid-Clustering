function [cidx_optimal_Test, silhouette_scores]=HCM(data)
    % Minimize and seperate datasets
    subset_size = 5000;
    test_size = 10000;
    if size(data,1) > subset_size
        trainData = datasample(data, subset_size, 'Replace', false);
        testData = datasample(data, test_size, 'Replace', false);
    end

    % Approximate train distance, Nystrom method
    K = 500; % Number of landmark points
    idx_landmarks = randsample(size(trainData,1), K);
    landmarks = trainData(idx_landmarks,:);
    D_landmarks = pdist(landmarks);
    W_landmarks = squareform(exp(-D_landmarks.^2 / (2 * median(D_landmarks)^2)));
    [V, Sigma] = eig(W_landmarks);
    U = V(:,1:K);
    Sigma_half = sqrt(Sigma(1:K,1:K));
    S = U * Sigma_half;
    S_inv = pinv(S);
    W_approx = S_inv * W_landmarks * S_inv';
    D_approx = sqrt(bsxfun(@plus,diag(W_approx)',diag(W_approx))-2*W_approx);
    D_approx = abs(D_approx);
    
    tree = linkage(D_approx, 'ward');
    
    % Possible clusters
    krange = 13:25;
    
    % Find silhouettes
    silhouette_scores = zeros(length(krange), 1);
    fprintf("sil length: %d\n", length(silhouette_scores))

    for i = 1:length(krange)
        t = krange(i);
        cidx = cluster(tree, 'maxclust', t);
        if length(cidx) ~= size(trainData, 1)
            % If number of cluster assignments don't match the number of
            % data points, resize cidx to match data
            cidx = imresize(cidx, [size(trainData, 1) 1], 'nearest');
        end
        silhouette_scores(i) = mean(silhouette(trainData, cidx));
        fprintf("cluster %d\n", t)
    end
    fprintf("%d\n", silhouette_scores)
    
    % Plot the scores
    figure;
    plot(krange, silhouette_scores, '-o');
    xlim([13 25])
    xlabel('Number of Clusters');
    ylabel('Average Silhouette Score');
    figure()
    
    % Choose cluster amount with max silo
    [~, k_optimal] = max(silhouette_scores);
    fprintf('Optimal number of clusters: %d\n', krange(k_optimal));
    cidx_optimal = cluster(tree, 'maxclust', krange(k_optimal));
    
%     % Approximate testdistance, Nystrom method
%     K = 5000; % Number of landmark points
%     idx_landmarks = randsample(size(testData,1), K);
%     landmarks = testData(idx_landmarks,:);
%     D_landmarks = pdist(landmarks);
%     W_landmarks = squareform(exp(-D_landmarks.^2 / (2 * median(D_landmarks)^2)));
%     [V, Sigma] = eig(W_landmarks);
%     U = V(:,1:K);
%     Sigma_half = sqrt(Sigma(1:K,1:K));
%     S = U * Sigma_half;
%     S_inv = pinv(S);
%     W_approx = S_inv * W_landmarks * S_inv';
%     D_approxTest = sqrt(bsxfun(@plus,diag(W_approx)',diag(W_approx))-2*W_approx);
%     D_approxTest = abs(D_approxTest);
%     
%     treeTest = linkage(D_approxTest, 'ward');
%     cidx_optimal_Test = cluster(treeTest, 'maxclust', 15);
% 
%     % Plot cluster on testData
%     figure;
%     scatter3(landmarks(:,1), landmarks(:,2), landmarks(:,3), 10, cidx_optimal_Test);
%     xlabel('Eccentricity');
%     ylabel('Semi-major axis');
%     zlabel('Inclination');
%     title('HCM Clustering');
%     xlabel('Eccentricity');
%     ylabel('Semi-major axis');
%     zlabel('Inclination');
%     axis padded
%     ylim([0 200]);
%     figure()
    
end
