%% Function to Compute 2-Wasserstein Distance Between Two 2D Point Clouds
% Author: JoeyBG.
% Improved By: JoeyBG.
% Time: 2025/5/28
% Affiliation: Beijing Institute of Technology.
%
% Information:
% This function computes the 2-Wasserstein distance between two 2D point clouds
% represented as 2xN and 2xM matrices, where the first dimension represents the
% 2D coordinates and the second dimension represents the number of points.
% The Wasserstein distance is computed using the linear programming formulation
% of the optimal transport problem, with the cost being the squared Euclidean
% distance between points. The point clouds are treated as empirical uniform
% distributions, with masses 1/N and 1/M for each point in the respective clouds.
%
% Inputs:
%   pointCloud1 - 2xN matrix representing the first point cloud.
%   pointCloud2 - 2xM matrix representing the second point cloud.
%
% Outputs:
%   similarity - scalar representing the 2-Wasserstein distance between the two point clouds.
%
% Reference:
%   [1] Rubner, Y., Tomasi, C., & Guibas, L. J. (2000). 
%       The Earth Mover’s Distance as a Metric for Image Retrieval. 
%       International Journal of Computer Vision, 40(2), 99-121.

%% Function Body
function similarity = Wasserstein_Similarity(pointCloud1, pointCloud2)
    % Check inputs for correct dimensions.
    if size(pointCloud1,1) ~= 2 || size(pointCloud2,1) ~= 2
        error('Point clouds must be 2xN and 2xM matrices.');
    end
    N = size(pointCloud1,2);
    M = size(pointCloud2,2);
    
    % Compute pairwise squared Euclidean distances as cost matrix.
    C = pdist2(pointCloud1', pointCloud2', 'euclidean').^2;
    c = C(:); % Vectorize cost matrix for linear programming.
    
    % Create equality constraint matrices using sparse representation.
    row_idx1 = kron(1:N, ones(1,M));
    col_idx1 = repmat((0:M-1)*N, 1, N) + repmat(1:N,1,M);
    col_idx1 = col_idx1(:)';
    A_eq1 = sparse(row_idx1, col_idx1, 1, N, N*M);
    row_idx2 = kron(1:M, ones(1,N));
    col_idx2 = repmat(1:N,1,M) + kron((0:M-1)*N, ones(1,N));
    col_idx2 = col_idx2(:)';
    A_eq2 = sparse(row_idx2, col_idx2, 1, M, N*M);
    
    % Combine equality constraints.
    A_eq = [A_eq1; A_eq2];
    b_eq = [ones(N,1)/N; ones(M,1)/M]; % Marginal constraints for uniform distributions.
    lb = zeros(N*M,1);
    
    % Solve the optimal transport problem using linear programming.
    options = optimoptions('linprog','Display','none');
    [~, fval] = linprog(c, [], [], A_eq, b_eq, lb, [], options);
    similarity = sqrt(fval);
end