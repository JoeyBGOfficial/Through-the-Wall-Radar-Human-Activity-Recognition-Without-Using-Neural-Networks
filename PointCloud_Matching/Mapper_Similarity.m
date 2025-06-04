%% Function for Computing Topological Similarity between Two Point Clouds using Mapper Algorithm
% Author: JoeyBG.
% Improved By: JoeyBG.
% Time: 2025/5/28.
% Affiliation: Beijing Institute of Technology.
%
% Information:
% This function computes the topological similarity between two 2D point clouds, PC and
% PC_Class, using a simplified Mapper algorithm. The point clouds are represented as 2xN
% and 2xM matrices, where the first dimension contains x and y coordinates, and the
% second dimension represents the number of points. The filter function is interpreted
% as using the 2D coordinates directly, with a covering of overlapping squares in R^2.
% The Mapper graph for each point cloud is constructed by defining nodes as non-empty
% intersections with covering squares and edges where points exist in the overlap of
% adjacent squares. Similarity is measured using the Jaccard index of the edge sets.
%
% Inputs:
%   PC - 2xN matrix representing the first point cloud, where each column is [x; y].
%   PC_Class - 2xM matrix representing the second point cloud.
%   nx - integer, number of squares along the x-axis in the covering grid (default: 100).
%   ny - integer, number of squares along the y-axis in the covering grid (default: 100).
%   overlap_factor - scalar, factor determining square size relative to step size for overlap (default: 1.5).
%
% Outputs:
%   similarity - scalar, Jaccard similarity between the edge sets of the two Mapper graphs.
%
% Note:
% The filter function specification ("flattening 2D matrices into 1D vectors") is ambiguous
% for point-wise application in Mapper. Here, it is implemented as using the 2D coordinates
% with a grid covering, which aligns with typical Mapper usage. The grid size and overlap
% can be adjusted via nx, ny, and overlap_factor to tune the topological resolution.

%% Function Body
function similarity = Mapper_Similarity(PC, PC_Class, nx, ny, overlap_factor)
    % Set default parameters if not provided.
    if nargin < 3
        nx = 100;        % Default number of grid squares along x-axis.
    end
    if nargin < 4
        ny = 100;        % Default number of grid squares along y-axis.
    end
    if nargin < 5
        overlap_factor = 1.5;  % Default overlap factor for square size.
    end
    
    % Compute bounding box encompassing both point clouds.
    min_x = min([PC(1,:), PC_Class(1,:)]);
    max_x = max([PC(1,:), PC_Class(1,:)]);
    min_y = min([PC(2,:), PC_Class(2,:)]);
    max_y = max([PC(2,:), PC_Class(2,:)]);
    
    % Calculate step sizes for grid spacing.
    step_x = (max_x - min_x) / (nx - 1);
    step_y = (max_y - min_y) / (ny - 1);
    
    % Determine square sizes to ensure overlap.
    s_x = step_x * overlap_factor;
    s_y = step_y * overlap_factor;
    
    % Initialize edge list for Mapper graph construction.
    edges = {};
    edge_index = 1;
    
    % Define horizontal edges: connect U_{i,j} to U_{i+1,j}.
    for i = 0:nx-2
        for j = 0:ny-1
            % Compute overlap region between adjacent horizontal squares.
            xmin = min_x + (i+1)*step_x - s_x/2;
            xmax = min_x + i*step_x + s_x/2;
            ymin = min_y + j*step_y - s_y/2;
            ymax = min_y + j*step_y + s_y/2;
            edges{edge_index} = struct('type', 'horizontal', ...
                                       'indices', [i,j], ...
                                       'overlap', [xmin, xmax, ymin, ymax]);
            edge_index = edge_index + 1;
        end
    end
    
    % Define vertical edges: connect U_{i,j} to U_{i,j+1}.
    for i = 0:nx-1
        for j = 0:ny-2
            % Compute overlap region between adjacent vertical squares.
            xmin = min_x + i*step_x - s_x/2;
            xmax = min_x + i*step_x + s_x/2;
            ymin = min_y + (j+1)*step_y - s_y/2;
            ymax = min_y + j*step_y + s_y/2;
            edges{edge_index} = struct('type', 'vertical', ...
                                       'indices', [i,j], ...
                                       'overlap', [xmin, xmax, ymin, ymax]);
            edge_index = edge_index + 1;
        end
    end
    
    % Total number of possible edges.
    num_edges = length(edges);
    
    % Compute edge presence for PC.
    E_PC = false(1, num_edges);
    for k = 1:num_edges
        overlap = edges{k}.overlap;
        % Check if any point in PC lies within the overlap region.
        in_overlap = PC(1,:) >= overlap(1) & PC(1,:) <= overlap(2) & ...
                     PC(2,:) >= overlap(3) & PC(2,:) <= overlap(4);
        if any(in_overlap)
            E_PC(k) = true;
        end
    end
    
    % Compute edge presence for PC_Class.
    E_PC_Class = false(1, num_edges);
    for k = 1:num_edges
        overlap = edges{k}.overlap;
        % Check if any point in PC_Class lies within the overlap region.
        in_overlap = PC_Class(1,:) >= overlap(1) & PC_Class(1,:) <= overlap(2) & ...
                     PC_Class(2,:) >= overlap(3) & PC_Class(2,:) <= overlap(4);
        if any(in_overlap)
            E_PC_Class(k) = true;
        end
    end
    
    % Compute Jaccard similarity between edge sets.
    intersection = sum(E_PC & E_PC_Class);
    union = sum(E_PC | E_PC_Class);
    if union > 0
        similarity = intersection / union;
    else
        similarity = 0;  % Handle case where no edges exist in either graph.
    end
end