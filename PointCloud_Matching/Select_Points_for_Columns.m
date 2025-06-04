%% Function to Find Top Points per Column in a Matrix
% Author: JoeyBG.
% Improved By: JoeyBG.
% Date: 2025/5/27.
% Affiliation: Beijing Institute of Technology.
%
% Information:
% This function takes a matrix phi_1_Normalized of size Estimation_Resolution x Estimation_Resolution,
% and for each column, finds the six points with the largest values, recording their row indices.
% The row indices are determined by sorting each column in descending order and selecting the
% top six indices.
%
% Inputs:
%   phi_1_Normalized - A 2D matrix of size Estimation_Resolution x Estimation_Resolution,
%                      typically containing normalized values.
%
% Outputs:
%   row_indices - A matrix of size Estimation_Resolution x 6, where each row contains the row
%                 indices of the six largest values in the corresponding column of phi_1_Normalized.
%
% Note:
%   - It is assumed that Estimation_Resolution >= 6. If Estimation_Resolution < 6, the function
%     will need adjustments to handle fewer elements per column.
%   - The function assumes phi_1_Normalized is a valid 2D matrix with real-valued elements.

%% Function Body.
function points = Select_Points_for_Columns(phi_1_Normalized,Points_Num_Per_Column)
    % Determine the size of the input matrix.
    Estimation_Resolution = size(phi_1_Normalized, 1);
    points = zeros(Points_Num_Per_Column*Estimation_Resolution, 2);
    
    % Loop over each column of the matrix.
    for i = 1:Estimation_Resolution
        col = phi_1_Normalized(:,i);
        
        % Sort the column in descending order and get the indices of the top six values.
        [~, sorted_indices] = sort(col, 'descend');
        top_indices = sorted_indices(1:Points_Num_Per_Column);
        points((i-1)*Points_Num_Per_Column+1:i*Points_Num_Per_Column,1) = top_indices;      
        points((i-1)*Points_Num_Per_Column+1:i*Points_Num_Per_Column,2) = i; 
    end
    
    points = points';
end