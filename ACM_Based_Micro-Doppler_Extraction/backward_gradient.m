%% Function of Backward Gradient Calculation
% Author: Chunming Li.
% Improved By: JoeyBG.
% Time: 2025/5/26.
% Affiliation: Beijing Institute of Technology.
%
% Information:
% This function computes the backward differences of a 2D matrix f, 
% which serve as approximations to the partial derivatives along its rows and columns.
% It takes a 2D matrix f as input and returns two matrices: bdy and bdx, both of the same size as f. 
% Specifically, bdx calculates the row-wise backward differences, 
% defined as bdx(i,j) = f(i,j) - f(i-1,j) for rows 2 to the number of rows (nr), 
% with the first row set to zero since there’s no previous row to subtract. 
% This approximates the partial derivative with respect to the row index (e.g., ∂f/∂x if rows correspond to the x-direction). 
% Similarly, bdy computes the column-wise backward differences, 
% defined as bdy(i,j) = f(i,j) - f(i,j-1) for columns 2 to the number of columns (nc), 
% with the first column set to zero, 
% approximating the partial derivative with respect to the column index (e.g., ∂f/∂y if columns correspond to the y-direction). 
%
% Inputs:
%   f - 2D matrix representing the level set function or image.
%
% Outputs:
%   bdy - ∂f/∂y if columns correspond to the y-direction.
%   bdx - ∂f/∂x if rows correspond to the x-direction.
%
% Reference:
% [1] A Multiphase Level Set Framework for Image Segmentation Using the Mumford and Shah Model, IJCV 2002.
% [2] Level Set Evolution Without Reinitialization: A New Variational Formulation", CVPR 2005.

%% Function Body
function [bdy,bdx] = backward_gradient(f)
    [nr,nc] = size(f);
    bdx = zeros(nr,nc);
    bdy = zeros(nr,nc);

    bdx(2:nr,:) = f(2:nr,:) - f(1:nr-1,:);
    bdy(:,2:nc) = f(:,2:nc) - f(:,1:nc-1);
end