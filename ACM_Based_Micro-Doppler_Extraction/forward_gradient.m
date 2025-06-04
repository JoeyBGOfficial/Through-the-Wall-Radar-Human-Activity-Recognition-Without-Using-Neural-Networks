%% Function of Forward Gradient Calculation
% Author: Chunming Li.
% Improved By: JoeyBG.
% Time: 2025/5/26.
% Affiliation: Beijing Institute of Technology.
%
% Information:
% This function computes the forward differences of a 2D matrix f, 
% which serve as approximations to the partial derivatives along its rows and columns.
% It takes a 2D matrix f as input and returns two matrices: fdy and fdx, both of the same size as f. 
% Specifically, fdx calculates the row-wise forward differences, 
% defined as fdx(i,j) = f(i+1,j) - f(i,j) for rows 1 to nr-1, 
% with the last row set to zero since there’s no next row to subtract from. 
% This approximates the partial derivative with respect to the row index (e.g., ∂f/∂x if rows correspond to the x-direction). 
% Similarly, fdy computes the column-wise forward differences, 
% defined as fdy(i,j) = f(i,j+1) - f(i,j) for columns 1 to nc-1, 
% with the last column set to zero, 
% approximating the partial derivative with respect to the column index (e.g., ∂f/∂y if columns correspond to the y-direction). 
%
% Inputs:
%   f - 2D matrix representing the input data.
%
% Outputs:
%   fdy - 2D matrix representing the forward differences along the columns.
%   fdx - 2D matrix representing the forward differences along the rows.
%
% Reference:
% [1] A Multiphase Level Set Framework for Image Segmentation Using the Mumford and Shah Model, IJCV 2002.
% [2] Level Set Evolution Without Reinitialization: A New Variational Formulation, CVPR 2005.

%% Function Body
function [fdy,fdx]=forward_gradient(f)
    [nr,nc] = size(f);           
    fdx = zeros(nr,nc);          
    fdy = zeros(nr,nc);           

    a = f(2:nr,:) - f(1:nr-1,:);   % Compute row-wise forward differences: a(i,j) = f(i+1,j) - f(i,j).
    fdx(1:nr-1,:) = a;           
    b = f(:,2:nc) - f(:,1:nc-1);   % Compute column-wise forward differences: b(i,j) = f(i,j+1) - f(i,j).
    fdy(:,1:nc-1) = b;          
end