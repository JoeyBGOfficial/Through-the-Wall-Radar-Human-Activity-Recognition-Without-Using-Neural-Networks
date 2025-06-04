%% Function to Initialize Level Set Functions for Two Circles
% Author: Chunming Li.
% Improved By: JoeyBG.
% Time: 2025/5/26.
% Affiliation: Beijing Institute of Technology.
%
% Information:
% This function generates initial level set functions for two circles within a grid of size nrow by ncol.
% Each level set function has its zero level set defining a circle of radius r/2 centered at specified pixel locations.
% The level set functions are positive inside the circles and negative outside.
%
% Inputs:
%   nrow - integer, number of rows in the grid.
%   ncol - integer, number of columns in the grid.
%   ic, jc - (unused) possibly intended for center coordinates.
%   r - scalar, determines the radius of the initial circles as r/2.
%   fun_n - (unused) possibly intended for another purpose.
%   far_pixel - 1x2 vector, [col, row] coordinates of the center for the second circle.
%   near_pixel - 1x2 vector, [col, row] coordinates of the center for the first circle.
%
% Outputs:
%   f - 3D matrix (nrow x ncol x 2), containing the two level set functions.
%
% Note: The parameters ic, jc, and fun_n are included in the function signature but are not used in the computation.
%       The level set functions are computed as f = -sqrt((X-cx)^2 + (Y-cy)^2) + r/2 for each circle center (cx, cy).
%       Variables a and b are defined but not used in the current implementation.

%% Function Body
function f = initial_sdf2circle(nrow, ncol, ic, jc, r, fun_n, far_pixel, near_pixel)
    [X, Y] = meshgrid(1:ncol, 1:nrow);
    
    a = 5;
    b = 6;
    
    f = ones(nrow, ncol, 2);
    
    % Compute the level set function for the first circle centered at near_pixel.
    f(:,:,1) = -sqrt((X - near_pixel(2)).^2 + (Y - near_pixel(1)).^2) + r/2;    
    % Compute the level set function for the second circle centered at far_pixel.
    f(:,:,2) = -sqrt((X - far_pixel(2)).^2 + (Y - far_pixel(1)).^2) + r/2;
    
    % Commented-out example usage or test code (not executed).
    % f = sdf2circle(100, 50, 51, 25, 10); figure; imagesc(f);
end