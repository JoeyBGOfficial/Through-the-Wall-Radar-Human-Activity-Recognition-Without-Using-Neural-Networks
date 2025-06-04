%% Function to Generate Contour Points for Initial Level Set
% Author: Chunming Li.
% Improved By: JoeyBG.
% Time: 2025/5/26.
% Affiliation: Beijing Institute of Technology.
%
% Information:
% This function generates a sequence of points that lie on the perimeter of a rectangular region
% within a grid of size nrow by ncol, inset by margin pixels from the edges. The points are ordered
% in a counter-clockwise direction starting from the top-left corner of the inset rectangle at
% (margin, margin). The contour traces the left side from top to bottom, the bottom side from
% left to right, the right side from bottom to top, and the top side from right to left. The
% contour is not closed; it does not repeat the starting point at the end.
%
% Inputs:
%   I - 2D matrix representing the image (not used in this function).
%   nrow - integer, number of rows in the grid.
%   ncol - integer, number of columns in the grid.
%   margin - integer, the number of pixels to inset from the edges of the grid.
%
% Outputs:
%   xcontour - 1D array of x-coordinates (column indices) of the contour points.
%   ycontour - 1D array of y-coordinates (row indices) of the contour points.
%
% Note: The input I is included in the function signature but is not utilized in the computation.
%       This function may be intended for initializing a level set function in a broader context.

%% Function Body
function [xcontour, ycontour] = get_phi(I, nrow, ncol, margin)
    count = 1;
    
    % Define left side of the contour, moving from top to bottom.
    x = margin;  
    for y = margin : nrow - margin + 1
        xcontour(count) = x;    
        ycontour(count) = y;    
        count = count + 1;      
    end
    
    % Define bottom side of the contour, moving from left to right.
    y = nrow - margin + 1;  
    for x = margin + 1 : ncol - margin + 1
        xcontour(count) = x;    
        ycontour(count) = y;    
        count = count + 1;      
    end
    
    % Define right side of the contour, moving from bottom to top.
    x = ncol - margin + 1; 
    for y = nrow - margin : -1 : margin
        xcontour(count) = x;    
        ycontour(count) = y;   
        count = count + 1;      
    end
    
    % Define top side of the contour, moving from right to left.
    y = margin;  
    for x = ncol - margin : -1 : margin + 1
        xcontour(count) = x;    
        ycontour(count) = y;   
        count = count + 1;    
    end
end