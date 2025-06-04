%% Function to Detect SIFT Corner Points and Compute Farthest and Nearest Pixels to the Centroid
% Author: JoeyBG.
% Improved By: JoeyBG.
% Time: 2025/5/26.
% Affiliation: Beijing Institute of Technology.
%
% Information:
% This function detects SIFT (Scale-Invariant Feature Transform) corner points in a grayscale image,
% computes the centroid of these points, and identifies the pixels farthest and nearest to this centroid.
%
% Inputs:
%   I - 2D matrix representing the grayscale image.
%
% Outputs:
%   farPixel - 1x2 vector [row, col], the pixel farthest from the average position of the corner points.
%   nearPixel - 1x2 vector [row, col], the pixel nearest to the centroid of the corner points.
%
% Note: This function requires the Computer Vision Toolbox for SIFT feature detection.

%% Function Body
function [farPixel, nearPixel, II_Normalized, locations, farRow, farCol, nearRow, nearCol] = Corner_Representation(I,Corner_Threshold_Ratio)
    % Image preprocessing for corner detection.
    II = double(I);  
    maxPixel = max(I(:));
    threshold = Corner_Threshold_Ratio * maxPixel;
    I(I < threshold) = 0;

    % Corner representation based on SIFT method.
    points = detectSIFTFeatures(I);
    strongPoints = selectStrongest(points, 30); % 30 strongest points are used based on human motion kinematic model.

    % Find the nearest and farest central points.
    locations = strongPoints.Location;
    numPoints = size(locations, 1);
    centroid = mean(locations, 1);
    [N, M] = size(I);
    [X, Y] = meshgrid(1:M, 1:N);
    avgDist = zeros(N, M);
    for k = 1:numPoints
        dist_k = sqrt((X - locations(k,1)).^2 + (Y - locations(k,2)).^2);
        avgDist = avgDist + dist_k;
    end
    avgDist = avgDist / numPoints;
    [~, maxIdx] = max(avgDist(:));
    [farRow, farCol] = ind2sub([N, M], maxIdx);
    distToCentroid = sqrt((X - centroid(1)).^2 + (Y - centroid(2)).^2);
    [~, minIdx] = min(distToCentroid(:));
    [nearRow, nearCol] = ind2sub([N, M], minIdx);
    farPixel = [farRow, farCol];
    nearPixel = [nearRow, nearCol];

    % Prepare image data for computations.
    II_max = max(max(II));
    II_min = min(min(II));
    II_Normalized = (II-II_min)/(II_max-II_min);
end