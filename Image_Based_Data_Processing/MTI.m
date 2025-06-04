%% Function for MTI Filtering of Radar Range-Time Image (RTM)
% Author: JoeyBG.
% Improved By: JoeyBG.
% Time: 2025/05/25.
% Affiliation: Beijing Institute of Technology.
%
% Information:
% This function applies a simple moving target indication (MTI) filter to a radar
% range-time image (RTM) by subtracting adjacent columns. The input
% image I has rows representing range and columns representing time. The output
% is the filtered image where stationary targets are suppressed, and moving targets
% are highlighted.
%
% Inputs:
%   I - 2D matrix representing the RTM, with rows as range and
%       columns as time.
%
% Outputs:
%   filtered_I - 2D matrix representing the MTI filtered image.
%
% Reference:
% None.

%% Function Body.
function filtered_I = MTI(I)
    % Ensure the input is a 2D matrix.
    if ndims(I) ~= 2
        error('Input must be a 2D matrix.');
    end
    
    % Initialize the filtered image with zeros.
    [num_range, num_time] = size(I);
    filtered_I = zeros(num_range, num_time - 1);
    
    % Apply MTI filter: subtract adjacent columns.
    for t = 1:num_time - 1
        filtered_I(:, t) = I(:, t + 1) - I(:, t);
    end
end