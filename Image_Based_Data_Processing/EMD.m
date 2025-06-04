%% Function of RTM Denoising using EMD method
% Author: JoeyBG.
% Improved By: JoeyBG.
% Time: 2025/5/25.
% Affiliation: Beijing Institute of Technology.
%
% Information:
% EMD_DENOISE_RTM Denoise radar range-time image using Empirical Mode Decomposition
% denoised_RTM = EMDM(RTM, num_discard) applies EMD to each column
% of the range-time image RTM, discards the first num_discard IMFs, and 
% reconstructs the signal using the remaining IMFs and the residual to produce
% a denoised version of the RTM.
%
% Inputs:
%   RTM - 2D matrix, range-time image where rows represent range and columns represent time.
%   num_discard - integer, number of initial IMFs to discard as noise (default: 1).
%
% Outputs:
%   denoised_RTM - 2D matrix, denoised range-time image of the same size as RTM.
%
% Note: 
%   - This function requires the Signal Processing Toolbox for the 'emd' function.
%   - The parameter num_discard should be chosen based on the noise characteristics
%     of the data; typically, the first few IMFs contain high-frequency noise.
%   - If the number of IMFs in a column is less than or equal to num_discard,
%     only the residual is kept.

%% Function Body
function denoised_RTM = EMD(RTM, num_discard)
    % Set default value for num_discard if not provided.
    if nargin < 2
        num_discard = 1;
    end
    
    % Get dimensions of the RTM.
    [num_range, num_time] = size(RTM);
    
    % Initialize the denoised RTM matrix.
    denoised_RTM = zeros(num_range, num_time);
    
    % Process each column (time series) of the RTM.
    for t = 1:num_time
        column = RTM(:, t);
        [imf, residual] = emd(column);
        num_imfs = size(imf, 2);
        
        % Reconstruct the denoised column.
        if num_imfs <= num_discard
            denoised_column = residual;
        else
            denoised_column = sum(imf(:, num_discard+1:end), 2) + residual;
        end
        
        % Store the denoised column in the output matrix.
        denoised_RTM(:, t) = denoised_column;
    end
end