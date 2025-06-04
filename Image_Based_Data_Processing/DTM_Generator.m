%% Function for Transforming Radar Range-Time Image to Doppler-Time Image using STFT
% Author: JoeyBG.
% Improved By: JoeyBG.
% Time: 2025/5/25.
% Affiliation: Beijing Institute of Technology.
%
% Information:
% This function takes a radar range-time image (RTM) as input, where the rows represent
% range and the columns represent time. It sums all range cells for each time step to
% obtain a time series, and then applies the Short-Time Fourier Transform (STFT) using
% the `stft` function to this time series to produce the Doppler-time image (DTM). The DTM
% is returned as a matrix where rows correspond to Doppler frequencies and columns correspond
% to time segments.
%
% Inputs:
%   I - 2D matrix representing the RTM, with rows as range and columns as time.
%   fs - (Optional) Sampling frequency of the time series. Default is 1.
%   window - (Optional) Window function or scalar window length for STFT. Default is 256 in Hamming window.
%   noverlap - (Optional) Number of overlapping samples between windows. Default is 128.
%   nfft - (Optional) Number of FFT points. Default is 512.
%
% Outputs:
%   DTM - 2D matrix representing the Doppler-time image, with rows as Doppler frequencies
%         and columns as time segments. The values are the power spectrum (abs(S).^2).
%
% Notes:
%   - The sampling frequency `fs` is set to 1 by default, assuming the time steps in RTM are
%     uniformly spaced with a step size of 1 unit. If the actual sampling frequency is known,
%     it should be provided for accurate frequency scaling.
%   - The window is set to a Hamming window of length 256 by default, but can be customized.
%
% Reference:
% None.

%% Function Body
function DTM = DTM_Generator(I, fs, window, noverlap, nfft)
    % Ensure the input is a 2D matrix.
    if ndims(I) ~= 2
        error('Input I must be a 2D matrix.');
    end
    
    % Set default values for STFT parameters if not provided (4s slow time for RTM in [256, 256] scale).
    if nargin < 2
        fs = 64;  
    end
    if nargin < 3
        window = 32;  
    end
    if nargin < 4
        noverlap = 29;  
    end
    if nargin < 5
        nfft = 256;  
    end
    
    % Sum all range cells for each time step to get the time series.
    signal = sum(I, 1);  
    if isscalar(window)
        window = hamming(window); % If window is a scalar, create a Hamming window of that length.
    end
    
    % Compute the STFT using the stft function.
    S = stft(signal, fs, 'Window', window, 'OverlapLength', noverlap, 'FFTLength', nfft);
    DTM = abs(S).^2;
end