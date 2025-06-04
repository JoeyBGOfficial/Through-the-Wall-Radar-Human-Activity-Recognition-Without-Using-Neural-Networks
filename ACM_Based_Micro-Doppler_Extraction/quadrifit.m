%% Function to Compute Regional Constants for Four-Phase Level Set Segmentation
% Author: Chunming Li.
% Improved By: JoeyBG.
% Time: 2025/5/26.
% Affiliation: Beijing Institute of Technology.
%
% Information:
% This function computes the constants C that best fit the image U in the four regions
% defined by the signs of two level set functions in phi. It uses a smooth approximation
% via the Heaviside function with parameter epsilon. The constants C are computed as
% the averages of U over each of the four regions defined by the combinations of the
% level set functions.
%
% Inputs:
%   phi - 3D matrix (size nrow x ncol x 2), containing the two level set functions.
%   U - 2D matrix (size nrow x ncol), the image to be segmented.
%   epsilon - scalar, parameter for the Heaviside function approximation.
%   fun_n - integer, number of level set functions (expected to be 2).
%
% Outputs:
%   C - vector of length 4, containing the constants for each of the four regions.
%   mult - 4D matrix (size nrow x ncol x fun_n x M/2*(fun_n-1)), containing intermediate
%          multiplicative factors (not used in this computation but returned).
%
% Note: This function is specific to fun_n=2, computing averages over the four regions
%       defined by the two level set functions. The variable mult is set but not used
%       in the computation of C.

%% Function Body
function [C,mult]= quadrifit(phi,U,epsilon,fun_n) 
    H = Heaviside(phi,epsilon);
    s = size(H);
    M = 2^fun_n;
    C = zeros(2^fun_n,1);
    a = zeros(size(H));
    mult = ones(s(1),s(2),fun_n,M/2*(fun_n-1));
    
    mult(:,:,1,1) = 1 - H(:,:,2);
    mult(:,:,1,2) = H(:,:,2);
    mult(:,:,2,1) = 1 - H(:,:,1);
    mult(:,:,2,2) = H(:,:,1);
    
    mult2 = zeros(s(1),s(2),4);
    mult2(:,:,1) = (1-H(:,:,1)).*(1-H(:,:,2));
    mult2(:,:,2) = H(:,:,1).*(1-H(:,:,2));
    mult2(:,:,3) = (1-H(:,:,1)).*(H(:,:,2));
    mult2(:,:,4) = H(:,:,1).*H(:,:,2);
    mult3 = mult2;

    for k = 1:M
        mult2(:,:,k) = mult2(:,:,k).*U;
    end

    for k = 1:M
        tmp1 = mult3(:,:,k);
        tmp2 = mult2(:,:,k);
        denum = sum(tmp1(:));
        num = sum(tmp2(:));
        C(k) = num/denum;
    end
end