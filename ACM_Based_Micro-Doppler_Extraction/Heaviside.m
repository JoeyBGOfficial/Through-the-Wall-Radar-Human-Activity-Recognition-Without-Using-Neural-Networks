%% Function to compute a smooth approximation of the Heaviside step function.
% Author: Chunming Li.
% Improved By: JoeyBG.
% Time: 2025/5/26.
% Affiliation: Beijing Institute of Technology.
%
% Information:
% This function calculates a smooth version of the Heaviside step function
% using the arctangent function. It is useful in numerical methods where
% a differentiable approximation of the step function is needed, such as
% in level set methods for image segmentation.
%
% Inputs:
%   phi - 2D matrix or scalar input values.
%   epsilon - scalar, controls the smoothness of the approximation.
%             Smaller epsilon makes the transition sharper.
%
% Outputs:
%   H - 2D matrix or scalar, the smooth Heaviside function values.
%
% The approximation is given by H = 0.5 * (1 + (2/pi) * atan(phi / epsilon)),
% which approaches 0 as phi -> -infty, 1 as phi -> +infty, and is 0.5 at phi=0.

%% Function Body
function H = Heaviside(phi, epsilon) 
    H = 0.5 * (1 + (2/pi) * atan(phi ./ epsilon)); 
end