%% Function of Dirac Delta Approximation for Level Set Methods
% Author: Chunming Li.
% Improved By: JoeyBG.
% Time: 2025/5/26.
% Affiliation: Beijing Institute of Technology.
%
% Information:
% This function computes an approximation to the Dirac delta function, which is widely used in 
% level set methods for applications such as image segmentation. The Dirac delta function, which 
% is zero everywhere except at the origin with an integral of one, is approximated here using a 
% smooth function. The approximation is defined as:
%   Delta_h(phi, epsilon) = (epsilon / pi) / (epsilon^2 + phi^2)
% This function accepts a 2D matrix phi, representing the level set function, and a scalar epsilon, 
% which controls the width of the approximation. It returns a matrix Delta_h of the same size as 
% phi, approximating the Dirac delta function applied to phi. This effectively concentrates the 
% function's effect near the zero level set of phi. A smaller epsilon yields a sharper 
% approximation, though it may demand a finer computational grid for stability and accuracy.
%
% Inputs:
%   phi - 2D matrix representing the level set function.
%   epsilon - scalar value determining the width of the delta function approximation.
%
% Outputs:
%   Delta_h - 2D matrix approximating the Dirac delta function applied to phi.
%
% Reference:
% [1] A Multiphase Level Set Framework for Image Segmentation Using the Mumford and Shah Model, IJCV 2002.
% [2] Level Set Evolution Without Reinitialization: A New Variational Formulation, CVPR 2005.

%% Function Body
function Delta_h = Delta(phi, epsilon)
    Delta_h = (epsilon/pi) ./ (epsilon^2 + phi.^2);
end