%% Function of Curvature Calculation using Different Difference Schemes
% Author: Chunming Li.
% Improved By: JoeyBG.
% Time: 2025/5/26.
% Affiliation: Beijing Institute of Technology.
%
% Information:
% This function computes the curvature of a 2D matrix f using different finite difference schemes
% specified by the input diff_scheme. The curvature is approximated by calculating the divergence
% of the unit normal vector to the level sets of f. It takes a 2D matrix f as input and returns
% a matrix K of the same size, representing the curvature. The function employs three distinct
% methods based on diff_scheme: 
%   - 0 uses Matlab's gradient function for initial derivatives and a combination of forward and 
%     backward gradients for the divergence, though the use of backward_gradient may only capture 
%     column-wise differences.
%   - 1 utilizes custom forward difference functions (Dx_forward, Dy_forward) and backward 
%     difference functions (Dx_backward, Dy_backward) to compute derivatives in specific directions.
%   - 2 relies on Matlab's gradient function for both initial and second-order derivatives, 
%     providing a central difference approximation.
% A small constant epsilon is used to prevent division by zero during normalization.
%
% Inputs:
%   f - 2D matrix representing the level set function or image.
%   diff_scheme - integer (0, 1, or 2) selecting the difference scheme to use.
%
% Outputs:
%   K - 2D matrix representing the approximated curvature of f.
%
% Reference:
% [1] A Multiphase Level Set Framework for Image Segmentation Using the Mumford and Shah Model, IJCV 2002.
% [2] Level Set Evolution Without Reinitialization: A New Variational Formulation", CVPR 2005.

%% Function Body
function K = CURVATURE_CV(f,diff_scheme)
    epsilon = 1e-10;
    
    if diff_scheme == 0
        [fx,fy] = gradient(f);
        fx_f = forward_gradient(f);
        ax = fx_f./sqrt(fx_f.^2+ fy.^2+epsilon);
        axx = backward_gradient(ax);
        fy_f = forward_gradient(f);
        ay = fy_f./sqrt(fx.^2 + fy_f.^2 + epsilon);
        ayy = backward_gradient(ay);
        K = axx + ayy;
    
    elseif diff_scheme == 1    
        fx_f = Dx_forward(f);
        fy_f = Dy_forward(f);
        ax = fx_f./sqrt(fx_f.^2+ fy_f.^2+epsilon);
        ay = fy_f./sqrt(fx_f.^2 + fy_f.^2 + epsilon);
        axx = Dx_backward(ax);
        ayy = Dy_backward(ay);
        K = axx + ayy;
    
    elseif diff_scheme == 2    
        [fx, fy] = gradient(f);
        ax = fx./sqrt(fx.^2+ fy.^2+epsilon);
        ay = fy./sqrt(fx.^2 + fy.^2 + epsilon);
        [axx, axy] = gradient(ax);
        [ayx, ayy] = gradient(ay);
        K = axx + ayy;
    
    else    
        disp('Wrong difference scheme: CURVATURE_CV.m');
        return;    
    end
end