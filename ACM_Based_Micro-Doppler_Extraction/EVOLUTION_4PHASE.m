%% Function for Evolution of Four-Phase Level Set Segmentation
% Author: Chunming Li.
% Improved By: JoeyBG.
% Time: 2025/5/26.
% Affiliation: Beijing Institute of Technology.
%
% Information:
% This function implements the evolution step for a four-phase level set segmentation method.
% It takes an image I and initial level set functions phi0 (with two components), and evolves
% the level set functions over a specified number of iterations using a combination of curvature
% and data fitting terms. The evolution is governed by parameters nu, lambda_1, lambda_2,
% with time step delta_t and epsilon for the Dirac delta approximation.
%
% The function uses Neumann boundary conditions to handle the domain boundaries and computes
% the curvature using a finite difference scheme (set to 0 in this implementation). It also
% relies on an external function quadrifit to compute coefficients for the data fitting term.
%
% Inputs:
%   I - 2D matrix representing the image to be segmented.
%   phi0 - 3D matrix (size of I x 2) containing the initial level set functions.
%   nu - scalar weight for the curvature term.
%   lambda_1, lambda_2 - scalar parameters for the data fitting term.
%   delta_t - scalar time step for the evolution.
%   epsilon - scalar for the Dirac delta approximation.
%   numIter - integer number of iterations to perform.
%
% Outputs:
%   phi - 3D matrix (size of I x 2) containing the evolved level set functions.
%
% Reference:
% [1] A Multiphase Level Set Framework for Image Segmentation Using the Mumford and Shah Model, IJCV 2002.
% [2] Level Set Evolution Without Reinitialization: A New Variational Formulation, CVPR 2005.
%
% See also: NeumannBoundCond, Delta, CURVATURE_CV, quadrifit.

%% Function Body
function phi = EVOLUTION_4PHASE(I, phi0, nu, lambda_1, lambda_2, delta_t, epsilon, numIter)
    phi(:,:,1) = phi0(:,:,1);
    phi(:,:,2) = phi0(:,:,2);
    
    for m = 1:numIter
        for k = 1:2
            tmp_phi = phi(:,:,k);
            tmp_phi = NeumannBoundCond(tmp_phi);
            delta_h = Delta(tmp_phi,epsilon);  
            differenceScheme = 0; % Use 0, 1, 2 for different schemes to compute curvature.
            Curv = CURVATURE_CV(tmp_phi,differenceScheme); 
            [C,mult] = quadrifit(phi,I,epsilon,2);
            b = zeros(size(Curv));
    
            if k==1
                b = (C(1)-C(2))*(C(1)+C(2)-2*I).*mult(:,:,1,1);
                b = b + (C(3)-C(4))*(C(3)+C(4)-2*I).*mult(:,:,1,2);
            else
                b = (C(1)-C(3))*(C(1)+C(3)-2*I).*mult(:,:,2,1);
                b = b + (C(2)-C(4))*(C(2)+C(4)-2*I).*mult(:,:,2,2);     
            end

            tmp_phi = tmp_phi + delta_t * (delta_h .* (nu * Curv + b));    
            phi(:,:,k) = tmp_phi;     
        end
    end
end

%% Helper Function for Neumann Boundary Conditions
% This function applies Neumann boundary conditions to a 2D matrix f.
% It sets the boundary values to be equal to the values two pixels inward,
% effectively making the gradient at the boundary zero.
%
% Input:
%   f - 2D matrix.
%
% Output:
%   g - 2D matrix with Neumann boundary conditions applied.

function g = NeumannBoundCond(f)
    [nrow,ncol] = size(f);
    g = f;
    g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
    g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
    g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);
end