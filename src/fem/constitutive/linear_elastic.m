%***********************************************
% FEMAT — Finite Element Analysis in MATLAB
%***********************************************

function mat = linear_elastic(stress_state, E, nu)
% Linear elastic material factory
%   mat.compute_stress(gradu, ndim) -> [sigma_voigt, D]

mat.type = 'linear_elastic';
mat.stress_state = stress_state;
mat.E = E;
mat.nu = nu;
mat.compute_stress = @(gradu, ndim) compute_stress(gradu, ndim, E, nu, stress_state);

end


function [sigma_voigt, D] = compute_stress(gradu, ndim, E, nu, stress_state)
% Compute Cauchy stress and small-strain tangent from displacement gradient

D = D_matrix(E, nu, ndim, stress_state);
epsilon = 0.5 * (gradu + gradu');

if ndim == 2
    epsilon_voigt = [epsilon(1,1); epsilon(2,2); 2 * epsilon(1,2)];
else
    epsilon_voigt = [epsilon(1,1); epsilon(2,2); epsilon(3,3);
                     2 * epsilon(1,2); 2 * epsilon(2,3); 2 * epsilon(1,3)];
end

sigma_voigt = D * epsilon_voigt;

end


function D = D_matrix(E, nu, ndim, stress_state)
% Elasticity matrix D in Voigt notation

if ndim == 2
    if strcmp(stress_state, 'plane_strain')
        c = E / ((1 + nu) * (1 - 2 * nu));
        D = c * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2 * nu) / 2];
    else
        c = E / (1 - nu ^ 2);
        D = c * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
    end
else
    c = E / ((1 + nu) * (1 - 2 * nu));
    D = c * [1 - nu, nu,   nu,   0,   0,   0;
             nu,   1 - nu, nu,   0,   0,   0;
             nu,   nu,   1 - nu, 0,   0,   0;
             0,    0,    0,     (1 - 2*nu)/2, 0, 0;
             0,    0,    0,     0, (1 - 2*nu)/2, 0;
             0,    0,    0,     0,   0, (1 - 2*nu)/2];
end

end
