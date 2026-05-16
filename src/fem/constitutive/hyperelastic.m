%***********************************************
% FEMAT — Finite Element Analysis in MATLAB
%***********************************************

function mat = hyperelastic(model, E, nu)
% Hyperelastic material factory
%   mat.compute_stress(gradu, ndim) -> [S_voigt, C_mat]
% Supported models: 'svk', 'neo_hookean'

mat.type = 'hyperelastic';
mat.model = model;
mat.E = E;
mat.nu = nu;

switch model
    case 'svk'
        mat.compute_stress = @(gradu, ndim) svk_stress(gradu, ndim, E, nu);
    case 'neo_hookean'
        mat.compute_stress = @(gradu, ndim) nh_stress(gradu, ndim, E, nu);
    otherwise
        error('Unknown hyperelastic model: %s', model)
end

end


function [S_voigt, C_mat] = svk_stress(gradu, ndim, E, nu)
% Saint Venant–Kirchhoff: 2nd PK stress and constant material tangent

F = eye(ndim) + gradu;
C_mat = D_svk(E, nu, ndim);
E_gl = 0.5 * (F' * F - eye(ndim));

if ndim == 2
    E_voigt = [E_gl(1,1); E_gl(2,2); 2 * E_gl(1,2)];
else
    E_voigt = [E_gl(1,1); E_gl(2,2); E_gl(3,3);
               2 * E_gl(1,2); 2 * E_gl(2,3); 2 * E_gl(1,3)];
end

S_voigt = C_mat * E_voigt;

end


function [S_voigt, C_mat] = nh_stress(gradu, ndim, E, nu)
% Compressible Neo-Hookean: W = (G/2)(J^{-2/3}I1 - 3) + (kappa/2)(J-1)^2

G = E / (2 * (1 + nu));
kappa = E / (3 * (1 - 2 * nu));

F = eye(ndim) + gradu;
C = F' * F;
CI = inv(C);
J = det(F);
I1 = trace(C);
if ndim == 2
    I1 = I1 + 1;  % plane strain: include F_33 = 1
end

% 2nd PK stress (tensor)
S = G * J^(-2/3) * (eye(ndim) - (1/3) * I1 * CI) + kappa * J * (J - 1) * CI;

% Convert to Voigt
if ndim == 2
    S_voigt = [S(1,1); S(2,2); S(1,2)];
else
    S_voigt = [S(1,1); S(2,2); S(3,3);
               S(1,2); S(2,3); S(1,3)];
end

% Material tangent: C_IJKL = 2 * dS_IJ / dC_KL in Voigt form
nvoigt = ndim * (ndim - 1) / 2 + ndim;
C_mat = zeros(nvoigt);
for p = 1:nvoigt
    for q = 1:nvoigt
        [i,j] = voigt2ij(p, ndim);
        [k,l] = voigt2ij(q, ndim);
        C_mat(p,q) = material_tangent(i,j,k,l, G, kappa, J, I1, CI);
    end
end

end


function C_val = material_tangent(i,j,k,l, G, kappa, J, I1, CI, ndim)

CI_ij = CI(i,j);
CI_kl = CI(k,l);
d_ij  = double(i==j);
d_kl  = double(k==l);

% Deviatoric part: C_dev = 2 * dS_dev / dC
apb = CI(i,k)*CI(j,l) + CI(i,l)*CI(j,k);

t1 = -(1/3) * CI_kl * d_ij;
t2 =  (1/9) * I1 * CI_kl * CI_ij;
t3 = -(1/3) * d_kl * CI_ij;
t4 =  (1/6) * I1 * apb;

C_dev = 2 * G * J^(-2/3) * (t1 + t2 + t3 + t4);

% Volumetric part: C_vol = 2 * dS_vol / dC
C_vol = kappa * (J * (2*J - 1) * CI_kl * CI_ij - J * (J - 1) * apb);

C_val = C_dev + C_vol;

end


function [i,j] = voigt2ij(p, ndim)
% Map Voigt index to tensor index pair

if ndim == 2
    map = [1 1; 2 2; 1 2];
else
    map = [1 1; 2 2; 3 3; 1 2; 2 3; 1 3];
end
i = map(p, 1);
j = map(p, 2);

end


function D = D_svk(E, nu, ndim)
% Material elasticity matrix (constant for SVK)

c = E / ((1 + nu) * (1 - 2 * nu));
if ndim == 2
    D = c * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2 * nu) / 2];
else
    D = c * [1 - nu, nu,   nu,   0,   0,   0;
             nu,   1 - nu, nu,   0,   0,   0;
             nu,   nu,   1 - nu, 0,   0,   0;
             0,    0,    0,     (1 - 2*nu)/2, 0, 0;
             0,    0,    0,     0, (1 - 2*nu)/2, 0;
             0,    0,    0,     0,   0, (1 - 2*nu)/2];
end

end
