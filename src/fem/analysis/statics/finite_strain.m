%***********************************************
% FEMAT — Finite Element Analysis in MATLAB
%***********************************************

function [eleKt, Fint] = finite_strain(FEModel, irho, idofs, iQuad, iEle)
% Total Lagrangian finite strain: Bm'*Sm*Bm + Bm'*Fm'*C*Fm*Bm

nEleNode = size(FEModel.eleNode, 2);
ndim     = size(FEModel.nodeCoords, 2);
nEleDof  = nEleNode * ndim;

detJxW = FEModel.detJxW(iQuad, :, iEle);
dNdX   = FEModel.shape_grads(:, :, iQuad, iEle);

gradu = reshape(idofs, ndim, []) * dNdX;
[S_voigt, C_mat] = FEModel.Mat.compute_stress(gradu, ndim);
S_voigt = irho * S_voigt;
C_mat   = irho * C_mat;

F = eye(ndim) + gradu;

Bm = zeros(ndim ^ 2, nEleDof);
for i = 1:ndim
    Bm((i - 1) * ndim + 1:i * ndim, i:ndim:end) = dNdX';
end

Fm = zeros(ndim * (ndim - 1) / 2 + ndim, ndim ^ 2);
if ndim == 2
    Fm(1, 1:2:end) = F(:, 1);
    Fm(2, 2:2:end) = F(:, 2);
    Fm(3, 1:2:end) = F(:, 2);
    Fm(3, 2:2:end) = F(:, 1);
else
    Fm(1, 1:3:end) = F(:, 1);
    Fm(2, 2:3:end) = F(:, 2);
    Fm(3, 3:3:end) = F(:, 3);
    Fm(4, 1:3:end) = F(:, 2);
    Fm(4, 2:3:end) = F(:, 1);
    Fm(5, 2:3:end) = F(:, 3);
    Fm(5, 3:3:end) = F(:, 2);
    Fm(6, 1:3:end) = F(:, 3);
    Fm(6, 3:3:end) = F(:, 1);
end

if ndim == 2
    S = [S_voigt(1), S_voigt(3); S_voigt(3), S_voigt(2)];
else
    S = [S_voigt(1), S_voigt(4), S_voigt(6);
         S_voigt(4), S_voigt(2), S_voigt(5);
         S_voigt(6), S_voigt(5), S_voigt(3)];
end
Sm = kron(eye(ndim), S);

eleKt = Bm' * Sm * Bm * detJxW + Bm' * Fm' * C_mat * Fm * Bm * detJxW;
eleKt = reshape(eleKt', [], 1)';
Fint  = Bm' * Fm' * S_voigt * detJxW;

end
