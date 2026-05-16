%***********************************************
% FEMAT — Finite Element Analysis in MATLAB
%***********************************************

function [eleKt, Fint] = small_strain(FEModel, irho, idofs, iQuad, iEle)
% Small deformation element: B'*D*B, Cauchy stress

nEleNode = size(FEModel.eleNode, 2);
ndim     = size(FEModel.nodeCoords, 2);
nEleDof  = nEleNode * ndim;

detJxW = FEModel.detJxW(iQuad, :, iEle);
dNdX   = FEModel.shape_grads(:, :, iQuad, iEle);

gradu = reshape(idofs, ndim, []) * dNdX;
[sigma_voigt, D] = FEModel.Mat.compute_stress(gradu, ndim);
D = irho * D;
sigma_voigt = irho * sigma_voigt;

B = zeros(ndim * (ndim - 1) / 2 + ndim, nEleDof);

if ndim == 2
    B(1, 1:2:end) = dNdX(:, 1)';
    B(2, 2:2:end) = dNdX(:, 2)';
    B(3, 1:2:end) = dNdX(:, 2)';
    B(3, 2:2:end) = dNdX(:, 1)';
else
    B(1, 1:3:end) = dNdX(:, 1)';
    B(2, 2:3:end) = dNdX(:, 2)';
    B(3, 3:3:end) = dNdX(:, 3)';
    B(4, 1:3:end) = dNdX(:, 2)';
    B(4, 2:3:end) = dNdX(:, 1)';
    B(5, 2:3:end) = dNdX(:, 3)';
    B(5, 3:3:end) = dNdX(:, 2)';
    B(6, 1:3:end) = dNdX(:, 3)';
    B(6, 3:3:end) = dNdX(:, 1)';
end

eleKt = B' * D * B * detJxW;
eleKt = reshape(eleKt', [], 1)';
Fint  = B' * sigma_voigt * detJxW;

end
