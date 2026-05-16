%***********************************************
% FEMAT — Finite Element Analysis in MATLAB
%***********************************************

function [Kt, R] = assembly(FEModel, ele_func, rho, penal, dofs)
% Assemble global stiffness matrix and internal force vector

eleNode    = FEModel.eleNode;
nodeCoords = FEModel.nodeCoords;
nEle       = size(eleNode, 1);
nEleNode   = size(eleNode, 2);
ndim       = size(nodeCoords, 2);
nEleDof    = nEleNode * ndim;
nQuad      = size(FEModel.quadrature_points, 1);

eDof = FEModel.eDof;
iK   = FEModel.iK;
jK   = FEModel.jK;

E = FEModel.Mat.E;
if isfield(FEModel.Mat, 'Emin')
    Emin = FEModel.Mat.Emin;
else
    Emin = 0;
end

sKt = zeros(nEle, nEleDof ^ 2);
Fint = zeros(numel(nodeCoords), 1);

for iEle = 1:nEle
    iEleDof = eDof(iEle, :)';
    irho = ((E - Emin) * rho(iEle) ^ penal + Emin) / E;

    for iQuad = 1:nQuad
        [isKt, iFint] = ele_func(FEModel, irho, dofs(iEleDof), iQuad, iEle);
        sKt(iEle, :)  = sKt(iEle, :) + isKt;
        Fint(iEleDof) = Fint(iEleDof) + iFint;
    end
end

sKt = reshape(sKt', [], 1);
Kt = sparse(iK, jK, sKt(:));
R = FEModel.F - Fint;

end
