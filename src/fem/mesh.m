%***********************************************
% FEMAT — Finite Element Analysis in MATLAB
%***********************************************

function FEModel = mesh(eleNode, nodeCoords)
% Preprocess mesh: Jacobian, shape grads, sparse DOF indexing

nEle     = size(eleNode, 1);
nEleNode = size(eleNode, 2);
ndim     = size(nodeCoords, 2);
nEleDof  = nEleNode * ndim;

if nEleNode == 4 && ndim == 2
    etype = 'QUAD4';
elseif nEleNode == 3
    etype = 'TRI3';
elseif nEleNode == 8 && ndim == 3
    etype = 'HEX8';
elseif nEleNode == 4 && ndim == 3
    etype = 'TET4';
else
    error('Unknown element type')
end

[quad_pts, weights, shape_func, dN_dxi] = basis(etype);
nQuad = size(quad_pts, 1);

shape_grads = zeros(nEleNode, ndim, nQuad, nEle);
J = zeros(ndim, ndim, nQuad, nEle);
detJ = zeros(nQuad, 1, nEle);
detJxW = zeros(nQuad, 1, nEle);

for iEle = 1:nEle
    eCoords = nodeCoords(eleNode(iEle, :), :);
    for iQuad = 1:nQuad
        for j = 1:ndim
            J(j, :, iQuad, iEle) = sum(dN_dxi(:, :, iQuad) .* eCoords(:, j), 1);
        end
        detJ(iQuad, :, iEle) = det(J(:, :, iQuad, iEle));
        detJxW(iQuad, :, iEle) = detJ(iQuad, :, iEle) * weights(iQuad);
        shape_grads(:, :, iQuad, iEle) = dN_dxi(:, :, iQuad) / J(:, :, iQuad, iEle);
    end
end

eDof = bsxfun(@plus, kron(eleNode, ndim * ones(1, ndim)), ...
              repmat((-(ndim - 1)):0, 1, nEleNode));
iK = repelem(eDof, 1, nEleDof * ones(nEleDof, 1));
iK = reshape(iK', [], 1);
jK = repmat(eDof, 1, nEleDof);
jK = reshape(jK', [], 1);

FEModel = struct();
FEModel.etype       = etype;
FEModel.eleNode     = eleNode;
FEModel.nodeCoords  = nodeCoords;
FEModel.J           = J;
FEModel.detJxW      = detJxW;
FEModel.shape_func  = shape_func;
FEModel.shape_grads = shape_grads;
FEModel.quadrature_points = quad_pts;
FEModel.eDof        = eDof;
FEModel.iK          = iK;
FEModel.jK          = jK;

end
