%***********************************************
% FEMAT — Finite Element Analysis in MATLAB
%***********************************************

function func = iotools
% I/O utilities: VTK export

func.mat2para = @matlab2paraview;

end


function matlab2paraview(FEModel, dofs, var, rho, filename)
% Export FEM results to VTK unstructured grid format

eleNode    = FEModel.eleNode;
nodeCoords = FEModel.nodeCoords;

if size(nodeCoords, 2) == 2
    dofs = reshape(dofs, 2, [])';
    dofs = [dofs, zeros(size(dofs, 1), 1)];
    nodeCoords = [nodeCoords, zeros(size(nodeCoords, 1), 1)];
else
    dofs = reshape(dofs, 3, [])';
end

nEle     = size(eleNode, 1);
nEleNode = size(eleNode, 2);
nNode    = size(nodeCoords, 1);
eleNodeInfo = [nEleNode * ones(nEle, 1), eleNode - 1];

if nEleNode == 3;      ct = 5;
elseif nEleNode == 4
    if size(nodeCoords, 2) == 2; ct = 9;  else; ct = 10; end
elseif nEleNode == 8;  ct = 12;
end

fid = fopen(filename, 'w');
fprintf(fid, '# vtk DataFile Version 2.0\nVTK from Matlab\nASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n\n');
fprintf(fid, 'POINTS %d double\n', nNode);
fprintf(fid, '%.6f %.6f %.6f\n', nodeCoords');
fprintf(fid, '\nCELLS %d %d\n', nEle, (nEleNode + 1) * nEle);
fprintf(fid, '%d %d %d %d %d\n', eleNodeInfo');
fprintf(fid, '\nCELL_TYPES %d\n', nEle);
fprintf(fid, '%d\n', ct * ones(nEle, 1));
fprintf(fid, '\nPOINT_DATA %d\n', nNode);
fprintf(fid, 'SCALARS var double 1\nLOOKUP_TABLE default\n');
fprintf(fid, '%.6f\n', var');
fprintf(fid, 'SCALARS dofs double 3\nLOOKUP_TABLE default\n');
fprintf(fid, '%.6f %.6f %.6f\n', dofs');
fprintf(fid, '\nCELL_DATA %d\n', nEle);
fprintf(fid, 'SCALARS rho double 1\nLOOKUP_TABLE default\n');
fprintf(fid, '%.6f\n', rho');
fclose(fid);

end
