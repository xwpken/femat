%***********************************************
% FEMAT — Finite Element Analysis in MATLAB
%***********************************************

function varargout = meshtools(varargin)
% Structured mesh generation utilities
%   meshtools()       -> struct of function handles
%   meshtools('rect', Lx, Ly, nx, ny) -> [ele, coo]
%   meshtools('tri',  Lx, Ly, nx, ny, pat) -> [ele, coo]
%   meshtools('hex',  Lx, Ly, Lz, nx, ny, nz) -> [ele, coo]
%   meshtools('tet',  Lx, Ly, Lz, nx, ny, nz) -> [ele, coo]

if nargin == 0
    func.rectangle_mesh = @rectangle_mesh;
    func.triangle_mesh  = @triangle_mesh;
    func.hex_mesh       = @hex_mesh;
    func.tet_mesh       = @tet_mesh;
    varargout = {func};
    return
end

switch varargin{1}
    case 'rect'
        [a, b] = rectangle_mesh(varargin{2:end});
    case 'tri'
        [a, b] = triangle_mesh(varargin{2:end});
    case 'hex'
        [a, b] = hex_mesh(varargin{2:end});
    case 'tet'
        [a, b] = tet_mesh(varargin{2:end});
    otherwise
        error('Unknown mesh type: %s', varargin{1})
end
varargout = {a, b};

end


function [eleNode, nodeCoords] = rectangle_mesh(Lx, Ly, Nx, Ny)
% Generate QUAD4 mesh for a 2D rectangular domain

ex = Lx / Nx;
ey = Ly / Ny;
[x, y] = meshgrid(0:ex:Lx, 0:ey:Ly);
nodeCoords = [x(:), y(:)];

eleN1 = repmat((1:Ny)', 1, Nx) + kron(0:Nx - 1, (Ny + 1) * ones(Ny, 1));
eleNode = repmat(eleN1(:), 1, 4) + repmat([0, Ny + [1, 2], 1], Nx * Ny, 1);

end


function [eleNode, nodeCoords] = triangle_mesh(Lx, Ly, Nx, Ny, pattern)
% Generate TRI3 mesh by splitting QUAD4 (pattern: 'left' or 'right')

[eleNode, nodeCoords] = rectangle_mesh(Lx, Ly, Nx, Ny);

if pattern == "left"
    order = [1, 2, 4, 2, 3, 4];
    eleNode = reshape(eleNode(:, order)', 3, [])';
elseif pattern == "right"
    order = [1, 2, 3, 1, 3, 4];
    eleNode = reshape(eleNode(:, order)', 3, [])';
else
    error('Unknown triangle pattern: %s', pattern)
end

end


function [eleNode, nodeCoords] = hex_mesh(Lx, Ly, Lz, Nx, Ny, Nz)
% Generate HEX8 mesh for a 3D brick domain

[X, Y, Z] = meshgrid(0:Lx/Nx:Lx, 0:Ly/Ny:Ly, 0:Lz/Nz:Lz);
nodeCoords = [X(:), Y(:), Z(:)];

faceN1  = repmat((1:Ny)', 1, Nx) + kron(0:Nx - 1, (Ny + 1) * ones(Ny, 1));
nfaceNode = (Ny + 1) * (Nx + 1);
faceEleNode = repmat(faceN1(:), 1, 4) + repmat([0, Ny + [1, 2], 1], Nx * Ny, 1);
nfaceEle = size(faceEleNode, 1);

off = reshape([(0:nfaceNode:(Nz-1)*nfaceNode);
               (nfaceNode:nfaceNode:nfaceNode*Nz)], [], 1);
off = reshape(off, [1, 1, 2, Nz]);
eleNode = bsxfun(@plus, repmat(faceEleNode, [1, 1, 2, Nz]), off);
eleNode = reshape(eleNode, nfaceEle, 8, []);
eleNode = reshape(permute(eleNode, [1, 3, 2]), [], 8);

end


function [eleNode, nodeCoords] = tet_mesh(Lx, Ly, Lz, Nx, Ny, Nz)
% Generate TET4 mesh by subdividing HEX8 into 6 tets each

[eleNode, nodeCoords] = meshtools().hex_mesh(Lx, Ly, Lz, Nx, Ny, Nz);
order = [2, 4, 1, 6; 2, 3, 4, 6; 3, 7, 4, 6;
         7, 8, 4, 6; 8, 5, 4, 6; 1, 4, 5, 6]';
eleNode = reshape(eleNode(:, order(:))', 4, [])';

end
