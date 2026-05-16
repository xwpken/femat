%***********************************************
% FEMAT — Finite Element Analysis in MATLAB
%***********************************************

%% Setup: run 'setup' from project root

% ===== 2D QUAD4 =====
Lx = 8; Ly = 1; Nx = 32; Ny = 4;
[ele, coo] = meshtools().rectangle_mesh(Lx, Ly, Nx, Ny);
FEModel = mesh(ele, coo);
FEModel.Atype = 'finite_strain';
FEModel.Mat = hyperelastic('neo_hookean', 1, 0.3);
FEModel.fixeddofs = [2 * find(coo(:,1) < 1e-8) - 1; 2 * find(coo(:,1) < 1e-8)];

r = find(abs(coo(:,1) - Lx) < 1e-8);
v = ones(length(r), 1) * (-1e-3) * (Ly/Ny); v(1) = v(1)/2; v(end) = v(end)/2;
FEModel.F = sparse(2 * r, 1, v, numel(coo), 1);

sol = solver(FEModel, [], struct('verbose', true));
ref = readmatrix(fullfile('output', 'sol2d.txt'));
err = abs(sol.dofs - ref);
fprintf('=== 2D (QUAD4) ===\nMax |error|: %.2e\n', max(err));
fprintf('Bottom-right: ux = %.6f (ref %.6f), uy = %.6f (ref %.6f)\n\n', ...
        sol.dofs(2*r(1)-1), ref(2*r(1)-1), sol.dofs(2*r(1)), ref(2*r(1)));

% ===== 3D HEX8 =====
[ele, coo] = meshtools().hex_mesh(8, 1, 1, 80, 10, 10);
FEModel = mesh(ele, coo);
FEModel.Atype = 'finite_strain';
FEModel.Mat = hyperelastic('neo_hookean', 1, 0.3);
left = find(coo(:,1) < 1e-8);
FEModel.fixeddofs = [3*left-2; 3*left-1; 3*left];

r = find(abs(coo(:,1) - 8) < 1e-8);
af = -1e-3 * (1/10) * (1/10);
FEModel.F = zeros(numel(coo), 1);
for i = 1:length(r)
    n = r(i); y = coo(n,2); z = coo(n,3);
    w = 1;
    if abs(y) < 1e-8 || abs(y-1) < 1e-8; w = w/2; end
    if abs(z) < 1e-8 || abs(z-1) < 1e-8; w = w/2; end
    FEModel.F(3*n-1) = af * w;
end

sol = solver(FEModel, [], struct('verbose', true));
ref = readmatrix(fullfile('output', 'sol3d.txt'));
err = abs(sol.dofs - ref);
fprintf('=== 3D (HEX8) ===\nMax |error|: %.2e\n', max(err));
fprintf('Bottom-right: ux = %.6f (ref %.6f), uy = %.6f (ref %.6f), uz = %.6f (ref %.6f)\n', ...
        sol.dofs(3*r(1)-2), ref(3*r(1)-2), sol.dofs(3*r(1)-1), ref(3*r(1)-1), sol.dofs(3*r(1)), ref(3*r(1)));

% VTK
iotools().mat2para(FEModel, sol.dofs, ones(size(coo,1),1), ones(size(ele,1),1), 'output/beam_3d.vtk');
fprintf('VTK saved to output/beam_3d.vtk\n');
