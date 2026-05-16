%***********************************************
% FEMAT — Finite Element Analysis in MATLAB
%***********************************************

function sol = solver(FEModel, init_dofs, opts, rho, penal)
% Newton-Raphson solver with adaptive load stepping

nDof = numel(FEModel.nodeCoords);

if nargin < 2 || isempty(init_dofs)
    init_dofs = zeros(nDof, 1);
end
if nargin < 3 || isempty(opts)
    opts = struct();
end
if nargin < 4 || isempty(rho)
    rho = ones(size(FEModel.eleNode, 1), 1);
end
if nargin < 5 || isempty(penal)
    penal = 1;
end

if ~isfield(opts, 'res_tol');  opts.res_tol = 1e-6;  end
if ~isfield(opts, 'max_iter'); opts.max_iter = 200;  end
if ~isfield(opts, 'max_cut');  opts.max_cut = 10;    end
if ~isfield(opts, 'verbose');  opts.verbose = true;   end

atype = FEModel.Atype;
if atype == "small_strain"
    ele_func = @small_strain;
elseif atype == "finite_strain"
    ele_func = @finite_strain;
else
    error('Unknown analysis type: %s', atype)
end

if opts.verbose
    nEle  = size(FEModel.eleNode, 1);
    nNode = size(FEModel.nodeCoords, 1);
    ndim  = size(FEModel.nodeCoords, 2);
    fprintf('[femat] Elements %d (%s) | Nodes %d | DOFs %d\n', ...
            nEle, FEModel.etype, nNode, nNode * ndim);
    fprintf('[femat] Newton-Raphson iteration\n')
end

fixeddofs = FEModel.fixeddofs;
freedofs  = setdiff((1:nDof)', fixeddofs);

res_tol = opts.res_tol;
max_cut = opts.max_cut;
max_iter = opts.max_iter;

dofs = init_dofs;
F_total = FEModel.F;
F_sub = zeros(nDof, 1);
F_sub_inc = F_total;
con_flag = false;
sub = 0;
cut = 1;

while ~con_flag && cut <= max_cut
    sub = sub + 1;
    loop_sub = 0;
    res_sub = 1;
    dofs_sub = dofs;
    dofs_inc_sub = zeros(nDof, 1);
    FEModel.F = F_sub + F_sub_inc;
    inc_ratio = sum(F_sub + F_sub_inc) / sum(F_total);

    while res_sub > res_tol && loop_sub < max_iter
        loop_sub = loop_sub + 1;
        [Kt, R] = assembly(FEModel, ele_func, rho, penal, dofs_sub);
        dofs_inc_sub(freedofs, :) = Kt(freedofs, freedofs) \ R(freedofs, :);
        dofs_inc_sub(fixeddofs, :) = 0;
        dofs_sub = dofs_sub + dofs_inc_sub;
        res_sub = norm(R(freedofs, :), 2) / norm(FEModel.F(freedofs, :), 2);

        if opts.verbose
            fprintf('[femat] Sub %.2d, Cut %.4f, Iter %.3d, res = %.2e\n', ...
                    sub, inc_ratio, loop_sub, res_sub);
        end
    end

    if res_sub < res_tol
        dofs = dofs_sub;
        F_sub = F_sub + F_sub_inc;
    else
        F_sub_inc = F_sub_inc / 2;
        cut = cut + 1;
    end

    if norm(F_total - F_sub) < 1e-8
        con_flag = true;
    end
end

if con_flag
    if opts.verbose
        fprintf('[femat] Converged in %d substeps\n', cut);
    end
else
    error('FEA fails to converge')
end

sol.dofs = dofs;
sol.K = Kt;
sol.R = R;

end
