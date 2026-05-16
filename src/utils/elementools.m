%***********************************************
% FEMAT — Finite Element Analysis in MATLAB
%***********************************************

function varargout = elementools(what, varargin)
% Element-level utilities
%   KE = elementools('ke', FEModel, iEle, dofs)  -> element tangent stiffness

switch what
    case 'ke'
        varargout{1} = get_ke(varargin{:});
end

end


function KE = get_ke(FEModel, iEle, dofs)

atype = FEModel.Atype;
if atype == "small_strain"
    f = @small_strain;
elseif atype == "finite_strain"
    f = @finite_strain;
else
    error('elementools: unknown analysis type %s', atype)
end

nq = size(FEModel.quadrature_points, 1);
n  = size(FEModel.eDof, 2);
KE = zeros(n);

for i = 1:nq
    [k, ~] = f(FEModel, 1, dofs, i, iEle);
    KE = KE + reshape(k, n, n);
end

end
