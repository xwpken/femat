from dataclasses import dataclass
import numpy as np
from . import engine as _eng


@dataclass
class Solution:
    dofs: np.ndarray
    K: np.ndarray
    R: np.ndarray


class Model:
    """FEM model that proxies to a MATLAB FEModel struct."""

    def __init__(self, mesh):
        self._var = "m" + str(id(self))
        _eng.put(self._var + "_ele", mesh.eleNode)
        _eng.put(self._var + "_coo", mesh.nodeCoords)
        _eng.execute("{} = mesh({}_ele, {}_coo);".format(self._var, self._var, self._var))

    def set(self, field, value):
        """Set any field on the MATLAB FEModel struct.
        field: 'Mat.Emin', 'Atype', 'fixeddofs', etc.
        value: Python value (int, float, str, list, np.ndarray)
        """
        vname = self._var + "_val"
        if isinstance(value, str):
            _eng.execute("{}.{} = '{}';".format(self._var, field, value))
        else:
            _eng.put(vname, value)
            _eng.execute("{}.{} = {};".format(self._var, field, vname))

    def get(self, field):
        """Get any field from the MATLAB FEModel struct as numpy array."""
        vname = self._var + "_fld"
        _eng.execute("{} = {}.{};".format(vname, self._var, field))
        return np.asarray(_eng.get(vname))

    def get_ke(self, ele_idx=1, dofs=None):
        """Get element tangent stiffness matrix (nDOF x nDOF).
        ele_idx: 1-based element index
        dofs:    element displacement vector (default zeros)
        """
        ndof = self.get("eDof").shape[1]
        if dofs is None:
            dofs = np.zeros(ndof)
        _eng.put(self._var + "_gd", np.asarray(dofs))
        _eng.put(self._var + "_ge", float(ele_idx))
        _eng.execute("{}_ke = elementools('ke', {}, {}_ge, {}_gd);".format(
            self._var, self._var, self._var, self._var))
        return np.asarray(_eng.get(self._var + "_ke"))

    def material(self, kind, **params):
        if kind == "linear_elastic":
            _eng.execute('{}.Mat = linear_elastic("{}", {}, {});'.format(
                self._var, params["stress_state"], params["E"], params["nu"]))
        elif kind == "hyperelastic":
            _eng.execute('{}.Mat = hyperelastic("{}", {}, {});'.format(
                self._var, params["model"], params["E"], params["nu"]))

    def fix(self, dofs):
        _eng.put(self._var + "_fd", np.asarray(dofs, dtype=float).reshape(-1, 1))
        _eng.execute("{}.fixeddofs = {}_fd;".format(self._var, self._var))

    def force(self, dofs, vals):
        _eng.put(self._var + "_ld", np.asarray(dofs, dtype=float).reshape(-1, 1))
        _eng.put(self._var + "_lv", np.asarray(vals, dtype=float).reshape(-1, 1))
        _eng.execute("{}_n = numel({}.nodeCoords);".format(self._var, self._var))
        _eng.execute("{}.F = sparse({}_ld, 1, {}_lv, {}_n, 1);".format(
            self._var, self._var, self._var, self._var))

    def analysis(self, atype):
        _eng.execute("{}.Atype = '{}';".format(self._var, atype))

    def solve(self, opts=None, rho=None, penal=None):
        if opts is None:
            opts = {}
        defaults = {"verbose": False, "res_tol": 1e-6, "max_iter": 200, "max_cut": 10}
        opts = {**defaults, **opts}
        _eng.put(self._var + "_opts", opts)

        if rho is not None:
            _eng.put(self._var + "_rho", np.asarray(rho, dtype=float).reshape(-1, 1))
            rho_arg = self._var + "_rho"
        else:
            rho_arg = "[]"

        if penal is not None:
            _eng.put(self._var + "_pen", float(penal))
            pen_arg = self._var + "_pen"
        else:
            pen_arg = "[]"

        _eng.execute("{}_sol = solver({}, [], {}_opts, {}, {});".format(
            self._var, self._var, self._var, rho_arg, pen_arg))
        _eng.execute("{}_dofs = {}_sol.dofs;".format(self._var, self._var))
        _eng.execute("{}_K = full({}_sol.K);".format(self._var, self._var))
        _eng.execute("{}_R = {}_sol.R;".format(self._var, self._var))
        raw = {
            "dofs": _eng.get(self._var + "_dofs"),
            "K": _eng.get(self._var + "_K"),
            "R": _eng.get(self._var + "_R"),
        }
        dofs = np.asarray(raw["dofs"]).ravel()
        K = np.asarray(raw["K"])
        R = np.asarray(raw["R"]).ravel()
        return Solution(dofs, K, R)

    def save_vtk(self, filename, sol):
        _eng.put(self._var + "_vtk_u", sol.dofs)
        _eng.execute("{}_vtk_f = '{}';".format(self._var, filename))
        _eng.execute("{}_vtk_v = ones(size({}.nodeCoords,1),1);".format(self._var, self._var))
        _eng.execute("{}_vtk_r = ones(size({}.eleNode,1),1);".format(self._var, self._var))
        _eng.execute("iotools().mat2para({}, {}_vtk_u, {}_vtk_v, {}_vtk_r, {}_vtk_f);".format(
            self._var, self._var, self._var, self._var, self._var))
