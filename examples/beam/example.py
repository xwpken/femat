"""Cantilever beam — Neo-Hookean finite strain. Compare with JAX-FEM reference."""

import numpy as np
from pathlib import Path
import pyfemat

OUT = Path(__file__).resolve().parent / "output"

# ===== 2D QUAD4 =====
mesh = pyfemat.QuadMesh(8.0, 1.0, 32, 4)
model = pyfemat.Model(mesh)
model.material("hyperelastic", model="neo_hookean", E=1, nu=0.3)
model.analysis("finite_strain")

coords = model.get("nodeCoords")
left = np.where(coords[:, 0] < 1e-8)[0] + 1
model.fix(np.concatenate([2 * left - 1, 2 * left]).astype(int).tolist())

right = np.where(np.abs(coords[:, 0] - 8.0) < 1e-8)[0]
v = np.full(len(right), -1e-3 * 0.25)
v[0] /= 2; v[-1] /= 2
model.force((2 * (right + 1)).tolist(), v.tolist())

sol = model.solve(opts={"verbose": True})
ref = np.loadtxt(OUT / "sol2d.txt")
diff = sol.dofs - ref
i = 2 * right[0]
print("=== 2D (QUAD4) ===")
print("Max |error|: {:.2e}".format(np.max(np.abs(diff))))
print("Bottom-right: ux = {:.6f} (ref {:.6f}), uy = {:.6f} (ref {:.6f})".format(
    sol.dofs[i], ref[i], sol.dofs[i + 1], ref[i + 1]))
print()

# ===== 3D HEX8 =====
mesh = pyfemat.HexMesh(8.0, 1.0, 1.0, 80, 10, 10)
model = pyfemat.Model(mesh)
model.material("hyperelastic", model="neo_hookean", E=1, nu=0.3)
model.analysis("finite_strain")

coords = model.get("nodeCoords")
left = np.where(coords[:, 0] < 1e-8)[0] + 1
model.fix(np.concatenate([3 * left - 2, 3 * left - 1, 3 * left]).astype(int).tolist())

right = np.where(np.abs(coords[:, 0] - 8.0) < 1e-8)[0]
af = -1e-3 * 0.1 * 0.1
ld, lv = [], []
for node in right:
    y, z = coords[node, 1], coords[node, 2]
    w = 1.0
    if abs(y) < 1e-8 or abs(y - 1.0) < 1e-8: w /= 2
    if abs(z) < 1e-8 or abs(z - 1.0) < 1e-8: w /= 2
    ld.append(3 * (node + 1) - 1)
    lv.append(af * w)
model.force(ld, lv)

sol = model.solve(opts={"verbose": True})
ref = np.loadtxt(OUT / "sol3d.txt")
diff = sol.dofs - ref
i = right[0]
print("=== 3D (HEX8) ===")
print("Max |error|: {:.2e}".format(np.max(np.abs(diff))))
print("Bottom-right: ux = {:.6f} (ref {:.6f}), uy = {:.6f} (ref {:.6f}), uz = {:.6f} (ref {:.6f})".format(
    sol.dofs[3 * i], ref[3 * i], sol.dofs[3 * i + 1], ref[3 * i + 1],
    sol.dofs[3 * i + 2], ref[3 * i + 2]))

model.save_vtk(str(OUT / "beam_3d.vtk"), sol)
print("VTK saved to output/beam_3d.vtk")
