"""JAX-FEM reference solutions for Neo-Hookean beam benchmark."""

import jax
import jax.numpy as np
import numpy as onp

from jax_fem.problem import Problem
from jax_fem.generate_mesh import rectangle_mesh, box_mesh_gmsh, Mesh, get_meshio_cell_type
from jax_fem.solver import solver


class Hyperelasticity(Problem):

    def psi(self, F):
        E = 1.
        nu = 0.3
        mu = E / (2. * (1. + nu))
        kappa = E / (3. * (1. - 2. * nu))

        if self.dim == 2:
            F = np.array([[F[0,0], F[0,1], 0.],
                          [F[1,0], F[1,1], 0.],
                          [0.,       0.,   1.]])

        J = np.linalg.det(F)
        Jinv = J**(-2. / 3.)
        I1 = np.trace(F.T @ F)
        energy = (mu / 2.) * (Jinv * I1 - 3.) + (kappa / 2.) * (J - 1.)**2.
        return energy

    def get_tensor_map(self):
        P_fn = jax.grad(self.psi)
        def first_PK_stress(u_grad):
            I = np.eye(self.dim)
            F = u_grad + I
            P = P_fn(F)
            return P
        return first_PK_stress

    def get_surface_maps(self):
        if self.dim == 2:
            def surface_map(u, x):
                return np.array([0., 1e-3])
        else:
            def surface_map(u, x):
                return np.array([0., 1e-3, 0.])
        return [surface_map]


from scipy.spatial import cKDTree


# ========== 2D QUAD4 ==========
Lx, Ly = 8., 1.
Nx, Ny = 32, 4
ele_type = 'QUAD4'
cell_type = get_meshio_cell_type(ele_type)
meshio_mesh = rectangle_mesh(Nx=Nx, Ny=Ny, domain_x=Lx, domain_y=Ly)
mesh = Mesh(meshio_mesh.points, meshio_mesh.cells_dict[cell_type], ele_type=ele_type)

def left(point):
    return np.isclose(point[0], 0., atol=1e-5)
def right(point):
    return np.isclose(point[0], Lx, atol=1e-5)
def dirichlet_val(point):
    return 0.

dirichlet_bc_info = [[left, left], [0, 1], [dirichlet_val, dirichlet_val]]

problem = Hyperelasticity(mesh, vec=2, dim=2, ele_type=ele_type,
                          dirichlet_bc_info=dirichlet_bc_info,
                          location_fns=[right])

sol_2d = solver(problem, solver_options={'spsolve_solver':{'tol':1e-8}})[0]
onp.savetxt('output/sol2d.txt', sol_2d.ravel())
print(f"2D solved, dofs min={sol_2d.min():.4f} max={sol_2d.max():.4f}")


# ========== 3D HEX8 ==========
Lx, Ly, Lz = 8., 1., 1.
Nx, Ny, Nz = 80, 10, 10
ele_type = 'HEX8'
cell_type = get_meshio_cell_type(ele_type)
meshio_mesh = box_mesh_gmsh(Nx=Nx, Ny=Ny, Nz=Nz,
                             domain_x=Lx, domain_y=Ly, domain_z=Lz,
                             data_dir='./output', ele_type=ele_type)
mesh = Mesh(meshio_mesh.points, meshio_mesh.cells_dict[cell_type])

dirichlet_bc_info = [[left, left, left], [0, 1, 2], [dirichlet_val, dirichlet_val, dirichlet_val]]

problem = Hyperelasticity(mesh, vec=3, dim=3, ele_type=ele_type,
                          dirichlet_bc_info=dirichlet_bc_info,
                          location_fns=[right])

sol_3d = solver(problem, solver_options={'spsolve_solver':{}})[0]
raw_3d = sol_3d.ravel()

# align to our node ordering
import pyfemat
our_mesh = pyfemat.HexMesh(Lx, Ly, Lz, Nx, Ny, Nz)
pts_our = our_mesh.nodeCoords
aligned = onp.zeros(raw_3d.shape)
tree = cKDTree(pts_our)
_, idx = tree.query(meshio_mesh.points)
for i, j in enumerate(idx):
    aligned[3*j:3*j+3] = raw_3d[3*i:3*i+3]
onp.savetxt('output/sol3d.txt', aligned)
print(f"3D solved, dofs min={sol_3d.min():.4f} max={sol_3d.max():.4f}")
