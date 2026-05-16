import numpy as np
from . import engine as _eng


class QuadMesh:
    def __init__(self, Lx, Ly, nx, ny):
        ele, coo = _eng.call("meshtools", "rect", Lx, Ly, nx, ny, nargout=2)
        self.eleNode = np.asarray(ele)
        self.nodeCoords = np.asarray(coo)


class TriMesh:
    def __init__(self, Lx, Ly, nx, ny, pattern="left"):
        ele, coo = _eng.call("meshtools", "tri", Lx, Ly, nx, ny, pattern, nargout=2)
        self.eleNode = np.asarray(ele)
        self.nodeCoords = np.asarray(coo)


class HexMesh:
    def __init__(self, Lx, Ly, Lz, nx, ny, nz):
        ele, coo = _eng.call("meshtools", "hex", Lx, Ly, Lz, nx, ny, nz, nargout=2)
        self.eleNode = np.asarray(ele)
        self.nodeCoords = np.asarray(coo)


class TetMesh:
    def __init__(self, Lx, Ly, Lz, nx, ny, nz):
        ele, coo = _eng.call("meshtools", "tet", Lx, Ly, Lz, nx, ny, nz, nargout=2)
        self.eleNode = np.asarray(ele)
        self.nodeCoords = np.asarray(coo)
