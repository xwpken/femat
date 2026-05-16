"""88-line topology optimization — using pyfemat as FE solver."""

import argparse
from pathlib import Path
import numpy as np
from scipy import sparse as sp
import matplotlib.pyplot as plt
import pyfemat

_OUT_DIR = Path(__file__).resolve().parent / "output"


def prepare_filter(nelx, nely, rmin):
    nEle = nelx * nely
    r = int(np.ceil(rmin)) - 1
    k, nm = 0, (2 * r + 1) ** 2
    iH, jH, sH = np.ones(nEle * nm, int), np.ones(nEle * nm, int), np.zeros(nEle * nm)
    for i1 in range(1, nelx + 1):
        for j1 in range(1, nely + 1):
            e1 = (i1 - 1) * nely + j1
            for i2 in range(max(i1 - r, 1), min(i1 + r, nelx) + 1):
                for j2 in range(max(j1 - r, 1), min(j1 + r, nely) + 1):
                    e2 = (i2 - 1) * nely + j2
                    iH[k], jH[k], sH[k] = e1, e2, max(0, rmin - np.sqrt((i1 - i2) ** 2 + (j1 - j2) ** 2))
                    k += 1
    H = sp.csr_matrix((sH[:k], (iH[:k] - 1, jH[:k] - 1)), shape=(nEle, nEle))
    Hs = np.asarray(H.sum(axis=1)).ravel()
    return H, Hs


def top88(nelx, nely, volfrac, penal, rmin, ft):
    E0, Emin, nu = 1.0, 1e-9, 0.3

    mesh = pyfemat.QuadMesh(Lx=float(nelx), Ly=float(nely), nx=nelx, ny=nely)
    model = pyfemat.Model(mesh)
    model.material("linear_elastic", stress_state="plane_stress", E=E0, nu=nu)
    model.set("Mat.Emin", Emin)
    model.analysis("small_strain")

    KE = model.get_ke(ele_idx=1)
    eDof = np.asarray(model.get("eDof"), dtype=int) - 1

    fixeddofs = list(range(1, 2 * (nely + 1) + 1, 2))
    fixeddofs.append(2 * (nelx + 1) * (nely + 1))
    model.fix(fixeddofs)
    model.force([2], [-1.0])

    H, Hs = prepare_filter(nelx, nely, rmin)

    x = np.full((nely, nelx), volfrac)
    xPhys = x.copy()
    loop, change = 0, 1.0

    fig, ax = plt.subplots()
    im = ax.imshow(1 - xPhys, cmap="gray", vmin=0, vmax=1)
    ax.axis("equal")
    ax.axis("off")
    plt.ion()
    fig.show()

    while change > 0.01:
        loop += 1
        sol = model.solve(rho=xPhys.ravel(order="F"), penal=penal)

        ue = sol.dofs[eDof]
        ce = np.sum((ue @ KE) * ue, axis=1).reshape((nely, nelx), order="F")
        c = np.sum((Emin + xPhys ** penal * (E0 - Emin)) * ce)
        dc = -penal * (E0 - Emin) * xPhys ** (penal - 1) * ce
        dv = np.ones((nely, nelx))

        if ft == 1:
            df = H @ (x.ravel(order="F") * dc.ravel(order="F")) / Hs / np.maximum(1e-3, x.ravel(order="F"))
            dc = df.reshape((nely, nelx), order="F")
        elif ft == 2:
            df = H @ (dc.ravel(order="F") / Hs)
            dv = H @ (dv.ravel(order="F") / Hs)
            dc = df.reshape((nely, nelx), order="F")
            dv = dv.reshape((nely, nelx), order="F")

        l1, l2, move = 0.0, 1e9, 0.2
        while (l2 - l1) / (l1 + l2) > 1e-3:
            lmid = 0.5 * (l2 + l1)
            xnew = np.maximum(0, np.maximum(x - move, np.minimum(1, np.minimum(x + move, x * np.sqrt(-dc / dv / lmid)))))
            if ft == 1:
                xPhys = xnew
            else:
                xPhys = (H @ xnew.ravel(order="F") / Hs).reshape((nely, nelx), order="F")
            l1, l2 = (lmid, l2) if np.sum(xPhys) > volfrac * nelx * nely else (l1, lmid)

        change = np.max(np.abs(xnew - x))
        x = xnew
        print(f" It.:{loop:5d} Obj.:{c:11.4f} Vol.:{np.mean(xPhys):7.3f} ch.:{change:7.3f}")

        im.set_data(1 - xPhys)
        fig.canvas.draw()
        plt.pause(0.001)

    _OUT_DIR.mkdir(exist_ok=True)
    fig.savefig(str(_OUT_DIR / "top88.png"), dpi=300, bbox_inches="tight", pad_inches=0)
    plt.ioff()
    plt.show(block=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for a in [("--nelx", 60), ("--nely", 20), ("--volfrac", 0.5), ("--penal", 3), ("--rmin", 2.4), ("--ft", 1)]:
        p.add_argument(a[0], type=type(a[1]), default=a[1])
    a = p.parse_args()
    top88(a.nelx, a.nely, a.volfrac, a.penal, a.rmin, a.ft)
