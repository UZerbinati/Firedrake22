from firedrake import *
from mpi4py import MPI
from netgen.geom2d import SplineGeometry
from netgen.occ import *
import netgen
import numpy as np
import pytest

def poisson3D(h, degree=2):
    comm = MPI.COMM_WORLD
    # Setting up Netgen geometry and mesh
    if comm.Get_rank() == 0:
        box = Box(Pnt(0, 0, 0), Pnt(np.pi, np.pi, np.pi))
        box.bc("bcs")
        geo = OCCGeometry(box)
        ngmesh = geo.GenerateMesh(maxh=h)
        labels = ngmesh.GetBCIDs("bcs")
    else:
        ngmesh = netgen.libngpy._meshing.Mesh(3)
        labels = None

    labels = comm.bcast(labels, root=0) 
    msh = Mesh(ngmesh)

    # Setting up the problem
    V = FunctionSpace(msh, "CG", degree)
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Function(V)
    x, y, z = SpatialCoordinate(msh)
    f.interpolate(3*sin(x)*sin(y)*sin(z))
    a = inner(grad(u), grad(v))*dx
    l = inner(f, v) * dx
    u = Function(V)
    bc = DirichletBC(V, 0.0, labels)

    # Assembling matrix
    A = assemble(a, bcs=bc)
    b = assemble(l)
    bc.apply(b)

    # Solving the problem
    solve(A, u, b, solver_parameters={"ksp_type": "preonly", "pc_type": "lu"})

    # Computing the error
    f.interpolate(sin(x)*sin(y)*sin(z))
    S = sqrt(assemble(inner(u - f, u - f) * dx))
    return S

diff = np.array([poisson3D(h) for h in [2,1,1/2]])
print("l2 error norms:", diff)
conv = np.log2(diff[:-1] / diff[1:])
print("convergence order:", conv)
assert (np.array(conv) > 2.8).all()
