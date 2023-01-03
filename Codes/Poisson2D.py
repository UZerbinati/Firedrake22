from firedrake import *
from mpi4py import MPI
from netgen.geom2d import SplineGeometry
from netgen.occ import *
import netgen
import numpy as np
comm = MPI.COMM_WORLD
# Setting up Netgen geometry and mesh
if comm.Get_rank() == 0:
    geo = SplineGeometry()
    geo.AddRectangle((0, 0), (np.pi, np.pi), bc="rect")
    ngmesh = geo.GenerateMesh(maxh=0.1)
    labels = ngmesh.GetBCIDs("rect")
else:
    ngmesh = netgen.libngpy._meshing.Mesh(2)
    labels = None
labels = comm.bcast(labels, root=0) 
msh = Mesh(ngmesh)
# Setting up the problem
V = FunctionSpace(msh, "CG", 2)
u = TrialFunction(V)
v = TestFunction(V)
f = Function(V)
x, y = SpatialCoordinate(msh)
f.interpolate(2*sin(x)*sin(y))
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
f.interpolate(sin(x)*sin(y))
File("VTK/Poisson2D.pvd").write(u)

