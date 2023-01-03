from numpy import pi
from firedrake import *
from mpi4py import MPI
from netgen.geom2d import SplineGeometry
import netgen
comm = MPI.COMM_WORLD
if comm.rank == 0:
    geo = SplineGeometry()
    geo.AddRectangle((0, 0), (pi, pi), bc="rect")
    ngmesh = geo.GenerateMesh(maxh=0.1)
    labels = ngmesh.GetBCIDs("rect")
else:
    ngmesh = netgen.libngpy._meshing.Mesh(2)
    labels = None
labels = comm.bcast(labels, root=0) 
msh = Mesh(ngmesh)
File("VTK/Poisson2DMesh.pvd").write(msh)

