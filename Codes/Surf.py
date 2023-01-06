from firedrake import *
from netgen.csg import *
from netgen.meshing import MeshingParameters
from netgen.meshing import MeshingStep
from irksome import GaussLegendre, Dt, TimeStepper
geo = CSGeometry()
geo.Add(Sphere(Pnt(0,0,0),1).bc("sphere"))
mp = MeshingParameters(maxh=0.1, perfstepsend = MeshingStep.MESHSURFACE)
ngmesh = geo.GenerateMesh(mp=mp)
msh = Mesh(ngmesh)
File("VTK/SurfMesh.pvd").write(msh)
outfile = File("VTK/Surf.pvd")
V = FunctionSpace(msh, 'Lagrange', 1)
p = Function(V, name="p")
phi = Function(V, name="phi")
u = TrialFunction(V)
v = TestFunction(V)
x, y,z = SpatialCoordinate(msh)
phi.interpolate(sin(10*z))
outfile.write(phi)
T = 10.
dt = 0.001
t=0.
step = 0
while t<T:
    print(t)
    t=t+dt
    phi -= dt / 2 * p
    solve(u * v * dx == v * p * dx + dt * inner(grad(v), grad(phi)) * dx, p,solver_parameters={'ksp_type': 'cg','pc_type': 'sor','pc_sor_symmetric':True})
    phi -= dt / 2 * p
    t += dt
    step = step +1
    if step % 10 == 0:
        outfile.write(phi, time=t)
"""
V = FunctionSpace(msh, "CG", 2)
u = TrialFunction(V)
v = TestFunction(V)
f = Function(V)
x, y,z = SpatialCoordinate(msh)
f.interpolate(sin(10*z))
a = inner(grad(u), grad(v))*dx
l = inner(f, v) * dx
u = Function(V)
# Assembling matrix
A = assemble(a)
b = assemble(l)
# Solving the problem
solve(A, u, b, solver_parameters={"ksp_type": "preonly", "pc_type": "lu"})
# Computing the error
File("VTK/Surf.pvd").write(f)
"""

