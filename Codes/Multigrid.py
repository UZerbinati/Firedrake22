from math import pi
from firedrake import *
from netgen.occ import *
cyl = Cylinder((0,0,0), Z, r=0.01, h=0.03).faces[0]
heli = Edge(Segment((0,0), (12*pi, 0.03)), cyl)
ps,vs = heli.start, heli.start_tangent
pe,ve = heli.end, heli.end_tangent
e1 = Segment((0,0,-0.03), (0,0,-0.01))
c1 = BezierCurve( [(0,0,-0.01), (0,0,0), ps-vs, ps])
e2 = Segment((0,0,0.04), (0,0,0.06))
c2 = BezierCurve( [pe, pe+ve, (0,0,0.03), (0,0,0.04)])
spiral = Wire([e1, c1, heli, c2, e2])
circ = Face(Wire([Circle((0,0,-0.03), Z, 0.001)]))
coil = Pipe(spiral, circ)
coil.faces.maxh, coil.faces.name = 0.1, "coilbnd"
coil.faces.Max(Z).name,coil.faces.Min(Z).name="I","O"
ngmsh = OCCGeometry(coil).GenerateMesh(maxh=0.3)
msh = Mesh(ngmsh)
hierarchy = MeshHierarchy(msh, 2)
V = FunctionSpace(hierarchy[-1], "CG", 1)
u,v  = TrialFunction(V), TestFunction(V)
a,L = dot(grad(u), grad(v))*dx, 1*v*dx
bcsI=DirichletBC(V,1,ngmsh.GetBCIDs("I"))
bcsO=DirichletBC(V,0.,ngmsh.GetBCIDs("O"))
u = Function(V)
parameters = {"ksp_type": "preonly", "pc_type": "mg",
   "pc_mg_type": "full", "mg_levels_ksp_type": "chebyshev",
   "mg_levels_ksp_max_it": 2,"mg_levels_pc_type": "jacobi"}
solve(a==L, u, bcs=[bcsI, bcsO],solver_parameters=par)
File("VTK/Coil.pvd").write(u)
