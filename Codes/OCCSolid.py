from firedrake import *
from netgen.occ import *
box = Box((0,0,0), (3,0.6,1))
cyl = sum( [Cylinder((0.5+i,0,0.5), Y, 0.25,0.8) for i in range(3)] )
box.faces.name,cyl.faces.name="outer","cyl"
geo = box-cyl
cylboxedges = geo.faces["outer"].edges * geo.faces["cyl"].edges
cylboxedges.name = "cylbox"
geo = geo.MakeChamfer(cylboxedges, 0.03)
geo.faces.Min(X).name = "fix"
geo.faces.Max(X).name = "force"
ngmesh = OCCGeometry(geo).GenerateMesh(maxh=0.1)
msh = Mesh(ngmesh)
File("VTK/OCCSolidMesh.pvd").write(msh)
V = VectorFunctionSpace(msh, "CG", 3)
W = FunctionSpace(msh, "CG", 3)
u,v = TrialFunction(V), TestFunction(V)
sol = Function(V)
lbFix = ngmesh.GetBCIDs("fix")
lbF = ngmesh.GetBCIDs("force")[0]
bc = DirichletBC(V, 0, lbFix)
f = as_vector([1e-3, 0, 0])
E, nu, Id = 210.,0.3, Identity(3)
mu = Constant((0.5/(1+nu))*E)
lambda_ = Constant(E*(nu/((1+nu)*(1-2*nu))))
epsilon = lambda u: 0.5*(grad(u) + grad(u).T)
sigma = lambda u: lambda_*div(u)*Id + 2*mu*epsilon(u)
a = inner(sigma(u), epsilon(v))*dx
L = inner(f, v)*ds(lbF)
uh = Function(V, name="Displacement")
sp = {"ksp_type": "preonly",
      "pc_type": "cholesky",
      "pc_factor_mat_solver_type": "mumps"}
solve(a == L, uh, bcs=bc, solver_parameters=sp)
strain = Function(W).interpolate(inner(sigma(uh),sigma(uh))**(1/2))
File("VTK/OCCSolid.pvd").write(uh,strain)
