from firedrake import *
from netgen.occ import *
H, W, T  = 70.,50.,30.
P1,P2=Pnt(-W/2,0,0),Pnt(-W/2,-T/4,0)
P3,P4=Pnt(0,-T/2.,0),Pnt(W/2.,-T/ 4.,0)
P5=Pnt(W/2.,0, 0)
S1=Segment(P1, P2)
A=ArcOfCircle(P2, P3, P4)
S2=Segment(P4, P5)
wire=Wire([S1, A, S2])
mwire=wire.Mirror(Axis((0,0,0),X))
wire=Wire([wire,mwire])
base=Face(wire).Extrude(0.3*Z)
msh=Mesh(OCCGeometry(base).GenerateMesh(maxh=0.5))
File("VTK/OCCBottleBase.pvd").write(msh)
f = Face(wire)
f.bc("bottom")
body=f.Extrude(H*Z)
body=body.MakeFillet(body.edges,T/12)
neckAx=Axes(body.faces.Max(Z).center,Z)
NeckR,NeckH=T/4,H/10
neck=Cylinder(neckAx,NeckR,NeckH)
body=body+neck
fmax=body.faces.Max(Z)
thickbody=body.MakeThickSolid([fmax],-T/50,1.e-3)
msh=Mesh(OCCGeometry(thickbody).GenerateMesh(maxh=3.))
File("VTK/OCCBottleBody.pvd").write(msh)
from math import pi
cyl1=Cylinder(neckAx, NeckR*0.99,1).faces[0]
cyl2=Cylinder(neckAx, NeckR*1.05,1).faces[0]
aPnt=Pnt(2*pi,NeckH/2.0)
aDir=Dir(2*pi,NeckH/4.0)
anAx2d=gp_Ax2d(aPnt, aDir)
aMajor,aMinor=2*pi,NeckH/10
arc1=Ellipse(anAx2d,aMajor,aMinor).Trim(0,pi)
arc2=Ellipse(anAx2d,aMajor,aMinor/4).Trim(0,pi)
seg=Segment(arc1.start, arc1.end)
wire1=Wire([Edge(arc1,cyl1),Edge(seg,cyl1)])
wire2=Wire([Edge(arc2,cyl2),Edge(seg,cyl2)])
threading=ThruSections([wire1,wire2])
bottle=thickbody+threading
ngmsh = OCCGeometry(bottle).GenerateMesh(maxh=0.05)
msh = Mesh(ngmsh)
print("Meshed")
File("VTK/OCCBottle.pvd").write(msh)
V = VectorFunctionSpace(msh, "CG", 2)
W = FunctionSpace(msh, "CG", 2)
u = TrialFunction(V)
v = TestFunction(V)
sol = Function(V)
labels = ngmsh.GetBCIDs("bottom")
bc = DirichletBC(V, 0, labels)
f = as_vector([0, 0, -2710*9.8])
mu = Constant(25.95e9)
lambda_ = Constant(55.27e9)
Id = Identity(3)
epsilon = lambda u: 0.5*(grad(u) + grad(u).T)
sigma = lambda u: lambda_*div(u)*Id + 2*mu*epsilon(u)
a = inner(sigma(u), epsilon(v))*dx
L = inner(f, v)*dx
uh = Function(V, name="Displacement")
sp = {"ksp_type": "preonly",
      "pc_type": "cholesky",
      "pc_factor_mat_solver_type": "mumps"}
solve(a == L, uh, bcs=bc, solver_parameters=sp)
strein = Function(W).interpolate(inner(epsilon(uh),epsilon(uh))**(1/2))
File("VTK/OCCStrein.pvd").write(strein)
