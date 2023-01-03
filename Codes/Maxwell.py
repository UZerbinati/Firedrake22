from firedrake import *
from netgen.csg import *

geom = CSGeometry()
cube1 = OrthoBrick(Pnt(-1,-1,-1), Pnt(1,1,1))
cube2 = OrthoBrick(Pnt(0,0,0), Pnt(2,2,2))
fichera = cube1-cube2
geom.Add (fichera)
geom.SingularEdge (cube2,cube2, 1)
geom.SingularPoint (cube2,cube2,cube2, 1)
ngmsh = geom.GenerateMesh(maxh=0.5)
msh = Mesh(ngmsh)
File("VTK/Maxwell.pvd").write(msh)
ngmsh.HPRefine(2,0.3)
msh = Mesh(ngmsh)
File("VTK/MaxwellRefined.pvd").write(msh)
