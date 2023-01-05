from firedrake import *
from netgen.occ import *
box = Box(Pnt(0,0,0), Pnt(1,1,1))
cyl = Cylinder(Pnt(1,0.5,0.5), X, r=0.3, h=0.5)
solid1 = box + cyl
solid2 = solid1.Rotate(Axis((0,0,0),Y),180).Move((2.5,0.,1.))
solid = solid2 + solid1
geo = OCCGeometry(solid)
ngmesh = geo.GenerateMesh(maxh=0.1)
msh = Mesh(ngmesh)
File("VTK/OCC.pvd").write(msh)

