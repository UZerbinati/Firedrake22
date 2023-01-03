from firedrake import *
from netgen.csg import *
geo = CSGeometry()
brick = OrthoBrick(Pnt(-1,-1,-1),Pnt(1,1,1))
sphere = Sphere(Pnt(0.5,0.5,0.5),1)
sphere.maxh(0.05)
geo.Add(brick-sphere)
ngmesh = geo.GenerateMesh(maxh=0.3)
msh = Mesh(ngmesh)
File("VTK/CSG3DMesh.pvd").write(msh)

