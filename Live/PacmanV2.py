from firedrake import *
from mpi4py import MPI
from netgen.geom2d import SplineGeometry
from netgen.occ import *
import netgen
import numpy as np
import pytest
from time import sleep

comm = MPI.COMM_WORLD
from petsc4py import PETSc
from slepc4py import SLEPc

print = lambda msg : PETSc.Sys.Print(msg)

def Solve(msh,labels):
    V = FunctionSpace(msh, "CG", 2)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(grad(u), grad(v))*dx
    m = (u*v)*dx
    uh = Function(V)
    bc = DirichletBC(V, 0, labels)
    A = assemble(a, bcs=bc)
    M = assemble(m)
    Asc, Msc = A.M.handle, M.M.handle
    E = SLEPc.EPS().create()
    E.setType(SLEPc.EPS.Type.ARNOLDI)
    E.setProblemType(SLEPc.EPS.ProblemType.GHEP)
    E.setDimensions(1, SLEPc.DECIDE)
    E.setOperators(Asc, Msc)
    ST = E.getST()
    ST.setType(SLEPc.ST.Type.SINVERT)
    PC = ST.getKSP().getPC()
    PC.setType("lu")
    PC.setFactorSolverType("mumps")
    E.setST(ST)
    E.solve()
    vr, vi = Asc.getVecs()
    with uh.dat.vec_wo as vr:
        lam = E.getEigenpair(0, vr, vi)
    return (lam, uh, V)

def Mark(msh, uh, lam):
    W = FunctionSpace(msh, "DG", 0)
    w = TestFunction(W)
    R_T = lam.real*uh + div(grad(uh))
    n = FacetNormal(V.mesh())
    h = CellDiameter(msh)
    R_dT = dot(grad(uh), n)
    eta = assemble(h**2*R_T**2*w*dx + (h("+")+h("-"))*(R_dT("+")-R_dT("-"))**2*(w("+")+w("-"))*dS)
    frac = .95
    delfrac = .05
    part = .2
    mark = Function(W)
    with mark.dat.vec as markedVec:
        with eta.dat.vec as etaVec:
            sum_eta = etaVec.sum()
            if sum_eta < tolerance:
                return markedVec
            eta_max = etaVec.max()[1]
            sct, etaVec0 = PETSc.Scatter.toZero(etaVec)
            markedVec0 = etaVec0.duplicate()
            sct(etaVec,etaVec0)
            if etaVec.getComm().getRank() == 0:
                eta = etaVec0.getArray()
                marked = np.zeros(eta.size, dtype='bool')
                sum_marked_eta = 0.
                while sum_marked_eta < part*sum_eta:
                    new_marked = (~marked) & (eta > frac*eta_max)
                    sum_marked_eta += sum(eta[new_marked])
                    marked += new_marked
                    frac -= delfrac
                markedVec0.getArray()[:] = 1.0*marked[:]
            sct(markedVec0,markedVec,mode=PETSc.Scatter.Mode.REVERSE)
    return mark

tolerance = 1e-16
max_iterations = 15
exact = 3.375610652693620492628**2
geo=SplineGeometry()
pnts=[(0, 0), (1, 0), (1, 1),(0, 1),
    (-1, 1), (-1, 0),(-1, -1), (0, -1)]
P = [geo.AppendPoint(*pnt) for pnt in pnts]
L,B=["line","spline3"],["line","curve"]
curves = [[[L[0],P[0],P[1]],B[0]],
          [[L[1],P[1],P[2],P[3]],B[1]],
          [[L[1],P[3],P[4],P[5]],B[1]],
          [[L[1],P[5],P[6],P[7]],B[1]],
          [[L[0],P[7],P[0]],B[0]]]
[geo.Append(c, bc=bc) for c, bc in curves]
if comm.rank == 0:
    ngmsh = geo.GenerateMesh(maxh=0.2)
    labels=sum([ngmsh.GetBCIDs(label) for label in ["line","curve"]], [])
else:
    ngmsh=netgen.libngpy._meshing.Mesh(2)
    labels = None
msh = Mesh(ngmsh)
labels = comm.bcast(labels, root=0)
for i in range(max_iterations):
    print("Refinement Level {}".format(i))
    lam, uh, V = Solve(msh,labels)
    mark = Mark(msh, uh, lam)
    msh = msh.Refine(mark)
    sleep(2)
    msh.topology_dm.viewFromOptions("-dm_viewer")
    File("VTK/PacManAdp.pvd").write(uh,mark)
assert(abs(lam-exact)<1e-2)

