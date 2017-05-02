from ngsolve import *
from netgen.csg import *
from ngsolve.internal import *


def brick2_geometry(Rminus, Rplus, Rext, Rpml, delta, hsample, hmax):
    geometry = CSGeometry()
    o_minus = (OrthoBrick(Pnt(-Rminus,-Rminus,-Rminus),Pnt(Rminus,Rminus,Rminus)).maxh(0.1))
    
    geometry.Add (o_minus.mat("ominus").maxh(0.3))    
    #geometry.Add ((box * pl1 * pl2).mat("olayer").maxh(hmax),bcmod=[(pl1,"crack"),(box,"sides"),(pl2,"top")])

    #slices = [2**(-i) for i in reversed(range(1,6))]
    #geometry.CloseSurfaces(pl1,pl2)#,slices)
    
    return geometry
def brick_geometry(Rminus, Rplus, Rext, Rpml, delta, hsample, hmax):
    geometry = CSGeometry()
    o_ext = (Sphere(Pnt(0,0,0), Rext)).bc("outer")

    pml = Sphere(Pnt(0,0,0),Rpml)

    o_plus = Sphere(Pnt(0,0,0), Rplus).bc("interface")

    #This is to define the two close surfaces for the thin layer:
    box = OrthoBrick(Pnt(-Rminus,-Rminus,-Rminus),Pnt(Rminus,Rminus,Rminus))
    pl1 = Plane(Pnt(0,0,Rminus),Vec(0,0,-1)).bc("crack")
    pl2 = Plane(Pnt(0,0,Rminus+delta),Vec(0,0,1))#.bc("top")
    o_minus = (box - pl1).maxh(hsample)

    geometry.Add (o_minus.mat("ominus").maxh(hmax))    
    geometry.Add ((o_ext - pml).mat("pml"))
    geometry.Add ((pml-o_plus).mat("air"))
    geometry.Add ((o_plus-box).mat("oplus").maxh(hmax))
    geometry.Add ((box * pl1 * pl2).mat("olayer").maxh(hmax))

    #slices = [2**(-i) for i in reversed(range(1,6))]
    geometry.CloseSurfaces(pl1,pl2)#,slices)
    
    return geometry

def sphere_geometry(Rminus, Rplus, Rext, Rpml, Rother, c, delta, hsample, hmax):
    geometry = CSGeometry()
    o_ext = (Sphere(Pnt(0,0,0), Rext)).bc("outer")
    pml = Sphere(Pnt(0,0,0),Rpml)
    o_plus = Sphere(Pnt(0,0,0), Rplus).bc("interface")

    #This is to define the two close surfaces for the thin layer:
    
    o_minus = (Sphere(Pnt(0,0,0), Rminus)).maxh(hmax).mat("ominus").maxh(hsample)
    other = Sphere(Pnt(0,c,0),Rother)
    withCrack = (o_minus * other)
    withoutCrack = (o_minus - other)

    geometry.Add ((o_ext - pml).mat("air"))
    geometry.Add ((pml-o_plus).mat("air"))
    geometry.Add ((o_plus-o_minus-other).mat("oplus"))
    geometry.Add ((other-o_minus).mat("olayer"))
    geometry.Add (withCrack,bcmod=[(o_minus,"crack")])
    geometry.Add (withoutCrack,bcmod=[(o_minus,"nocrack")])

    return geometry
