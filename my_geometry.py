from ngsolve import *
from netgen.csg import *
from ngsolve.internal import *

def MakeGeometry(Rext,Rpml,Rminus,Rplus,c,Rother,res,hmax):
    geometry = CSGeometry()
    o_ext = (Sphere(Pnt(0,0,0), Rext)).bc("outer")

    pml = Sphere(Pnt(0,0,0),Rpml)

    o_plus = Sphere(Pnt(0,0,0), Rplus).bc("interface")
    o_minus = (Sphere(Pnt(0,0,0), Rminus)).maxh(res).mat("ominus")
    other = Sphere(Pnt(0,c,0),Rother).maxh(res)
    #other = Cylinder(Pnt(0,0,0),Pnt(0,1,0),0.2) * OrthoBrick(Pnt(-0.4,0,-0.4),Pnt(0.4,0.6,0.4))
    
    geometry.Add ((o_ext - pml).mat("air"))
    geometry.Add ((pml-o_plus).mat("air"))
    geometry.Add (((o_plus-o_minus)-other).mat("oplus"))
    geometry.Add ((other-o_minus).mat("olayer"))
    geometry.Add (o_minus)
    ####
    mesh = Mesh(geometry.GenerateMesh(maxh=hmax))
    mesh.Curve(5)
    nv=specialcf.normal(mesh.dim) # normals for the mesh
    ####
    viewoptions.clipping.enable = 1
    viewoptions.clipping.ny = 0
    viewoptions.clipping.nz = -1
    Draw(mesh)
    collection='interface'

    return mesh,nv,collection

def sphere_geom(R,hD,hair,PMLinnerrad,PMLouterrad):
    geometry = CSGeometry()
    PMLOuterSphere = Sphere(Pnt(0,0,0),PMLouterrad).bc('outer') # Outermost sphere of comp domain
    PMLInnerSphere = Sphere(Pnt(0,0,0),PMLinnerrad) # Middle (PML) sphere
    ##
    DomainD = Sphere(Pnt(0,0,0),R).maxh(hD).bc('dD') # D is a sphere here with labeled boundary
    ##
    geometry.Add(DomainD.mat('D')) # Set material of inner sphere to 'D'
    geometry.Add((PMLInnerSphere-DomainD).mat('air')) # Set material of domain between PML and D to 'air'
    geometry.Add((PMLOuterSphere-PMLInnerSphere).mat('pml')) # This is the PML labeled "pml"
    ####
    mesh = Mesh(geometry.GenerateMesh(maxh=hair))
    mesh.Curve(5)
    nv=specialcf.normal(mesh.dim) # normals for the mesh
    Draw(mesh)
    collection='dD'
    # There is no need to put collection on the boundary of D but
    # outside must be air
    return mesh,nv,collection

def sphere_stek_geom(R,hair,PMLinnerrad,PMLouterrad):
    geometry = CSGeometry()
    sphere = Sphere(Pnt(0,0,0),PMLouterrad).bc('outer')
    pmlSphere = Sphere(Pnt(0,0,0),PMLinnerrad)
    innerSphere = Sphere(Pnt(0,0,0),R).bc('bndinner')
    geometry.Add((pmlSphere-innerSphere).mat('air'))
    geometry.Add((sphere-pmlSphere).mat('PML'))
    netgenMesh = geometry.GenerateMesh(maxh=hair)
    mesh = Mesh(netgenMesh)
    mesh.Curve(5)
    Draw(mesh)
    nv=specialcf.normal(mesh.dim)
    collection='bndinner'
    return mesh,nv,collection

def cube_geom(L,hD,hair,PMLinnerrad,PMLouterrad):
    geometry = CSGeometry()
    PMLOuterSphere = Sphere(Pnt(0,0,0),PMLouterrad).bc('outer') # Outermost sphere of comp domain
    PMLInnerSphere = Sphere(Pnt(0,0,0),PMLinnerrad) # Middle (PML) sphere
    ##
    DomainD = OrthoBrick(Pnt(-L/2,-L/2,-L/2),Pnt(L/2,L/2,L/2)).maxh(hD).bc('dD')
    ##
    geometry.Add(DomainD.mat('D')) # Set material of inner sphere to 'D'
    geometry.Add((PMLInnerSphere-DomainD).mat('air')) # Set material of domain between PML and D to 'air'
    geometry.Add((PMLOuterSphere-PMLInnerSphere).mat('pml')) # This is the PML labeled "pml"
    ####
    mesh = Mesh(geometry.GenerateMesh(maxh=hair))
    mesh.Curve(5)
    nv=specialcf.normal(mesh.dim) # normals for the mesh
    Draw(mesh)
    collection='dD'
    # There is no need to put collection on the boundary of D but
    # outside must be air
    return mesh,nv,collection


def L_geom(L,hD,hair,PMLinnerrad,PMLouterrad):
    geometry = CSGeometry()
    PMLOuterSphere = Sphere(Pnt(0,0,0),PMLouterrad).bc('outer') # Outermost sphere of comp domain
    PMLInnerSphere = Sphere(Pnt(0,0,0),PMLinnerrad) # Middle (PML) sphere
    ##
    DomainD1 = OrthoBrick(Pnt(-L,-L,-L),Pnt(0,0,0)) 
    DomainD2 = OrthoBrick(Pnt(-L,-L,0),Pnt(0,0,L))
    DomainD3 = OrthoBrick(Pnt(-L,0,0),Pnt(0,L,L))
    DomainD = DomainD1+DomainD2+DomainD3
    DomainD.bc('dD')
    ##
    geometry.Add(DomainD.mat('D')) # Set material of inner sphere to 'D'
    geometry.Add((PMLInnerSphere-DomainD).mat('air')) # Set material of domain between PML and D to 'air'
    geometry.Add((PMLOuterSphere-PMLInnerSphere).mat('pml')) # This is the PML labeled "pml"
    ####
    mesh = Mesh(geometry.GenerateMesh(maxh=hair))
    mesh.Curve(5)
    nv=specialcf.normal(mesh.dim) # normals for the mesh
    Draw(mesh)
    collection='dD'
    # There is no need to put collection on the boundary of D but
    # outside must be air
    return mesh,nv,collection

