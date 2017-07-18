################################################################################################################
# (Version: april 2017)
# Program that computes the L2 and Hcurl error between the 2nd order ATC model vs. the exact model to simulate
# the effect of a thin layer of constant thickness.
# Outside the obstacle the scattered field is computed, while inside the obstacle the total field is computed.
# To this end, a Nitsche's method was used.
# In the outer domain a circular PML (implemented by Christopher LACKNER) was used.
################################################################################################################

from netgen.csg import *
from ngsolve import *
from ngsolve.internal import *
from netgen.meshing import SetTestoutFile
from mygeometry_1 import *
from aux_functions import *
import os

## Load Physical, Geometrical, and Computational parameters
from parameters import *

################################################################################################################
#      Definition of the geometry
################################################################################################################
print('Meshing parameter',2*3.141/k/hmax*order)

if geom == 1:
    geometry = ATCerror_brick_geometry(Rminus, Rplus, Rext, Rpml, delta, hmax)
    ngmesh = geometry.GenerateMesh()
if geom == 4:
    geometry = ATCerror_halfsphere_geometry(Rminus, Rplus, Rext, Rpml, delta, hmax)
    ngmesh = geometry.GenerateMesh()
if geom == 2:
    def MakeGeometry():
        geometry = CSGeometry()
        o_ext = (Sphere(Pnt(0,0,0), Rext)).bc("outer")

        pml = Sphere(Pnt(0,0,0),Rpml)

        o_plus = Sphere(Pnt(0,0,0), Rplus).bc("interface")

        #This is to define the two close surfaces for the thin layer:
    
        o_minus = (Sphere(Pnt(0,0,0), Rminus)).maxh(res).mat("ominus")
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
    ngmesh = MakeGeometry().GenerateMesh(maxh=hmax)
    # 4 = number of boundary conditions:
    for i in range(6):
        print("BC: ",ngmesh.GetBCName(i))
        if ngmesh.GetBCName(i) == "crack":
            crackIndex = i+1

    #the crack should go into domain o_plus, there are 3 materials defined
    for i in range(1,6):
        print("Material: ",ngmesh.GetMaterial(i))
        if ngmesh.GetMaterial(i) == "oplus":
            crackDomIndex = i+1
    ngmesh = MakeGeometry().GenerateMesh(maxh=hmax)
    thickness = 0.01
        
    # create boundary layer with: surface index, thickness, index of domain into which the layer will go, new material name of layer
    ngmesh.BoundaryLayer(crackIndex,thickness,crackDomIndex,"olayer")


#ZRefinement(ngmesh,MakeGeometry())
mesh = Mesh(ngmesh)
print(mesh.GetBoundaries())
Draw(mesh)

#ngsglobals.testout = "test.out"
#ngsglobals.msg_level = 5

# First, set PML parameters for the ATC model
SetPMLParameters(rad=Rpml,alpha=.6)
# curve elements for geometry approximation
mesh.Curve(1)

Vext = HCurl(mesh, order=2, complex=True, definedon=mesh.Materials("air")+mesh.Materials("pml"), dirichlet="outer")
Vplus  = HCurl(mesh, order=2, complex=True, definedon=mesh.Materials("oplus")+mesh.Materials("olayer"))
Vminus = HCurl(mesh, order=2, complex=True, definedon=mesh.Materials("ominus"))
print(mesh.GetBoundaries())

defbound = [i+1 for i, bc in enumerate(mesh.GetBoundaries()) if bc == "crack"]
Vcrack = HCurl(mesh, order = 1,complex=True, definedon=[],flags={"definedonbound":defbound})

V=FESpace([Vext,Vplus,Vminus,Vcrack])
V_l=FESpace([Vext,Vplus,Vminus])
###############################################################################################################
#     Incident field
###############################################################################################################
Ein= exp(k*1J*(d[0]*x+d[1]*y+d[2]*z))*CoefficientFunction(p)
Einext=GridFunction(Vext,"gEin")
Einext.Set(Ein)
curlEin = k*1.J*Cross(d,p)*exp(k*1.J*(d[0]*x+d[1]*y+d[2]*z))

################################################################################################################
#     Sesquilinear forms & assemblement of the FEM matrix using a Nitsche's method on the 
#     exterior boundary of the inhomogeneity ("interface").
################################################################################################################

# Material properties:
nu, epsilon = materials(mu1,mu2,mu2,eps1,eps2,eps2,mesh) # material functions for the ATC model
nu_l, epsilon_l = materials(mu1,mu2,mu0,eps1,eps2,eps0,mesh) # material functions for the full

nv = specialcf.normal(mesh.dim) # normal vector
Cross = lambda u,v: CoefficientFunction((u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0])) # cross product

# Assemblement of FEM matrices and sources
a,f = ATC_model_sesq_form(alpha_pml,Rpml,nu,epsilon,k,mu2,delta,nv,mesh,V,gamma,hmax,eps0,mu0,Ein,curlEin)
a_l,f_l = full_model_sesq_form(alpha_pml,Rpml,nu_l,epsilon_l,k,mu2,nv,mesh,V_l,gamma,hmax,Ein,curlEin)

#Define solution functions
u_l = GridFunction(V_l)#Layer Solution
u = GridFunction(V) #ATC Solution

u.vec.data += a.mat.Inverse(V.FreeDofs()) * f.vec
u_l.vec.data += a_l.mat.Inverse(V_l.FreeDofs()) * f_l.vec


viewoptions.clipping.enable = 1
viewoptions.clipping.nx = -1
visoptions.clipsolution = "scal"

nopml = CoefficientFunction ([0 if mat=='pml' else 1 for mat in mesh.GetMaterials() ])

Draw(CoefficientFunction(u.components[2])+CoefficientFunction(u.components[1])+nopml*CoefficientFunction(u.components[0])+nopml*Einext,mesh,"E")
Draw(u.components[0],mesh,"uext")

Draw(CoefficientFunction(u_l.components[2])+CoefficientFunction(u_l.components[1])+nopml*CoefficientFunction(u_l.components[0])+nopml*Einext,mesh,"E_l")
Draw(u_l.components[0],mesh,"uext_l")

# Error
Error = (u.components[0]-u.components[0]) + (u.components[1] -u_l.components[1]) + (u_l.components[2]-u.components[2])
Draw(Error,mesh,"Error")


uexact = nopml*u_l.components[0] + u_l.components[1] + u_l.components[2]

L2error = Integrate((Error*Conj(Error)), mesh, VOL, 2*order)
L2normu = Integrate((uexact*Conj(uexact)), mesh, VOL, 2*order)
L2rel_error = sqrt(L2error)/sqrt(L2normu)

print('hmax =',hmax)
print('delta =',delta)
print('L2error =',L2error)
print('L2normu =',L2normu)
print('L2rel_error =',L2rel_error)


curlerror0=GridFunction(Vext,'curlerr0')
curlerror1=GridFunction(Vplus,'curlerr1')
curlerror2=GridFunction(Vminus,'curlerr2')
curlerror0.Set(nopml*curl(u.components[0])-nopml*curl(u_l.components[0]))
curlerror1.Set(curl(u.components[1]) -curl(u_l.components[1])) 
curlerror2.Set(curl(u_l.components[2])-curl(u.components[2]))
curlerror=curlerror0+curlerror1+curlerror2

curluexact0=GridFunction(Vext,'curluexact0')
curluexact1=GridFunction(Vplus,'curluexact1')
curluexact2=GridFunction(Vminus,'curluexact2')

curluexact0.Set(nopml*curl(u_l.components[0]))
curluexact1.Set(curl(u_l.components[1]))
curluexact2.Set(curl(u_l.components[2]))
curluexact=curluexact0+curluexact1+curluexact2

L2curlerror=Integrate((curlerror*Conj(curlerror)),mesh,VOL,2*order)

L2curlnormu = Integrate((curluexact*Conj(curluexact)), mesh, VOL, 2*order)

Hcurl_rel_error = sqrt(L2error + L2curlerror)/sqrt(L2normu+L2curlnormu)

L2curlerror=sqrt(L2curlerror)
print('L2curlerror',L2curlerror)
L2curlnormu = sqrt(L2curlnormu)
print('L2curlnormu',L2curlnormu)
print('Hcurl_rel_error =',Hcurl_rel_error)
 
Draw(Error,mesh,"error")

Redraw()

