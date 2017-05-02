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
from ffcode import ExportFFMesh, SetupFFMesh, findffvec
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
#ngmesh = MakeGeometry().GenerateMesh(maxh=0.3)
##ngmesh.Save("scatterer.vol")
#mesh = Mesh(ngmesh)

ngsglobals.testout = "test.out"
# First, set PML parameters for the ATC model
SetPMLParameters(rad=Rpml,alpha=.6)
# curve elements for geometry approximation
mesh.Curve(4)

ngsglobals.msg_level = 5

Vext = HCurl(mesh, order=2, complex=True, definedon=mesh.Materials("air")+mesh.Materials("pml"), dirichlet="outer")
Vplus  = HCurl(mesh, order=2, complex=True, definedon=mesh.Materials("oplus")+mesh.Materials("olayer"))
Vminus = HCurl(mesh, order=2, complex=True, definedon=mesh.Materials("ominus"))
print(mesh.GetBoundaries())

defbound = [i+1 for i, bc in enumerate(mesh.GetBoundaries()) if bc == "crack"]
Vcrack = HCurl(mesh, order = 1,complex=True, definedon=[],flags={"definedonbound":defbound})
#defbound2 = [i+1 for i, bc in enumerate(mesh.GetBoundaries()) if bc == "nocrack"]
#Vnocrack = HCurl(mesh, order = 1,complex=True, definedon=[],flags={"definedonbound":defbound2})
#Vplus_l = HCurl(mesh, order=2, complex=True, definedon=mesh.Materials("oplus")+mesh.Materials("olayer")+mesh.Materials("ominus"))

V=FESpace([Vext,Vplus,Vminus,Vcrack])
#V_l=FESpace([Vext,Vplus_l])
V_l=FESpace([Vext,Vplus,Vminus])
#print("V:",V)
# u and v refer to trial and test-functions in the definition of forms below
uext,uplus,uminus,ubnd = V.TrialFunction()
vext,vplus,vminus,vbnd = V.TestFunction()

#uext_l,uplus_l = V_l.TrialFunction()
#vext_l,vplus_l = V_l.TestFunction()

uext_l,uplus_l,uminus_l= V_l.TrialFunction()
vext_l,vplus_l,vminus_l = V_l.TestFunction()

#Material properties for the crack
mur = {"ominus" : mu1 , "oplus" : mu2 , "olayer" : mu2 , "air" : 1, "pml" : 1 }
nu_coef = [ 1/(mur[mat]) for mat in mesh.GetMaterials() ]
nu = CoefficientFunction ([ 1/(mur[mat]) for mat in mesh.GetMaterials() ])

epsilonr = {"ominus" : eps1 , "oplus" : eps2 , "olayer" : eps2 , "air" : 1, "pml" : 1 }
epsilon_coef = [ epsilonr[mat] for mat in mesh.GetMaterials() ]
epsilon = CoefficientFunction (epsilon_coef)

#Material properties for the layer
mur_l = {"ominus" : mu1 , "oplus" : mu2 , "olayer" : mu0 , "air" : 1 , "pml" : 1}
nu_l_coef = [ 1/(mur_l[mat]) for mat in mesh.GetMaterials() ]
nu_l = CoefficientFunction ([ 1/(mur_l[mat]) for mat in mesh.GetMaterials() ])

epsilonr_l = {"ominus" : eps1 , "oplus" : eps2 , "olayer" : eps0 , "air" : 1 , "pml" : 1}
epsilon_coef_l = [ epsilonr_l[mat] for mat in mesh.GetMaterials() ]
epsilon_l = CoefficientFunction (epsilon_coef_l)

nv = specialcf.normal(mesh.dim) # normal vector
Cross = lambda u,v: CoefficientFunction((u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0])) # cross product


#Incident field
Ein= exp(k*1J*(d[0]*x+d[1]*y+d[2]*z))*CoefficientFunction(p)
Einext=GridFunction(Vext,"gEin")
Einext.Set(Ein)
curlEin = k*1.J*Cross(d,p)*exp(k*1.J*(d[0]*x+d[1]*y+d[2]*z))

# Coefficients for the asymptotic model
alpha1 = 2*mu0
alpha2 = 2*k*k*eps0
beta1  = 2*1./(k*k*eps0)
beta2  = 2*1./mu0

#Durufle model
#alpha1 = (mu0-mu2)
#alpha2 = (k*k*eps0-k*k*eps2)
#beta1  = (1./(k*k*eps0)-1./(k*k*eps2))
#beta2  = (1./mu0-1./mu2)
################################################################################################################
#     Sesquilinear form for the second order ATC's model
# Here, outside the obstacle (in o_ext, we compute the scattered field, while in the obstacle (o_plus + o_minus)
# we compute the total field.
################################################################################################################

a = BilinearForm(V, symmetric=True, flags={"printelmat":True})

a.components[0] += BFI("PML_curlcurledge",coef=nu)
a.components[0] += BFI("PML_massedge",coef=-k*k*epsilon)

a += SymbolicBFI(nu*curl(uplus)*curl(vplus)   - k*k*epsilon*uplus*vplus)
a += SymbolicBFI(nu*curl(uminus)*curl(vminus) - k*k*epsilon*uminus*vminus)

# Nietsche's method for the transmission between the scattered field in the exterior and the total field inside the obstacle
a += SymbolicBFI(0.5*( (1./mu2)*curl(uplus).Trace() + curl(uext).Trace()) * Cross(nv,vext.Trace()-vplus.Trace() ) ,BND,definedon=mesh.Boundaries("interface"))
a += SymbolicBFI( Cross(nv,uext.Trace()-uplus.Trace()) * 0.5*( (1./mu2)*curl(vplus).Trace() + curl(vext).Trace() ) ,BND,definedon=mesh.Boundaries("interface"))
a += SymbolicBFI( gamma/hmax*1.J*k* (uext.Trace() - uplus.Trace()) * (vext.Trace() - vplus.Trace()),BND,definedon=mesh.Boundaries("interface"))

# Second order asymptotic model using ATC's to simulate the effect of a thin layer betweeen o_plus and o_minus
# NOTE (03/09/2017): needs to be updated, to account for the curvature terms popping out from the determinant. Otherwise, use the
# contrast model (although no well posedness result for the latter model has been proved).
a += SymbolicBFI( -delta*0.25* alpha2*( uplus.Trace()+uminus.Trace() )*( vplus.Trace()+vminus.Trace() ),BND,definedon=mesh.Boundaries("crack"))
a += SymbolicBFI(  delta*0.25* beta2 *( uplus.Trace().Deriv()+uminus.Trace().Deriv() )*( vplus.Trace().Deriv()+vminus.Trace().Deriv() ),BND,definedon=mesh.Boundaries("crack"))
a += SymbolicBFI(  ubnd.Trace() * Cross(nv,vplus.Trace()-vminus.Trace())  ,BND,definedon=mesh.Boundaries("crack"))
a += SymbolicBFI(  Conj(Cross(nv,uplus.Trace()-uminus.Trace())) * Conj(vbnd.Trace()) ,BND,definedon=mesh.Boundaries("crack"))
a += SymbolicBFI(-delta*alpha1*ubnd.Trace()*vbnd.Trace() + delta*beta1*( ubnd.Trace().Deriv() * vbnd.Trace().Deriv() ),BND,definedon=mesh.Boundaries("crack"))

# Nitsche's term in nocrack
a += SymbolicBFI( gamma/hmax*1.J*k* (uplus.Trace() - uminus.Trace()) * (vplus.Trace() - vminus.Trace()),BND,definedon=mesh.Boundaries("nocrack"))
# Nitsche's term in sides
#a += SymbolicBFI( gamma/hmax*1.J*k* uplus.Trace() * vplus.Trace(),BND,definedon=mesh.Boundaries("sides"))
#a += SymbolicBFI( gamma/hmax*1.J*k* ((1./mu2)*curl(uplus).Trace() - (1./mu1)*curl(uminus).Trace()) * ((1./mu2)*curl(vplus).Trace() - (1./mu1)*curl(vminus).Trace()),BND,definedon=mesh.Boundaries("nocrack"))
#a += SymbolicBFI(  ubnd2.Trace() * Cross(nv,vplus.Trace()-vminus.Trace())  ,BND,definedon=mesh.Boundaries("nocrack"))
#a += SymbolicBFI(  Conj(Cross(nv,uplus.Trace()-uminus.Trace())) * Conj(vbnd2.Trace()) ,BND,definedon=mesh.Boundaries("nocrack"))
#a += SymbolicBFI(-1.e-10*ubnd2.Trace()*vbnd2.Trace() + 1.e-10*( ubnd2.Trace().Deriv() * vbnd2.Trace().Deriv() ),BND,definedon=mesh.Boundaries("nocrack"))
#RHS for the crack
f = LinearForm(V)
f += SymbolicLFI( curlEin * 0.5 * (Cross(nv,vplus.Trace() + vext.Trace() )),BND,definedon=mesh.Boundaries("interface"))
f += SymbolicLFI(Ein * 0.5*Cross(nv,(1./mu2)*curl(vplus).Trace()+curl(vext).Trace() ),BND,definedon=mesh.Boundaries("interface"))
f += SymbolicLFI(-gamma/hmax*1.J*k* Ein * (vext.Trace()-vplus.Trace()),BND,definedon=mesh.Boundaries("interface"))

# Set PML parameters for the "exact" model
SetPMLParameters(rad=Rpml,alpha=0.6)

################################################################################################################
#     Sesquilinear form for the exact layer model
# Outside the obstacle (in o_ext) we compute the scattered field, while inside the obstacle (o_plus + o_minus + o_layer)
# we compute the total field.
################################################################################################################
a_l = BilinearForm(V_l, symmetric=True, flags={"printelmat":True})

a_l.components[0] += BFI("PML_curlcurledge",coef=nu_l)
a_l.components[0] += BFI("PML_massedge",coef=-k*k*epsilon_l)

a_l += SymbolicBFI(nu_l*curl(uplus_l)*curl(vplus_l)   - k*k*epsilon_l*uplus_l*vplus_l)
a_l += SymbolicBFI(nu_l*curl(uminus_l)*curl(vminus_l)   - k*k*epsilon_l*uminus_l*vminus_l)

# Nietsche's method for the transmission between the scattered field in the exterior and the total field inside the obstacle
a_l += SymbolicBFI(0.5*( (1./mu2)*curl(uplus_l).Trace() + curl(uext_l).Trace()) * Cross(nv,vext_l.Trace()-vplus_l.Trace() ) ,BND,definedon=mesh.Boundaries("interface"))
a_l += SymbolicBFI( Cross(nv,uext_l.Trace()-uplus_l.Trace()) * 0.5*( (1./mu2)*curl(vplus_l).Trace() + curl(vext_l).Trace() ) ,BND,definedon=mesh.Boundaries("interface"))
a_l += SymbolicBFI( gamma/hmax*1.J*k* (uext_l.Trace() - uplus_l.Trace()) * (vext_l.Trace() - vplus_l.Trace()),BND,definedon=mesh.Boundaries("interface"))

# Nitsche's term in nocrack
a_l += SymbolicBFI( gamma/hmax*1.J*k* (uplus_l.Trace() - uminus_l.Trace()) * (vplus_l.Trace() - vminus_l.Trace()),BND,definedon=mesh.Boundaries("nocrack"))
a_l += SymbolicBFI( gamma/hmax*1.J*k* (uplus_l.Trace() - uminus_l.Trace()) * (vplus_l.Trace() - vminus_l.Trace()),BND,definedon=mesh.Boundaries("crack"))

#RHS for the layer
f_l = LinearForm(V_l)
f_l += SymbolicLFI( curlEin * 0.5 * (Cross(nv,vplus_l.Trace() + vext_l.Trace() )),BND,definedon=mesh.Boundaries("interface"))
f_l += SymbolicLFI(Ein * 0.5*Cross(nv,(1./mu2)*curl(vplus_l).Trace()+curl(vext_l).Trace() ),BND,definedon=mesh.Boundaries("interface"))
f_l += SymbolicLFI(-gamma/hmax*1.J*k* Ein * (vext_l.Trace()-vplus_l.Trace()),BND,definedon=mesh.Boundaries("interface"))

#Define solution functions
u_l = GridFunction(V_l)#Layer Solution
u = GridFunction(V) #ATC Solution

f.Assemble()
a.Assemble()

f_l.Assemble()
a_l.Assemble()

u.vec.data += a.mat.Inverse(V.FreeDofs()) * f.vec
u_l.vec.data += a_l.mat.Inverse(V_l.FreeDofs()) * f_l.vec


viewoptions.clipping.enable = 1
viewoptions.clipping.nx = -1
#viewoptions.clipping.ny = 0
#viewoptions.clipping.nz = 0
visoptions.clipsolution = "scal"

nopml = CoefficientFunction ([0 if mat=='pml' else 1 for mat in mesh.GetMaterials() ])

#Draw(Einext,mesh,"Einext")
Draw(CoefficientFunction(u.components[2])+CoefficientFunction(u.components[1])+nopml*CoefficientFunction(u.components[0])+nopml*Einext,mesh,"E")
Draw(u.components[0],mesh,"uext")
#Draw(u.components[1],mesh,"uplus")
#Draw(u.components[2],mesh,"uminus")

Draw(CoefficientFunction(u_l.components[2])+CoefficientFunction(u_l.components[1])+nopml*CoefficientFunction(u_l.components[0])+nopml*Einext,mesh,"E_l")
Draw(u_l.components[0],mesh,"uext_l")
#Draw(u_l.components[1],mesh,"uplus_l")

#solver = CGSolver(mat=a.mat, pre=c.mat)
#u.vec.data = a.mat.Inverse() * f.vec
#solver = CGSolver(mat=a.mat, pre=c.mat)

# Error
Error = (u.components[0]-u.components[0]) + (u.components[1] -u_l.components[1]) + (u_l.components[2]-u.components[2])
Draw(Error,mesh,"Error")

# Old error
#error = nopml*(u.components[0]-u_l.components[0])
#uexterior = nopml*u.components[0]

uexact = nopml*u_l.components[0] + u_l.components[1] + u_l.components[2]

L2error = Integrate((Error*Conj(Error)), mesh, VOL, 2*order)

L2normu = Integrate((uexact*Conj(uexact)), mesh, VOL, 2*order)

L2rel_error = sqrt(L2error)/sqrt(L2normu)

print('hmax =',hmax)
print('delta =',delta)
print('L2error =',L2error)
print('L2normu =',L2normu)
print('L2rel_error =',L2rel_error)

# Hcurl Error
#curlerror=GridFunction(Vext,'curlerr')
#curlerror.Set(nopml*curl(u.components[0])-nopml*curl(u_l.components[0]))
#curluexterior=GridFunction(Vext,'curluext')
#curluexterior.Set(nopml*curl(u.components[0]))

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

