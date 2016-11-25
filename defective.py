from netgen.csg import *
from ngsolve import *
from ngsolve.internal import *
from netgen.meshing import SetTestoutFile

#### FUNCTIONS
############################ cross products
Ocross=lambda a,b: (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-b[0]*a[1])
cross=lambda a,b: CoefficientFunction((a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-b[0]*a[1]))
############################ Far field
pi=4.*atan(1.)
def findff6(ii,jj,ka,ll,Es,order,mesh,nv,dv,pv,collection,nang,nphi,fespace,fff):

    # Use a finite element cutoff function and linear forms

    # This version outputs the ff pattern as a function of azimuthal angle

    # for comparison to MatScat (two calls are needed)

    #print('Using findff6....')

    w=fespace.TestFunction()

    zer0=CoefficientFunction((0,0,0))

    vv=GridFunction(fespace)

    vv.Set(zer0,VOL)

    for kk in range(0,nphi):

        phi=pi*kk/(nphi-1)

        FF=[]

        ang=[] #theta

        print('Starting kk=', kk) 

        for na in range(0,nang):

            ff=[0,0,0]

            theta=2*pi*na/nang

            xhat=(sin(phi)*cos(theta),sin(phi)*sin(theta),cos(phi))

            CEout = exp(1J*ka*(x*xhat[0]+y*xhat[1]+z*xhat[2])) # Conjugate to allow for

                                                          # conjugation in linear form
            for j in range(0,3):

                e=[0.,0.,0.]

                e[j]=1.

                ## Functional for the full FF pattern

                func=cross(Ocross(e,xhat),nv)*CEout*(-1J*ka)

                gf = LinearForm(fespace)

                gf += SymbolicLFI(func*w.Trace(),BND,definedon=mesh.Boundaries(collection))

                ## Functional for the second part using finite element cutoff

                eT=cross(xhat,Ocross(e,xhat))

                vv.Set(eT*CEout,BND,definedon=mesh.Boundaries(collection))

                gf += SymbolicLFI(curl(vv)*curl(w)-ka*ka*vv*w,VOL,definedon=mesh.Materials('air'))

                with TaskManager():

                    gf.Assemble()

                # Conjugates are needed since inner product conjugates the second component(?)

                ff1a=gf.vec.InnerProduct(Es.vec)

                ff[j]=ff1a.conjugate()/4./pi

            FF.append(ff)

            ang.append(theta)

        print('Done ',na,' of ',nang-1)

        for mm in range(0,nang):

            F1=FF[mm][0]

            F2=FF[mm][1]

            F3=FF[mm][2]

            theta=ang[mm]

            print("%d %d %d %d %12.12e %12.12e %12.12e %12.12e %12.12e %12.12e"%(ii,jj,kk,mm,F1.real,F1.imag,F2.real,F2.imag,F3.real,F3.imag),file=fff)

    return()

############################ Writing the background solution for the rhs

def writebackground(ii,jj,E,nang,nphi,rhs):
    for kk in range(0,nphi):

        phi=pi*kk/(nphi-1)

        print('Starting kk=', kk) 

        for mm in range(0,nang):

            theta=2*pi*mm/nang

            v=E(sin(phi)*cos(theta),sin(phi)*sin(theta),cos(phi))

            print("%d %d %d %d %12.12e %12.12e %12.12e %12.12e %12.12e %12.12e"%(ii,jj,kk,mm,v[0].real,v[0].imag,v[1].real,v[1].imag,v[2].real,v[2].imag),file=rhs)

    return()

############################ Solver for the construction of the FF operator

#Geometrical parameters
pi=4.*atan(1.)
k = 3. # Wave number
Rext = 2.7 #Radius of external circle
Rpml = 2. #Radius of PML
Rplus =  1.5 #Radius of Oplus
Rminus = 1. #Radius of Ominus
delta = 0.2 #Layer thickness
ddelta =delta/10 #parameter to avoid edge of the delamination 
#Geometry of the circular delamination
c = (2*delta+delta*delta)/(2+2*delta-sqrt(3)) #center of other
Rother = 1-c+delta # Radius of other
#Rotherr = Rother + ddelta #to avoid edge of delamination

#Physical parameters
mu2 = 1 #mu in oplus
mu1 = 1 #mu in ominus
mu0 = 1 #mu in the layer

eps2 = 1 #1.  + 1J*1.01 #epsilon in oplus
eps1 = 1 #2.  + 1J*2.01 #epsilon in ominus
eps0 = 3.5 + 1J*3.51 #epsilon in the layer

#FE parameters
order = 2 #order of the integration in the error computation
gamma=-1.e2*1J
hmax=0.4
res=0.4
print('Meshing parameter',2*3.141/k/hmax*order)

#Parameter for the FF computation
collection='interface' # There is no need to put collection on the boundary of D but outside must be air
print('Computing FF using surface ',collection)
nphi=4 # has to be greater than 2
nang=nphi #has to be even

def MakeGeometry():
    geometry = CSGeometry()
    o_ext = (Sphere(Pnt(0,0,0), Rext)).bc("outer")

    pml = Sphere(Pnt(0,0,0),Rpml)

    o_plus = Sphere(Pnt(0,0,0), Rplus).bc("interface")
    o_minus = (Sphere(Pnt(0,0,0), Rminus)).maxh(0.2).mat("ominus")
    other = Sphere(Pnt(0,c,0),Rother)
    #other = Cylinder(Pnt(0,0,0),Pnt(0,1,0),0.2) * OrthoBrick(Pnt(-0.4,0,-0.4),Pnt(0.4,0.6,0.4))
    
    geometry.Add ((o_ext - pml).mat("air"))
    geometry.Add ((pml-o_plus).mat("air"))
    geometry.Add (((o_plus-o_minus)-other).mat("oplus"))
    geometry.Add ((other-o_minus).mat("olayer"))
    geometry.Add (o_minus)

    return geometry



ngmesh = MakeGeometry().GenerateMesh(maxh=hmax)
#ngmesh.Save("scatterer.vol")
mesh = Mesh(ngmesh)

ngsglobals.testout = "test.out"
# First, set PML parameters for the ATC model
SetPMLParameters(rad=Rpml,alpha=.6)
# curve elements for geometry approximation
mesh.Curve(4)

viewoptions.clipping.enable = 1
viewoptions.clipping.ny = 0
viewoptions.clipping.nz = -1

Draw(mesh)

ngsglobals.msg_level = 5

Vext = HCurl(mesh, order=2, complex=True, definedon=mesh.Materials("air")+mesh.Materials("pml"), dirichlet="outer")

Vplus_l = HCurl(mesh, order=2, complex=True, definedon=mesh.Materials("oplus")+mesh.Materials("olayer")+mesh.Materials("ominus"))

V_l=FESpace([Vext,Vplus_l])


uext_l,uplus_l = V_l.TrialFunction()
vext_l,vplus_l = V_l.TestFunction()

uext_b,uplus_b = V_l.TrialFunction()


#Material properties for the layer
mur_l = {"ominus" : mu1 , "oplus" : mu2 , "olayer" : mu0 , "air" : 1 , "pml" : 1}
nu_l_coef = [ 1/(mur_l[mat]) for mat in mesh.GetMaterials() ]
nu_l = CoefficientFunction ([ 1/(mur_l[mat]) for mat in mesh.GetMaterials() ])

epsilonr_l = {"ominus" : eps1 , "oplus" : eps2 , "olayer" : eps0 , "air" : 1 , "pml" : 1}
epsilon_coef_l = [ epsilonr_l[mat] for mat in mesh.GetMaterials() ]
epsilon_l = CoefficientFunction (epsilon_coef_l)

#Material properties for the layer
mur_b = {"ominus" : mu1 , "oplus" : mu2 , "olayer" : mu2 , "air" : 1 , "pml" : 1}
nu_b_coef = [ 1/(mur_l[mat]) for mat in mesh.GetMaterials() ]
nu_b = CoefficientFunction ([ 1/(mur_l[mat]) for mat in mesh.GetMaterials() ])

epsilonr_b = {"ominus" : eps1 , "oplus" : eps2 , "olayer" : eps2 , "air" : 1 , "pml" : 1}
epsilon_coef_b = [ epsilonr_l[mat] for mat in mesh.GetMaterials() ]
epsilon_b = CoefficientFunction (epsilon_coef_l)

nv = specialcf.normal(mesh.dim) # normal vector

# Set PML parameters for the "exact" model
SetPMLParameters(rad=Rpml,alpha=0.6)

#Sesquilinear form for the layer
a_l = BilinearForm(V_l, symmetric=True, flags={"printelmat":True})

a_l.components[0] += BFI("PML_curlcurledge",coef=nu_l)
a_l.components[0] += BFI("PML_massedge",coef=-k*k*epsilon_l)

a_l += SymbolicBFI(nu_l*curl(uplus_l)*curl(vplus_l)   - k*k*epsilon_l*uplus_l*vplus_l)


a_l += SymbolicBFI(0.5*( (1./mu2)*curl(uplus_l).Trace() + curl(uext_l).Trace()) * cross(nv,vext_l.Trace()-vplus_l.Trace() ) ,BND,definedon=mesh.Boundaries("interface"))
a_l += SymbolicBFI( cross(nv,uext_l.Trace()-uplus_l.Trace()) * 0.5*( (1./mu2)*curl(vplus_l).Trace() + curl(vext_l).Trace() ) ,BND,definedon=mesh.Boundaries("interface"))
a_l += SymbolicBFI( gamma/hmax*1.J*k* (uext_l.Trace() - uplus_l.Trace()) * (vext_l.Trace() - vplus_l.Trace()),BND,definedon=mesh.Boundaries("interface"))

a_l.Assemble()

# Set PML parameters for the "exact" model
#SetPMLParameters(rad=Rpml,alpha=0.6)

#Sesquilinear form for the layer
a_b = BilinearForm(V_l, symmetric=True, flags={"printelmat":True})

a_b.components[0] += BFI("PML_curlcurledge",coef=nu_b)
a_b.components[0] += BFI("PML_massedge",coef=-k*k*epsilon_b)

a_b += SymbolicBFI(nu_b*curl(uplus_b)*curl(vplus_l)   - k*k*epsilon_b*uplus_b*vplus_l)


a_b += SymbolicBFI(0.5*( (1./mu2)*curl(uplus_b).Trace() + curl(uext_b).Trace()) * cross(nv,vext_l.Trace()-vplus_l.Trace() ) ,BND,definedon=mesh.Boundaries("interface"))
a_b += SymbolicBFI( cross(nv,uext_b.Trace()-uplus_b.Trace()) * 0.5*( (1./mu2)*curl(vplus_l).Trace() + curl(vext_l).Trace() ) ,BND,definedon=mesh.Boundaries("interface"))
a_b += SymbolicBFI( gamma/hmax*1.J*k* (uext_b.Trace() - uplus_b.Trace()) * (vext_l.Trace() - vplus_l.Trace()),BND,definedon=mesh.Boundaries("interface"))

a_b.Assemble()


fff=open('FFdefect.txt','w')

rhs=open('Epl_back.txt','w')

print(k,file=fff)

print(nang,file=fff)

print(nphi,file=fff)

print(nang,file=rhs)

print(nphi,file=rhs)

for ii in range(0,nphi):

    phii = pi*ii/(nphi-1)

    if ii == 0:

        d=(0,0,1)
        p=(-1,1,0)
        ll=0
        #for ll in range(0,2):
            #if ll==1:
                #p=pp #polarization vector
            #print('ll=',ll, 'and p=(', p[0],',',p[1],',',p[2],')')
        #Incident field
        Ein= exp(k*1J*(d[0]*x+d[1]*y+d[2]*z))*CoefficientFunction(p)
        Einext=GridFunction(Vext,"gEin")
        Einext.Set(Ein)
        curlEin = k*1.J*cross(d,p)*exp(k*1.J*(d[0]*x+d[1]*y+d[2]*z))
        curlEinext = GridFunction(Vext,"gcurlEin")
        curlEinext.Set(curlEin)

        #RHS for the layer
        f_l = LinearForm(V_l)
        f_l += SymbolicLFI( curlEin * 0.5 * (cross(nv,vplus_l.Trace() + vext_l.Trace() )),BND,definedon=mesh.Boundaries("interface"))
        f_l += SymbolicLFI(Ein * 0.5*cross(nv,(1./mu2)*curl(vplus_l).Trace()+curl(vext_l).Trace() ),BND,definedon=mesh.Boundaries("interface"))
        f_l += SymbolicLFI(-gamma/hmax*1.J*k* Ein * (vext_l.Trace()-vplus_l.Trace()),BND,definedon=mesh.Boundaries("interface"))

        #Define solution functions
        u_l = GridFunction(V_l)#bubble Solution

        u_b = GridFunction(V_l)#bubble Solution

        f_l.Assemble()

        u_l.vec.data += a_l.mat.Inverse(V_l.FreeDofs()) * f_l.vec

        u_b.vec.data += a_b.mat.Inverse(V_l.FreeDofs()) * f_l.vec

        viewoptions.clipping.enable = 1
        viewoptions.clipping.ny = 0
        viewoptions.clipping.nz = -1
        visoptions.clipsolution = "scal"

        nopml = CoefficientFunction ([0 if mat=='pml' else 1 for mat in mesh.GetMaterials() ])

        #Draw(Einext,mesh,"Einext")
        #Draw(CoefficientFunction(u_l.components[1])+nopml*CoefficientFunction(u_l.components[0])+nopml*Einext,mesh,"E")
        #Draw(u_l.components[0],mesh,"uext")
        #Draw(Ein,mesh,"Ein")
        #Redraw()

        Es = u_l.components[0]
        Eb = u_b.components[1]
        findff6(ii,jj,k,ll,Es,order,mesh,nv,d,p,collection,nang,nphi,Vext,fff) # this version outputs to a file
        writebackground(ii,jj,Eb,nang,nphi,rhs) # this version outputs to a file
    elif ii == nphi-1:

        d=(0,0,-1)
        p=(-1,1,0)

        ll=0
        #for ll in range(0,2):
            #if ll==1:
                #p=pp #polarization vector
            #print('ll=',ll, 'and p=(', p[0],',',p[1],',',p[2],')')
        #Incident field
        Ein= exp(k*1J*(d[0]*x+d[1]*y+d[2]*z))*CoefficientFunction(p)
        Einext=GridFunction(Vext,"gEin")
        Einext.Set(Ein)
        curlEin = k*1.J*cross(d,p)*exp(k*1.J*(d[0]*x+d[1]*y+d[2]*z))
        curlEinext = GridFunction(Vext,"gcurlEin")
        curlEinext.Set(curlEin)

        #RHS for the layer
        f_l = LinearForm(V_l)
        f_l += SymbolicLFI( curlEin * 0.5 * (cross(nv,vplus_l.Trace() + vext_l.Trace() )),BND,definedon=mesh.Boundaries("interface"))
        f_l += SymbolicLFI(Ein * 0.5*cross(nv,(1./mu2)*curl(vplus_l).Trace()+curl(vext_l).Trace() ),BND,definedon=mesh.Boundaries("interface"))
        f_l += SymbolicLFI(-gamma/hmax*1.J*k* Ein * (vext_l.Trace()-vplus_l.Trace()),BND,definedon=mesh.Boundaries("interface"))

        #Define solution functions
        u_l = GridFunction(V_l)#bubble Solution

        u_b = GridFunction(V_l)#bubble Solution

        f_l.Assemble()

        u_l.vec.data += a_l.mat.Inverse(V_l.FreeDofs()) * f_l.vec

        u_b.vec.data += a_b.mat.Inverse(V_l.FreeDofs()) * f_l.vec

        viewoptions.clipping.enable = 1
        viewoptions.clipping.ny = 0
        viewoptions.clipping.nz = -1
        visoptions.clipsolution = "scal"

        nopml = CoefficientFunction ([0 if mat=='pml' else 1 for mat in mesh.GetMaterials() ])

        #Draw(Einext,mesh,"Einext")
        #Draw(CoefficientFunction(u_l.components[1])+nopml*CoefficientFunction(u_l.components[0])+nopml*Einext,mesh,"E")
        #Draw(u_l.components[0],mesh,"uext")
        #Draw(Ein,mesh,"Ein")
        #Redraw()

        Es = u_l.components[0]
        Eb = u_b.components[1]
        findff6(ii,jj,k,ll,Es,order,mesh,nv,d,p,collection,nang,nphi,Vext,fff) # this version outputs to a file
        writebackground(ii,jj,Eb,nang,nphi,rhs) # this version outputs to a file
    else:

        for jj in range(0,nang):

            thetajj = 2*pi*jj/(nang-1)

            d=(sin(phii)*cos(thetajj),sin(phii)*sin(thetajj),cos(phii))#incident direction
            p=(d[1]-d[2],d[2]-d[0],d[0]-d[1])
            #normp=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])
            #p=(p[0]/normp,p[1]/normp,p[2]/normp)
            #pp= cross(p,d)
            ll=0
            #for ll in range(0,2):
                #if ll==1:
                    #p=pp #polarization vector
                #print('ll=',ll, 'and p=(', p[0],',',p[1],',',p[2],')')
            #Incident field
            Ein= exp(k*1J*(d[0]*x+d[1]*y+d[2]*z))*CoefficientFunction(p)
            Einext=GridFunction(Vext,"gEin")
            Einext.Set(Ein)
            curlEin = k*1.J*cross(d,p)*exp(k*1.J*(d[0]*x+d[1]*y+d[2]*z))
            curlEinext = GridFunction(Vext,"gcurlEin")
            curlEinext.Set(curlEin)

            #RHS for the layer
            f_l = LinearForm(V_l)
            f_l += SymbolicLFI( curlEin * 0.5 * (cross(nv,vplus_l.Trace() + vext_l.Trace() )),BND,definedon=mesh.Boundaries("interface"))
            f_l += SymbolicLFI(Ein * 0.5*cross(nv,(1./mu2)*curl(vplus_l).Trace()+curl(vext_l).Trace() ),BND,definedon=mesh.Boundaries("interface"))
            f_l += SymbolicLFI(-gamma/hmax*1.J*k* Ein * (vext_l.Trace()-vplus_l.Trace()),BND,definedon=mesh.Boundaries("interface"))

            #Define solution functions
            u_l = GridFunction(V_l)#bubble Solution

            u_b = GridFunction(V_l)#bubble Solution

            f_l.Assemble()

            u_l.vec.data += a_l.mat.Inverse(V_l.FreeDofs()) * f_l.vec

            u_b.vec.data += a_b.mat.Inverse(V_l.FreeDofs()) * f_l.vec

            viewoptions.clipping.enable = 1
            viewoptions.clipping.ny = 0
            viewoptions.clipping.nz = -1
            visoptions.clipsolution = "scal"

            nopml = CoefficientFunction ([0 if mat=='pml' else 1 for mat in mesh.GetMaterials() ])

            #Draw(Einext,mesh,"Einext")
            #Draw(CoefficientFunction(u_l.components[1])+nopml*CoefficientFunction(u_l.components[0])+nopml*Einext,mesh,"E")
            #Draw(u_l.components[0],mesh,"uext")
            #Draw(Ein,mesh,"Ein")
            #Redraw()

            Es = u_l.components[0]
            Eb = u_b.components[1]
            findff6(ii,jj,k,ll,Es,order,mesh,nv,d,p,collection,nang,nphi,Vext,fff) # this version outputs to a file
            writebackground(ii,jj,Eb,nang,nphi,rhs) # this version outputs to a file

rhs.close()
fff.close()

