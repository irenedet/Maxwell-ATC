from netgen.csg import *
from ngsolve import *
from ngsolve.internal import *
from netgen.meshing import SetTestoutFile
from my_geometry import MakeGeometry,sphere_geom,cube_geom,L_geom
from functions import findff6,writebackground, cross, Ocross 
from parameters import *

##########################################################################################
#                      Setting the Geometry                                              #
##########################################################################################
mesh,nv,collection = MakeGeometry(Rext,Rpml,Rminus,Rplus,c,Rother,res,hmax) # Create the mesh
ngsglobals.msg_level = 5

NGP=12;
hair=lambda0/NGP #mesh size in the air
print('Mesh size in air is ',hair)
modeps0 = sqrt(eps0.real**2+eps0.imag**2)
#res=hair/modeps0 #mesh size in the internal layer
print('Mesh size in ominus',res)

##########################################################################################
#                 Definition of the Variational problem                                  #
##########################################################################################
Vext = HCurl(mesh, order=2, complex=True, definedon=mesh.Materials("air")+mesh.Materials("pml"), dirichlet="outer")
Vplus_l = HCurl(mesh, order=2, complex=True, definedon=mesh.Materials("oplus")+mesh.Materials("olayer")+mesh.Materials("ominus"))
V_l=FESpace([Vext,Vplus_l])


uext_l,uplus_l = V_l.TrialFunction()
vext_l,vplus_l = V_l.TestFunction()
uext_b,uplus_b = V_l.TrialFunction()


################ Functions of the material properties for the defective ("exact") problem
mur_l = {"ominus" : mu1 , "oplus" : mu2 , "olayer" : mu0 , "air" : 1 , "pml" : 1}
nu_l_coef = [ 1/(mur_l[mat]) for mat in mesh.GetMaterials() ]
nu_l = CoefficientFunction ([ 1/(mur_l[mat]) for mat in mesh.GetMaterials() ])

epsilonr_l = {"ominus" : eps1 , "oplus" : eps2 , "olayer" : eps0 , "air" : 1 , "pml" : 1}
epsilon_coef_l = [ epsilonr_l[mat] for mat in mesh.GetMaterials() ]
epsilon_l = CoefficientFunction (epsilon_coef_l)

################ Functions of the material properties for the background problem
mur_b = {"ominus" : mu1 , "oplus" : mu2 , "olayer" : mu2 , "air" : 1 , "pml" : 1}
nu_b_coef = [ 1/(mur_l[mat]) for mat in mesh.GetMaterials() ]
nu_b = CoefficientFunction ([ 1/(mur_l[mat]) for mat in mesh.GetMaterials() ])

epsilonr_b = {"ominus" : eps1 , "oplus" : eps2 , "olayer" : eps2 , "air" : 1 , "pml" : 1}
epsilon_coef_b = [ epsilonr_l[mat] for mat in mesh.GetMaterials() ]
epsilon_b = CoefficientFunction (epsilon_coef_l)

#nv = specialcf.normal(mesh.dim) # normal vector

SetPMLParameters(rad=Rpml,alpha=0.6) # Set PML parameters for the "exact" model

################ Sesquilinear form for the layer
a_l = BilinearForm(V_l, symmetric=True, flags={"printelmat":True})
a_l.components[0] += BFI("PML_curlcurledge",coef=nu_l)
a_l.components[0] += BFI("PML_massedge",coef=-k*k*epsilon_l)
a_l += SymbolicBFI(nu_l*curl(uplus_l)*curl(vplus_l)   - k*k*epsilon_l*uplus_l*vplus_l)
a_l += SymbolicBFI(0.5*( (1./mu2)*curl(uplus_l).Trace() + curl(uext_l).Trace()) * cross(nv,vext_l.Trace()-vplus_l.Trace() ) ,BND,definedon=mesh.Boundaries("interface"))
a_l += SymbolicBFI( cross(nv,uext_l.Trace()-uplus_l.Trace()) * 0.5*( (1./mu2)*curl(vplus_l).Trace() + curl(vext_l).Trace() ) ,BND,definedon=mesh.Boundaries("interface"))
a_l += SymbolicBFI( gamma/hmax*1.J*k* (uext_l.Trace() - uplus_l.Trace()) * (vext_l.Trace() - vplus_l.Trace()),BND,definedon=mesh.Boundaries("interface"))

# Set PML parameters for the "exact" model
# SetPMLParameters(rad=Rpml,alpha=0.6)

################ Sesquilinear form for the layer
a_b = BilinearForm(V_l, symmetric=True, flags={"printelmat":True})
a_b.components[0] += BFI("PML_curlcurledge",coef=nu_b)
a_b.components[0] += BFI("PML_massedge",coef=-k*k*epsilon_b)
a_b += SymbolicBFI(nu_b*curl(uplus_b)*curl(vplus_l)   - k*k*epsilon_b*uplus_b*vplus_l)
a_b += SymbolicBFI(0.5*( (1./mu2)*curl(uplus_b).Trace() + curl(uext_b).Trace()) * cross(nv,vext_l.Trace()-vplus_l.Trace() ) ,BND,definedon=mesh.Boundaries("interface"))
a_b += SymbolicBFI( cross(nv,uext_b.Trace()-uplus_b.Trace()) * 0.5*( (1./mu2)*curl(vplus_l).Trace() + curl(vext_l).Trace() ) ,BND,definedon=mesh.Boundaries("interface"))
a_b += SymbolicBFI( gamma/hmax*1.J*k* (uext_b.Trace() - uplus_b.Trace()) * (vext_l.Trace() - vplus_l.Trace()),BND,definedon=mesh.Boundaries("interface"))

################ Assembly of the FEM matrices
Es0b = GridFunction(V_l,name="Es0b") # Create discrete function for Es
Es0l = GridFunction(V_l,name="Es0l") # Create discrete function for Es

Es1b = GridFunction(V_l,name="Es1b") # Create discrete function for Es
Es1l = GridFunction(V_l,name="Es1l") # Create discrete function for Es



noPml = CoefficientFunction([0 if mat=='pml' else 1 for mat in mesh.GetMaterials()]) # Indicator function of the region without the pml


################ Setup and factor FEM matrices
with TaskManager():

    a_l.Assemble()
    a_b.Assemble()
    Minv_l=a_l.mat.Inverse(V_l.FreeDofs(),inverse='pardiso') # Invert global matrix
    Minv_b=a_b.mat.Inverse(V_l.FreeDofs(),inverse='pardiso') # Invert global matrix

################# Loop for the solution and computation of farfields, with many incident plane waves:
fff=open('FFdefect.txt','w') # file to write far field data for the defective problem
rhs=open('Epl_back.txt','w') # file to write the total field of the background problem for the rhs of the LSM

print(k,file=fff)
print(nang,file=fff)
print(nphi,file=fff)
print(nang,file=rhs)
print(nphi,file=rhs)

# Notation:
# (ii,jj) = dhat incident directions
# (kk,mm) = xhat directions of observation
 
for ii in range(0,nphi):
    print('starting ii=',ii)

    xii = 1-2*ii/(nphi-1)
    phii = -pi*(xii-1)*(xii**2+0.5*xii+0.5) #non homogeneous discretization of phi direction

    if ii == 0:
        jj = 0
        d=(0,0,1)
        p=(-1,1,0)
        pp=(-1,-1,0)
        for ll in range(0,2):
	# ll polarization label
            if ll==1:
                p=pp #polarization vector
            print('ll=',ll, 'and p=(', p[0],',',p[1],',',p[2],')')
            #Incident field
            Ein= exp(k*1J*(d[0]*x+d[1]*y+d[2]*z))*CoefficientFunction(p)
            Einext=GridFunction(Vext,"gEin")
            Einext.Set(Ein)
            curlEin = k*1.J*cross(d,p)*exp(k*1.J*(d[0]*x+d[1]*y+d[2]*z))
            curlEinext = GridFunction(Vext,"gcurlEin")
            curlEinext.Set(curlEin)

            #RHS for the layer
            f = LinearForm(V_l)
            f += SymbolicLFI( curlEin * 0.5 * (cross(nv,vplus_l.Trace() + vext_l.Trace() )),BND,definedon=mesh.Boundaries("interface"))
            f += SymbolicLFI(Ein*0.5*cross(nv,mu2inv*curl(vplus_l).Trace()+curl(vext_l).Trace() ),BND,definedon=mesh.Boundaries("interface"))
            f += SymbolicLFI(-gamma/hmax*1.J*k* Ein * (vext_l.Trace()-vplus_l.Trace()),BND,definedon=mesh.Boundaries("interface"))

            with TaskManager():
            	f.Assemble()
            	Es0l.vec.data += Minv_l * f.vec
            	Es0b.vec.data += Minv_b * f.vec

            #viewoptions.clipping.enable = 1
            #viewoptions.clipping.ny = 0
            #viewoptions.clipping.nz = -1
            #visoptions.clipsolution = "scal"

        
            #Draw(Einext,mesh,"Einext")
            #Draw(CoefficientFunction(u_l.components[1])+nopml*CoefficientFunction(u_l.components[0])+nopml*Einext,mesh,"E")
            #Draw(u_l.components[0],mesh,"uext")
            #Draw(Ein,mesh,"Ein")
            #Redraw()

            Es = Es0l.components[0]
            Eb = Es0b.components[1]
            findff6(ii,jj,k,ll,Es,order,mesh,nv,d,p,collection,nang,nphi,Vext,fff) # this version outputs to a file
            writebackground(ii,jj,Eb,nang,nphi,rhs) # this version outputs to a file
    
    elif ii == nphi-1:
        jj = 0
        d=(0,0,-1)
        p=(-1,1,0)
        pp=(-1,-1,0)
        for ll in range(0,2):
            if ll==1:
                p=pp #polarization vector
            print('ll=',ll, 'and p=(', p[0],',',p[1],',',p[2],')')
            #Incident field
            Ein= exp(k*1J*(d[0]*x+d[1]*y+d[2]*z))*CoefficientFunction(p)
       	    Einext=GridFunction(Vext,"gEin")
            Einext.Set(Ein)
            curlEin = k*1.J*cross(d,p)*exp(k*1.J*(d[0]*x+d[1]*y+d[2]*z))
            curlEinext = GridFunction(Vext,"gcurlEin")
            curlEinext.Set(curlEin)

            #RHS for the layer
       	    f = LinearForm(V_l)
            f += SymbolicLFI( curlEin * 0.5 * (cross(nv,vplus_l.Trace() + vext_l.Trace() )),BND,definedon=mesh.Boundaries("interface"))
            f += SymbolicLFI(Ein*0.5*cross(nv,mu2inv*curl(vplus_l).Trace()+curl(vext_l).Trace() ),BND,definedon=mesh.Boundaries("interface"))
            f += SymbolicLFI(-gamma/hmax*1.J*k* Ein * (vext_l.Trace()-vplus_l.Trace()),BND,definedon=mesh.Boundaries("interface"))

            with TaskManager():
            	f.Assemble()
            	Es0l.vec.data += Minv_l * f.vec
            	Es0b.vec.data += Minv_b * f.vec

            #Draw(Einext,mesh,"Einext")
            #Draw(CoefficientFunction(u_l.components[1])+nopml*CoefficientFunction(u_l.components[0])+nopml*Einext,mesh,"E")
            #Draw(u_l.components[0],mesh,"uext")
            #Draw(Ein,mesh,"Ein")
            #Redraw()

            Es = Es0l.components[0]
            Eb = Es0b.components[1]
            findff6(ii,jj,k,ll,Es,order,mesh,nv,d,p,collection,nang,nphi,Vext,fff) # this version outputs to a file
            writebackground(ii,jj,Eb,nang,nphi,rhs) # this version outputs to a file
    else:

        for jj in range(0,nang):

            thetajj = 2*pi*jj/(nang-1)

            d=(sin(phii)*cos(thetajj),sin(phii)*sin(thetajj),cos(phii))#incident direction
            p=(d[1]-d[2],d[2]-d[0],d[0]-d[1])
            #normp=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])
            #p=(p[0]/normp,p[1]/normp,p[2]/normp)
            pp= cross(p,d)
            for ll in range(0,2):
                if ll==1:
                    p=pp #polarization vector
                print('ll=',ll, 'and p=(', p[0],',',p[1],',',p[2],')')
                #Incident field
                Ein= exp(k*1J*(d[0]*x+d[1]*y+d[2]*z))*CoefficientFunction(p)
                Einext=GridFunction(Vext,"gEin")
                Einext.Set(Ein)
                curlEin = k*1.J*cross(d,p)*exp(k*1.J*(d[0]*x+d[1]*y+d[2]*z))
                curlEinext = GridFunction(Vext,"gcurlEin")
                curlEinext.Set(curlEin)
	
            	#RHS for the layer
                f = LinearForm(V_l)
                f += SymbolicLFI( curlEin * 0.5 * (cross(nv,vplus_l.Trace() + vext_l.Trace() )),BND,definedon=mesh.Boundaries("interface"))
                f += SymbolicLFI(Ein*0.5*cross(nv,mu2inv*curl(vplus_l).Trace()+curl(vext_l).Trace() ),BND,definedon=mesh.Boundaries("interface"))
                f += SymbolicLFI(-gamma/hmax*1.J*k* Ein * (vext_l.Trace()-vplus_l.Trace()),BND,definedon=mesh.Boundaries("interface"))

                with TaskManager():
                    f.Assemble()
                    Es0l.vec.data += Minv_l * f.vec
                    Es0b.vec.data += Minv_b * f.vec

          	#Draw(Einext,mesh,"Einext")
            	#Draw(CoefficientFunction(u_l.components[1])+nopml*CoefficientFunction(u_l.components[0])+nopml*Einext,mesh,"E")
            	#Draw(u_l.components[0],mesh,"uext")
            	#Draw(Ein,mesh,"Ein")
            	#Redraw()

                Es = Es0l.components[0]
                Eb = Es0b.components[1]
                findff6(ii,jj,k,ll,Es,order,mesh,nv,d,p,collection,nang,nphi,Vext,fff) # this version outputs to a file
                writebackground(ii,jj,Eb,nang,nphi,rhs) # this version outputs to a file

rhs.close()
fff.close()

