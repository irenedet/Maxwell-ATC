from ngsolve import *
from netgen.csg import *
from ngsolve.internal import *
import numpy as np

pi=4.*atan(1.)

############################ cross products
Ocross=lambda a,b: (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-b[0]*a[1])
cross=lambda a,b: CoefficientFunction((a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-b[0]*a[1]))
############################

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

def writeRHS(ka,FFP,Es0b,Es1b,order,mesh,nv,collection,fespace):
    #(ii,jj,ka,ll,Es,order,mesh,nv,dv,pv,collection,nang,nphi,fespace,fff):
    # This version computes the bistatic ff pattern and returns
    # it to the calling program, does both polarizations at once
    # print('Using findff6B....')
    w=fespace.TestFunction()
    FF0=[]
    FF1=[]
    np=0
    vv=GridFunction(fespace)
    zer0=CoefficientFunction((0.,0.,0.))
    vv.Set(zer0,VOL) 
    for p in FFP:
        #print('starting ',p.p)
        xhat=p.p
        nxhat=sqrt(xhat[0]**2+xhat[1]**2+xhat[2]**2)
        if nxhat>1-.0001:  # the ff mesh contains some interior points
            ff0=[0,0,0]
            ff1=[0,0,0]
            CEout = exp(1J*ka*(x*xhat[0]+y*xhat[1]+z*xhat[2]))
            # Conjugate to allow for conjugation in linear forms
            for j in range(0,3):
                e=[0.,0.,0.]
                e[j]=1.
                ## Functional for the full FF pattern
                func=cross(Ocross(e,xhat),nv)*CEout*(-1J*ka)
                gf = LinearForm(fespace)
                #gf.vec[:]=0
                gf += SymbolicLFI(func*w.Trace(),BND,definedon=
                                  mesh.Boundaries(collection))
                ## Functional for the second part using finite element cutoff
                eT=cross(xhat,Ocross(e,xhat))
                vv.Set(eT*CEout,BND,definedon=mesh.Boundaries(collection))
                gf += SymbolicLFI(curl(vv)*curl(w)-ka*ka*vv*w,VOL,
                                  definedon=mesh.Materials('air'))
                with TaskManager():
                    gf.Assemble()
                # Conjugates are needed since inner product conjugates
                # the second component(?)
                ff1a0=gf.vec.InnerProduct(Es0.vec)
                ff0[j]=ff1a0.conjugate()/4./pi
                ff1a1=gf.vec.InnerProduct(Es1.vec)
                ff1[j]=ff1a1.conjugate()/4./pi
            FF0.append(ff0)
            FF1.append(ff1)
        else:
            FF0.append([-9999,-9999,-9999])
            FF1.append([-9999,-9999,-9999])
        #print('FF Done ',np,' of ',len(FFP)-1)
        np=np+1
    return(FF0,FF1)
    
def findff6B2(ka,FFP,Es0,Es1,order,mesh,nv,collection,fespace):
    #(ii,jj,ka,ll,Es,order,mesh,nv,dv,pv,collection,nang,nphi,fespace,fff):
    # This version computes the bistatic ff pattern and returns
    # it to the calling program, does both polarizations at once
    # print('Using findff6B....')
    w=fespace.TestFunction()
    FF0=[]
    FF1=[]
    np=0
    vv=GridFunction(fespace)
    zer0=CoefficientFunction((0.,0.,0.))
    vv.Set(zer0,VOL) 
    for p in FFP:
        #print('starting ',p.p)
        xhat=p.p
        nxhat=sqrt(xhat[0]**2+xhat[1]**2+xhat[2]**2)
        if nxhat>1-.0001:  # the ff mesh contains some interior points
            ff0=[0,0,0]
            ff1=[0,0,0]
            CEout = exp(1J*ka*(x*xhat[0]+y*xhat[1]+z*xhat[2]))
            # Conjugate to allow for conjugation in linear forms
            for j in range(0,3):
                e=[0.,0.,0.]
                e[j]=1.
                ## Functional for the full FF pattern
                func=cross(Ocross(e,xhat),nv)*CEout*(-1J*ka)
                gf = LinearForm(fespace)
                #gf.vec[:]=0
                gf += SymbolicLFI(func*w.Trace(),BND,definedon=
                                  mesh.Boundaries(collection))
                ## Functional for the second part using finite element cutoff
                eT=cross(xhat,Ocross(e,xhat))
                vv.Set(eT*CEout,BND,definedon=mesh.Boundaries(collection))
                gf += SymbolicLFI(curl(vv)*curl(w)-ka*ka*vv*w,VOL,
                                  definedon=mesh.Materials('air'))
                with TaskManager():
                    gf.Assemble()
                # Conjugates are needed since inner product conjugates
                # the second component(?)
                ff1a0=gf.vec.InnerProduct(Es0.vec)
                ff0[j]=ff1a0.conjugate()/4./pi
                ff1a1=gf.vec.InnerProduct(Es1.vec)
                ff1[j]=ff1a1.conjugate()/4./pi
            FF0.append(ff0)
            FF1.append(ff1)
        else:
            FF0.append([-9999,-9999,-9999])
            FF1.append([-9999,-9999,-9999])
        #print('FF Done ',np,' of ',len(FFP)-1)
        np=np+1
    return(FF0,FF1)
