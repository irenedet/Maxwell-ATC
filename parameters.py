from netgen.csg import *
from ngsolve import *
from ngsolve.internal import *
from netgen.meshing import SetTestoutFile


############################ Parameters

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
nphi=3 # has to be greater than 2 and odd
nang=2 #has to be even

mu2inv=1./(mu2)
lambda0=2.*pi/k
print('Wavelength in background medium is ',lambda0)
print('Wave number k is ',k)
print('Epsilon in olayer is',eps0)
print('Mesh size h for FF mesh is',hmax)
data_dir='LData_'+str(k)+'_'+str(eps0)
print('Saving to '+ data_dir)

