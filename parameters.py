from ngsolve import *
from netgen.csg import *
from ngsolve.internal import *

#Geometrical parameters
k = 3. # Wave number
Rext = 2.7 #Radius of external circle
Rpml = 2. #Radius of PML
alpha_pml = 0.6 # absorbing coefficient of PML
Rplus =  1.3 #Radius of Oplus
Rminus = 0.6 #Radius of Ominus
delta = 0.0025 #Layer thickness

#Physical parameters
mu2 = 1. #mu in oplus
mu1 = 1.5 #1.5 #mu in ominus
mu0 = 1.5 #2 #mu in the layer

eps2 = 2  #+ 1J*1.01 #epsilon in oplus
eps1 = 3 #2.  + 1J*2.01 #epsilon in ominus
eps0 = 3 # + 1J*3.51 #epsilon in the layer

#FE parameters
order = 2 #order of the integration in the error computation
gamma=-1.e3
hmax=0.1#0.2#0.2
res=0.1
geom = 4 # 1-->Brick, 2-->Sphere, 4-->HalphSphere

d = (1,0,0)
p = (0,0,1)

#Geometry of the circular delamination
c = (2*delta+delta*delta)/(2+2*delta-sqrt(3)) #center of other
Rother = 1-c+delta # Radius of other

print('Meshing parameter',2*3.141/k/hmax*order)
# Far field pattern computation parameters
pi=4.*atan(1.)
FFh=0.4
data_dir='Data_'+str(k)+'_'+str(eps0)

# Visualization parameters 
plot_id = 1 # 1 = to plot the solutions

