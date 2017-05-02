# Maxwell-ATC

This code is part of my PhD project.

The aim is to compare two models to simulate the scattering of electromagnetic waves in the presence of a thin layer of compact thickness that is located on just one part of the interface of two materials. 

The first model also called "full model" consists of meshing the whole domain and compute the transmission problem directly.
The second or "reduced model", has been derived using the concept of ATCs (Approximate Transmission Conditions). In our work, we derive a second order ATC model, where thee original efect by a thin layer is no longer computed, and instead of it we define suitable transmission conditions on the part of the interface where the original thin layer was located. These new transmission conditions involve second order surface operators, which require a more technical analytic approach, but has the advantage of providing an approximation that no longer requires a high mesh refinement.

This code is written using the Netgen/Ngsolve FEM library precisely to compare the L2 and Hcurl relative errors associated to these two models. The code is divided in the following files:

1) Brick_ATCerror.py is the main file where the solutions and the errors are computed.

2) Auxiliary files: aux_functions.py - where the sesquilinear form and the source terms are defined
                    parameters.py - where the physical, geometrical and numerical parameters are defined
                    mygeometry1.py - where the geometry of the numerical experiment is defined
