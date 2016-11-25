# Maxwell-ATC

This code is part of my PhD project.
The aim is to, based on electromagnetic far-field data, test a numerical algorithm 
for the identification of thin flaws (delamination) in a layered media.

The numerical algorithm was developed adapting qualitative inverse scattering techniques (the linear sampling method) to our 
problem. This algorithm is based on the solution of a severely ill-posed linear system Ax=b, where:

1) A = A_b - A_d, where A_b is the farfield matrix of the problem without the defect (the background medium), and
A_d is the farfield matrix of the material to be tested (possibly defective).

2) b is the far field pattern of the Green's function associated to the background media.

This project consists of the following files:
1) assemble_farfield.m : is a Matlab code that reads the electromagnetic farfield information computed by the finite element 
code written in the Netgen/Ngsolve library. With this function we assemble both A_b and A_d.

2) assemble_rhs.m : is a Matlab code that assembles the right-hand-side b, by reading
the total electromagnetic field data of the background media from a 
finite element code developed also in Netgen/Ngsolve. The total field of the background media is a fast way to compute b, 
thanks to the mixed reciprocity principle proven for this setting of study.
