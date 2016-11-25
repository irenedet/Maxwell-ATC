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
thanks to the mixed reciprocity principle proven for this setting of study, which relates the farfield pattern of
the Green's function of the background to the total field due to the background when the incident fild are plane waves.

3) defective.py : is a code in python that uses the FEM library Netgen/Ngsolve, to compute the farfield data of the problem with and without the defect. Its output are FFdefect.txt, FFbackground.txt and Epl_back.txt. These files that can be read in Matlab to generate A_d and A_b and b, using A_d = assemble_farfield('FFdefect.txt'), A_b = assemble_farfield('FFbackground.txt') and
b = assemble_rhs('Epl_back.txt').

****NOTE: In the current version the background is homogeneous so we only output FFdefect.txt and Epl_back.txt

