Quick User Manual
=================

Date: July 11, 2010
By: Teepanis Chachiyo

1) Preparation
==============

To run Siam Quantum (or SQ for short) we would need a few files.

a) the exectuable, "sq.exe"
b) lapack and blass DLL library in case you are running in WindowsXP
c) basis set file such as 631g.txt or sto3g.txt. These are in GAMESS-US
   format and can be downloaded from the internet
d) the molecular structure in XYZ format. See several files like benzene.xyz
   in the "examples" folder

Suppose we save all the files above in the same directory, then all you
need to do is to execute.

# ./sq.exe

This will give you all available options. SQ is a command-line driven program.


2) Basic Calculation
====================

For example, you would like to compare the hartree-fock energy of a molecular 
oxygen between two states: singlet and triplet. Using the file o2.xyz in
"examples" directory, simply execute

# ./sq.exe o2.xyz 321g.txt -M=1

and 

# ./sq.exe o2.xyz 321g.txt -M=3

Here, option -M is for specifying the spin multiplicity of the molecule. There
are other options that you can play with. For example, use "-Q" to set the total
charge of the molecule.


3) Convergence Problems
=======================

In some cases, it might be difficult to get the solution to converge, specially
when the energy levels are highly degenerate like in metal complexes. Siam Quantum
uses a very simple convergence scheme: a weighted updating scheme. You can control
this by using option "-SCFDRAG=[REAL NUMBER BETWEEN 0..1]"

For example, 

-SCFDRAG=0.1   will usually always converge but it will take a lot of iteration
-SCFDRAG=0.9   will either converge very rapidly, or the solution goes into an
               oscillation

Try to set it as high as possible, but if it is "too" high the solution will go
into an oscillation.


4) View Electron Density or Molecular Oribitals
===============================================

Simple options like "-DENSITY" will create a file called "density.xsf" which can
be visualized using VMD program. It contains information about the electron density.

Use option "-MO=INTEGER" to create molecular orbital instead.


5) Other options you might want to experiment on
================================================

-SCFMAX=INTEGER   : set the number of maximum iteration
-SDMATRIX         : save density matrix file to restart the calculation
-LDMATRIX         : load density matrix file before starting the first iteration
-FDMATRIX=STRING  : set the density matrix file name (default "dmatrix.txt")


6) Future
=========

This project started out as a simple program for teaching students about 
electronic structure theory. However, the code is getting larger with more and
more features and is getting faster. 

Theerapon Khamlar, a graduate student in our research group, is working on
how to compute the gradient of energy. Once that is done, it will be possible
to compute the equilibrium structure of molecules.

We also need some help on how to improve the speed of 2-electron integrals;
and how to perform numerical integration of electron density, in which case,
it will be possible to implement DFT calculation too.

You might contact us at: teepanis@kku.ac.th   
or visit us on the web: http://www.physics.kku.ac.th/sq


Best wishes,

Teepanis Chachiyo.

