## Introduction
This open source python package MDIOSVCOR (Molecular Dynamics IOn SolVation CORrections) provides routines for calculating the correction terms for solvation parameters and partial molar variables of ion solvation calculated with molecular dynamics using the lattice summation treatment of electrostatic interactions.

For installation instructions see below

The program is used for calculating the correction terms applied in the publication:

   *Björn Dahlgren, Maria M. Reif, Philippe H. Hünenberger & Niels Hansen
   "Calculation of derivative thermodynamic hydration and aqueous partial
   molar properties of ions based on atomistic simulations"
   J. Chem. Theory Comput.*

The program was implemented for use in a research project in the IGC group at ETH Zürich.
The theory behind the correction terms were developed by several persons in the IGC group (see references in the article above)

If you make use the code in your own research, please cite the article referenced above.

## Installation
This package is written in python, i.e. no compilation is required.
However, in order to speed up calculations the default behavior of the script is to compile functions via FORTRAN to binaries. This is done on the fly the first time the program needs the binaries.

Therefore you just need to execute "calc_correction_terms.py" (or "calc_water_eps.py"/"calc_water_rho.py") with appropriate flags. E.g:

     ./calc_correction_terms.py -T 298.15 -P 101.3e3 -N 1024 -Y G -c all -i sod

(for further information on invocation see help by executing e.g. `./calc_correction_terms.py --help`)

## Prerequisites
This package relies on:

- [Python](http://python.org) 2.7.x
- [Sympy](http://sympy.org) 0.7.1 to treat algebraic expressions symbolically (Ubuntu package name: python-sympy)

Furthermore, for full functionality, additional requirements apply:

- f2py (provided by [NumPy](http://numpy.scipy.org/))
- a FORTRAN compiler (if compiled versions are wanted)

## More advanced usage
To calculate a large number of correction terms takes a rather long
while due to the symbolic treatment of the expressions and use of
explicit units on the numerical parameters. Therefore it might be of
interest to run a batch job and store the output in an easily parsed
format.

An example on how to run this kind of batch calculation for large number of different parameters is given below:

     python -c "import calc_correction_terms as cct; cct.batch_calc(['G','H'], [101.3e3], [273.15, 298.15, 323.15], [1024, 512], ['sod','cls'],['B','C1','C2','D'])"

## Additional information
Beyond the core equations used in calculating the correction terms specified in the article above,
this project re-implements the [IAPWS95](http://iapws.org) formulation of thermodynamic properties
of ordinary water substance. It also uses a sophisticated method the retrieve the density of water,
and its temperature and/or pressure derivative of arbitrary order with known numerical precision.
This is however due to the large number of non-linear terms computationally expensive for higher
order derivatives since the manipulations are done symbolically up to the last step of numerical
solution of the resulting implicit relation.

## Caveats
Cached compiled binaries are stored in cache/ directory. When moving between hosts of different architecture/OS
one must purge this directory of "*.so" files.

## Possible future extensions
Feel free to write them and make a pull request at github:
- Implement expressions not only for LS but also for CT
- Optimize performance

## Participate
Any comments and/or improvements of the code are greatly appreciated.
Improvements are best done by making a pull request at github
It is also the place to raise issues and filing bug reports.

Latest version is available at github.com/bjodah/mdiosvcor


## Contact information
- Author of python package: Björn Ingvar Dahlgren
- Email (@gmail.com): bjodah

## License information
This work is open source and is released under the 2-clause BSD license (see LICENSE.txt for further information)

Copyright (c) 2011, 2012 Björn Ingvar Dahlgren
