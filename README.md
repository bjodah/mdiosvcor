## Introduction
This open source python package MDIOSVCOR (Molcular Dynamics IOn SolVation CORrections) provides routines for calculating the correction terms for solvation parameters and partial molar variables of ion solvation calculated with molecular dynamics using the lattice summmation treatment of electrostatic interactions.

For installation instructions see below

Program for calculating the correction terms applied in the publication:

   Björn Dahlgren, Maria M. Reif, Philippe H. Hünenberger & Niels Hansen
   "Calculation of derivative thermodynamic hydration and aqueous partial
   molar properties of ions based on atomistic simulations"
   J. Chem. Theory Comput.

The program was implemented for use in a research project in the IGC group at ETH Zürich

If you make use the code in your own research, please cite the article referenced above.

## Installation
This package is written in python, i.e. no compilation is required, just execute "get_correction_terms.py" with appropriate flags. E.g:

     ./get_correction_terms.py -T 298.15 -P 101.3e3 -N 1024 -Y G -c all -i sod

(for further information on invocation see help by executing `./get_correction_terms.py --help`)

## Prerequisities
This package relies on:

- [Python](http://python.org) 2.7.x
- [Sympy](http://sympy.org) to treat algebraic expressions symbolically

## Additional information
Beyond the core equations used in calculating the correction terms specified in the article above,
this project reimplements the [IAPWS95](http://iapws.org) formaulation of thermodynamic properties
of ordinary water substance. It also uses a sophisticated method the retrieve the denisty of water,
and its temperature and/or pressure derivative of arbitrary order with arbitrary precision.
This is however due to the large number of non-linear terms compuationally expensive for higher
derivatives since the manipulations are done symbolically up to the last step of numerical solution
of the resulting implicit relation.

## Possible future extensions
Feel free to write them and make a pull request at github:
- Implement expressions not only for LS but also for CT
- Ideas?

## Participate
Any comments and/or imporvements of the code are greatly appreciated.
Improvements are best done by making a pull request at github
It is also the place to raise issues and filing bug reports.

Latest version is available at github.com/bjodah/mdiosvcor


## Contact information
Author: Björn Ingvar Dahlgren
Email (@gmail.com): bjodah

## License information
This work is open source and is released under the 2-clause BSD license (see LICENSE.txt for further information)

Copyright (c) 2011, 2012 Björn Ingvar Dahlgren
