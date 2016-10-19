# EAM7Routine
Embedded Atom Method Hydrogen Nickel potential

https://doi.org/10.5281/zenodo.161855


[![GPL licensed](https://img.shields.io/aur/license/yaourt.svg)](http://github.com/jbr36/EAM7Routine/blob/master/license.md)

## Synopsis

This is an implementation of EAM7 Ni(100)/H for the interaction of a single H atom with a Ni(100). 
It is part of a comparison study with an improved Embedded Atom Method (EAM) potential. Literature:

M. S. Daw and M. I. Baskes, Phys. Rev. B 29, 6443 (1984)

Wonchoba, S. E., Truhlar, D. G. (1996) Phys. Rev. B 53:11222-11241 and REV B 51, 9985 (1995)

James T. Kindt and John C. Tully, J. Chem. Phys. 111, 11060 (1999)

## Code (Example and tests)

The input is an xyz structure of a Nickel(100) surface and a Hydrogen atom.
The output is the potential energy of the input structure and its analytical gradients.

Examples and tests will be published with the corresponding journal publication.

## Installation

The code provides a subroutine which can easily be linked into main MD, RPMD, OPTIM etc. drivers.
