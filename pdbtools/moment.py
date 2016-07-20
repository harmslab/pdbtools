#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ = \
"""
pdb_moment.py

Calculate the approximate dipole moment of a protein.  This program needs polar
hydrogens to do its calculation.  The pdb file must be in the slightly non-
standard (UHBD) format written out by pdb_addH.  Calculation is done assuming
that Asp, Glu, Arg, Lys, and His are fully charged.  Charges are placed on the
CD, CG, CZ, NZ, and NE2 atoms respectively (see pdb_data/polar_param.txt).
"""
__author__ = "Michael J. Harms"
__date__ = ""

import os, sys
from .helper import geometry
from .data.polar import POLAR_CHARGE_DICT
from math import pi

class PdbMomentError(Exception):
    """
    General error class for this module.
    """

    pass

def calcMoment(charges,coord):
    """
    Actually do moment calculation on centered coordinates.
    """

    # unit conversion
    elementary_q = 1.602176487e-19 # elementary charge in C
    debye_conv = 3.33564e-30       # 1 Debye in C*m
    k = elementary_q*1e-10/debye_conv

    num_atoms = len(charges)

    # Calculate the moment by sum(q(i)*r(i)) for x, y, and z
    moment = [sum([coord[j][i]*charges[j] for j in range(num_atoms)])
            for i in range(3)]

    moment = [m*k for m in moment]

    return moment


def pdbMoment(pdb):
    """
    Calculate the dipole moment of a pdb file using the charge set in
    pdb_data.polar.
    """

    # Parse pdb file
    pdb = [l for l in pdb if l[0:4] == "ATOM"]

    # Grab fractional charge of each atom and coordinates from pdb file
    try:
        charges = [POLAR_CHARGE_DICT[l[17:21]][l[12:16]] for l in pdb]
    except KeyError:
        err = l
        err += "Residue/atom pair not recognized"
        raise PdbMomentError(err)

    coord = [[float(line[30+i*8:38+i*8])for i in range(3)] for line in pdb]

    num_atoms = len(charges)

    # Center coordinates on center of mass
    center = [sum([c[i] for c in coord])/num_atoms for i in range(3)]
    coord = [[c[i] - center[i] for i in range(3)] for c in coord]

    return calcMoment(charges,coord)
