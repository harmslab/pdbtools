#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ = \
"""
pdb_coulomb.py

Calculates the total coulomb energy (kcal/mol) of a protein structure in a pdb
file assuming that all groups titrate with model compound pKa values.
"""

__author__ = "Michael J. Harms"
__date__ = "070520"

import sys, os, copy
from math import sqrt, exp
from .data.common import *

def readPDB(pdb_file):
    """
    Takes a pdb file and reads in the coordinates of each titratable group.
    Assigns pka and charge state of each group.
    """


    # Open pdb_file and read each line into pdb (a list of lines)
    f = open(pdb_file,'r')
    pdb = f.readlines()
    f.close()

    # Grab only ATOM entries that are titratable
    pdb = [l for l in pdb if l[0:4] == "ATOM" and
                             l[17:20] in list(TITR_ATOM.keys()) and
                             l[13:16] == TITR_ATOM[l[17:20]]]

    # Initalize lists to hold coordinates, pkas, and charge
    coord, pKa, charge = [], [], []

    # Go through each line in the pdb file
    for line in pdb:
        amino_acid = line[17:20]

        # Grab the xyz coordinates
        coord.append([float(line[30+8*i:38+8*i]) for i in range(3)])

        pKa.append(PKA_DICT[amino_acid])

        # Look up charge
        charge.append(CHARGE_DICT[amino_acid])

    # Return the coordinates, pka, and charge
    return coord, pKa, charge

def hendersonHasselbach(pKa,charge,pH):
    """
    Calculate the fractional charge on a group with pKa and charge at some
    pH value.
    """

    return charge/(1 + 10**(charge*(pH-pKa)))


def pdbCoulomb(coord,pKa,charge,dielec_const,ionic_str,pH,temperature):
    """
    Calculates the energy of a structure given the coordinates of each
    charged atom, their fractional charge, the dielectric constant, and the
    ionic strength.
    """

    ionic_str = ionic_str/1000

    # Initialize variables
    kappa = 50.29*sqrt(ionic_str/(dielec_const*temperature))
    num_groups = len(coord)
    energy = 0.

    hh_chg = [hendersonHasselbach(pKa[i],charge[i],pH)
              for i in range(num_groups)]

    # Calculate energy of interaction of every ij interaction (making sure not
    # to double count; note we start j at i + 1).
    for i in range(num_groups):
        for j in range(i+1,num_groups):

            # Calculate distance between atom i and atom j
            r = sqrt(sum([(coord[i][k]-coord[j][k])**2 for k in range(3)]))

            # Add the energy of this interaction to the total energy
            energy += 332*hh_chg[i]*hh_chg[j]/(r*dielec_const)*exp(-kappa*r)

    # Return energy
    return energy
