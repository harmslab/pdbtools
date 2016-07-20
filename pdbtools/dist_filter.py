#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__author__ = "Michael J. Harms"
__date__ = "070709"
__description__ = \
"""
pdb_dist-filter.py

Takes pdb file and calculates the distance between "residue" "atom" and all
other atoms of this type in the pdb file with "column" matching "select".
"""

from math import sqrt

def extractCoor(line):
    """
    Take a line out of a pdb and extract coordinates.
    """

    return [float(line[30+8*i:38+8*i]) for i in range(3)]

def distFilter(pdb,residue,atom="N",column=[60,66],select="    NA"):
    """
    Calculate the distance between "residue" "atom" and the "atom" of residues
    that have column defined by "column" == "select".  Default is to compare
    nitrogen of "residue" to nitrogens or residues with b-factor == NA.
    """

    # Make sure atom entry is in correct format
    atom = "%-3s" % (atom.strip())

    # Grab only "atom" lines out of pdb
    pdb = [l for l in pdb if l[0:4] == "ATOM" and l[13:16] == atom]

    if len(pdb) == 0:
        err = "pdb file does not contain any atoms of type \"%s\"" % atom
        raise IOError(err)

    # Pull residue coordinates
    res_coord = [extractCoor(l) for l in pdb if l[22:26] == "%4i" % residue][0]

    # Pull selected residue coordinates
    try:
        na_coord = [extractCoor(l) for l in pdb
                    if l[column[0]:column[1]] == select]
    except IndexError:
        err = "Invalid column defined by line[%i:%i]" % (column[0],column[1])
        raise IOErro(err)

    if len(na_coord) == 0:
        err = "Column line[%i:%i] does not contain any \"%s\" entries!" % \
              (column[0],column[1],select)
        raise IOError(err)


    # Calculate distances
    dist = []
    for c in na_coord:
        dist.append(sqrt(sum([(c[i]-res_coord[i])**2 for i in range(3)])))

    return dist
