#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ = \
"""
pdb_water-contact.py

Determine which residues are in "contact" with other residues based on sum of
vdw radii and water radius.
"""
__author__ = "Michael J. Harms"
__date__ = "080214"

from .helper import geometry
from .data import common


def pdbWaterContact(pdb,water_radii=1.4):
    """
    Residues are considered neighbors if *any* atom of theirs approaches closer
    than vdw1 + vd2 + water.
    """

    # Grab atoms
    atoms = [l for l in pdb if l[0:4] == "ATOM"]

    # Grab coordinates and calculate distance matrix
    coord = [[float(l[30+8*i:38+8*i]) for i in range(3)] for l in atoms]
    dist = geometry.calcDistances(coord)

    # Grab list of unique residues
    residues = []
    for line in atoms:
        if line[21:26] not in residues:
            residues.append(line[21:26])
    neighbors = dict([(r,[]) for r in residues])

    wat_sep = 2*water_radii
    for i in range(len(atoms)):
        for j in range(i+1,len(atoms)):
            max_sep = wat_sep + \
                      common.VDW_DICT[atoms[i][13:16]] + \
                      common.VDW_DICT[atoms[j][13:16]]

            if dist[i][j] <= max_sep:
                if atoms[j][21:26] not in neighbors[atoms[i][21:26]]:
                    neighbors[atoms[i][21:26]].append(atoms[j][21:26])
                    neighbors[atoms[j][21:26]].append(atoms[i][21:26])

    return neighbors
