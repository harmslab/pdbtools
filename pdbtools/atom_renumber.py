#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

# Thanks to Blatter Markus for adding the capability to deal with HETATM and
# NMR models

__description__ = \
"""
pdb_atom_renumber.py

Renumbers atom records in a pdb without touching residue numbers.
"""

__author__ = "Michael J. Harms"
__date__ = "070529"

import os, sys

def pdbAtomRenumber(pdb,renumber_het=True):
    """
    Renumber all atoms in pdb file, starting from 1.
    """

    entries_to_renumber = ["ATOM  ","TER   ","ANISOU"]
    if renumber_het == True:
        entries_to_renumber.append("HETATM")

    out = []
    counter = 1
    for line in pdb:
        # For and ATOM record, update residue number
        if line[0:6] in entries_to_renumber:
            out.append("%s%5s%s" % (line[0:6],counter,line[11:]))
            counter += 1
        else:

            # reset the counter for a new model
            if line[0:6] == "ENDMDL":
                counter = 1

            out.append(line)

    return out
