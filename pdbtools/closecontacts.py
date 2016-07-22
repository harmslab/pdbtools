#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ = \
"""
pdb_close-contacts.py

Create a list of atoms that are within some distance cutoff of each other in
a pdb file.
"""

__author__ = "Michael J. Harms"
__date__ = ""

import os, sys
from .helper import geometry
from .contactplot import pdbContact

def pdbCloseContacts(pdb,ca_only,distance_cutoff):
    """
    Return a list of atom-atom distances that are less than distance_cutoff.
    It ignores atom-atom distances within a given residue.  If ca_only is True,
    it will only look at CA atoms.
    """

    pdb = [l for l in pdb if l[:4] == "ATOM"]

    distances = pdbContact(pdb,not ca_only)

    if ca_only:
        pdb = [l for l in pdb if l[13:16] == "CA "]

    output = []
    for i in range(len(pdb)):
        for j in range(i+1,len(pdb)):

            # Don't look at distances within the same residue
            if pdb[i][17:26] != pdb[j][17:26]:
                if distances[i][j] < distance_cutoff:

                    # Also skip consecutive backbone bonds
                    if int(pdb[j][22:26]) == int(pdb[i][22:26]) + 1:
                        if (pdb[i][13] == "C" and pdb[j][13] == "N") or \
                           (pdb[i][13] == "O" and pdb[j][13] == "N") or \
                           (pdb[i][13] == "C" and pdb[j][13:15] == "CA"):
                            continue


                    output.append((pdb[i][13:26],pdb[j][13:26],distances[i][j]))

    return output
