#!/usr/bin/env python

# Copyright 2008, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ = \
"""
pdb_ligand.py

Spits out a list of all ligands in a pdb.
"""
__author__ = "Michael J. Harms"
__date__ = ""

import os

BORING_LIGANDS = ["HOH","CA","SO4","IOD","NA","CL","GOL","PO4"]

def pdbLigand(pdb,skip_boring=False):
    """
    """

    ligands = [l[7:10].strip() for l in pdb if l.startswith("HET   ")]
    ligands = list(dict([(l,[]) for l in ligands if l not in BORING_LIGANDS]).keys())

    return ";".join(ligands)
