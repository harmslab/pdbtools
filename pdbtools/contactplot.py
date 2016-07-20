#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ = \
"""
pdb_contact.py

Generate a distance-based contact plot for a protein.
"""

__author__ = "Michael J. Harms"
__date__ = ""

import os, sys
from .helper import geometry

def pdbContact(pdb,all=False):
    """
    Generate array of all atom-atom distances in a pdb file.  Defaults to only
    look at CA atoms
    """

    # Remove non-CA atoms
    if not all:
        pdb = [l for l in pdb if l[0:4] == "ATOM" and l[13:16] == "CA "]
    else:
        pdb = [l for l in pdb if l[0:4] == "ATOM"]

    coord = [[float(l[30+8*i:38+8*i]) for i in range(3)] for l in pdb]

    dist = geometry.calcDistances(coord)

    return dist
