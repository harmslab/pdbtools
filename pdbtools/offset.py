#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ =\
"""
pdb_offset.py

Adds an offset to the residue numbers of a pdb without touching anything else.
"""

__author__ = "Michael J. Harms"
__date__ = "070131"


import os, sys

def pdbOffset(pdb_file,offset):
    """
    Adds an offset to the residue column of a pdb file without touching anything
    else.
    """

    # Read in the pdb file
    f = open(pdb_file,'r')
    pdb = f.readlines()
    f.close()

    out = []
    for line in pdb:
        # For and ATOM record, update residue number
        if line[0:6] == "ATOM  " or line[0:6] == "TER   ":
            num = offset + int(line[22:26])
            out.append("%s%4i%s" % (line[0:22],num,line[26:]))
        else:
            out.append(line)

    return "".join(out)
