#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ = \
"""
pdb_splitnmr.py

Takes an nmr pdb file and splits each model into its own pdb file.
"""
__author__ = "Michael J. Harms"
__date__ = "080204"

import sys, os

def splitNMR(pdb):
    """
    Split each model in an NMR pdb file into its own pdb.
    """

    to_strip = ["ENDMDL","MASTER"]
    pdb = [l for l in pdb if l[0:6] not in to_strip]
    pdb_splitter = [(l[0:6],i) for i, l in enumerate(pdb)]
    pdb_splitter = [x[1] for x in pdb_splitter if x[0] == "MODEL "]
    model_numbers = [pdb[i].split()[1].strip() for i in pdb_splitter]
    pdb_splitter.append(len(pdb))

    all_models = []
    for i in range(1,len(pdb_splitter)):
        all_models.append((model_numbers[i-1],
                           pdb[pdb_splitter[i-1]:pdb_splitter[i]]))

    return all_models
