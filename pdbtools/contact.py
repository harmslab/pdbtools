#!/usr/bin/env python

# Copyright 2008, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ = \
"""
pdb_contacting-residues.py

Find all residues within some cutoff of a target residue in a pdb file.
"""
__author__ = "Michael J. Harms"
__date__ = "111116"

import os

def pdbContacting(pdb,target,cutoff,target_type="resname"):
    """
    """

    cutoff_sq = (cutoff)**2

    to_take = ["ATOM  ","HETATM"]
    all_coord = [[l[12:26],[float(l[30+8*i:38+8*i]) for i in range(3)]]
                 for l in pdb if l[0:6] in to_take]

    if target_type == "resname":
        target_list = [a for a in all_coord if a[0][5:8].strip() == target]
    else:
        target_list = [a for a in all_coord if int(a[0][10:14]) == target]

    out = []
    for t in target_list:
        contacts = []
        for a in all_coord:
            if sum([(a[1][i]-t[1][i])**2 for i in range(3)]) < cutoff_sq:

                # ignore self
                if t[0] == a[0]:
                    continue
                contacts.append(a[0][5:].strip())

        # Grab only unique contacts
        contacts = dict([(c,()) for c in contacts]).keys()
        out.append("%s\t%s\n" % (t[0],("\t".join(contacts))))

    return out 
