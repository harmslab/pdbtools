#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ = \
"""
pdb_ion-dist.py

Calculates the distribution of distances between titratable residues (from
titratable atoms in pdb_data.common) for acid/acid, acid/base, and base/base
interactions.
"""
__author__ = "Michael J. Harms"
__date__ = "080303"

from .data.common import CHARGE_DICT, TITR_ATOM
from .helper.geometry import calcDistances

def pdbIonDist(pdb,hist_step,remove_resid="TYR"):
    """
    Calculate ij distances between all titratable atoms in pdb (except residues
    in skip_resid) then bin according to hist_step.  The output histogram is
    a 3 x N nested list where N is the length required to cover all distances
    in pdb using hist_step and the 3 overarching lists hold acid/acid,
    acid/base, and base/base interactions in 0, 1, and 2 respectively.
    """

    pdb = [l for l in pdb if l[0:4] == "ATOM"]
    titr_atom = [l for l in pdb
                 if l[17:20] not in remove_resid and
                    l[17:20] in list(TITR_ATOM.keys()) and
                    l[13:16] == TITR_ATOM[l[17:20]]]

    # Extract coordinates and charges for each titratable atom
    coord = [[float(l[30+8*i:38+8*i]) for i in range(3)] for l in titr_atom]
    charge = [CHARGE_DICT[l[17:20]] for l in titr_atom]

    # Calculate all ij distances
    dist = calcDistances(coord)
    N = len(coord)

    # Initialize histogram by finding maximum distance that will have to be
    # counted.
    max_dist_bin = max([max([dist[i][j] for j in range(N)])
                        for i in range(N)])
    max_dist_bin = int(round(max_dist_bin/hist_step)) + 1
    histogram = [[0 for j in range(max_dist_bin)] for i in range(3)]

    # Populate histogram
    #   Interaction types:
    #   acid/acid            --> 0 (-1 + -1)/2 + 1
    #   acid/base, base/acid --> 1 (-1 +  1)/2 + 1
    #   base/base            --> 2 ( 1 +  1)/2 + 2
    for i in range(N):
        for j in range(i+1,N):
            interaction_type = int((charge[i] + charge[j])/2 + 1)
            dist_bin = int(round(dist[i][j]/hist_step))
            histogram[interaction_type][dist_bin] += 1

    return histogram
