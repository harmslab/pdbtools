# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
polar.py

Module with dictionaries for polar charge and radii.  Dictionaries are 
generated from polar_param.txt.
"""
__author__ = "Michael J. Harms"
__date__ = "070204"

import os

PARAM_FILE = os.path.join(os.path.split(__file__)[0],"polar_param.txt")

def _readParam(param_file):
    """
    Read a polar parameter file, returning charge and radii keyed to residue
    and atom type.
    """

    # Read in file
    f = open(param_file,'r')
    param = f.readlines()
    f.close()

    # Parse file
    param = [l for l in param if l[0] != "#" and l.strip() != ""]
    param = [(l[0:4],l[5:9],float(l[10:16]),float(l[30:37])) for l in param]

    # Fill charge and radii dictionaries keyed to aa, then atom type
    charge = dict([(x[0],[]) for x in param])
    radii  = dict([(x[0],[]) for x in param])
    for k in charge.keys():
        charge[k] = dict([(x[1],x[2]) for x in param if x[0] == k])
        radii[k]  = dict([(x[1],x[3]) for x in param if x[0] == k])

    return charge, radii


POLAR_CHARGE_DICT, POLAR_RADII_DICT = _readParam(PARAM_FILE)


