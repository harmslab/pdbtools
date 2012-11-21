# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
common.py

Common data for all scripts in pdb_tools.
"""
__author__ = "Michael J. Harms"
__date__ = "080123"

# --------------------------------------------------------------------------- #
# Chemical data
# --------------------------------------------------------------------------- #

# Dictionary of pKa values and unprotonated charge state.
PKA_DICT = {"ASP":     4.0,
            "CYS":     8.5,
            "GLU":     4.4,
            "TYR":    10.0,
            "CTERM":   3.1,
            "ARG":    12.0,
            "HIS":     6.5,
            "LYS":    10.4,
            "NTERM":   8.0}

CHARGE_DICT = {"ASP":   -1.,
               "CYS":   -1.,
               "GLU":   -1.,
               "TYR":   -1.,
               "CTERM": -1.,
               "ARG":    1.,
               "HIS":    1.,
               "LYS":    1.,
               "NTERM":  1.}

# Atom on which to place charge
TITR_ATOM = {"ASP":"CG ",
             "GLU":"CD ",
             "TYR":"OH ",
             "ARG":"CZ ",
             "HIS":"NE2",
             "LYS":"NZ "}

# Dictionary of amino acid molecular weights.  The the molecular weight of water
# should be subtracted for each peptide bond to calculate a protein molecular
# weight.
MW_DICT = {"ALA" : 89.09,
           "ARG" : 174.20,
           "ASN" : 132.12,
           "ASP" : 133.10,
           "ASX" : 132.61,
           "CYS" : 121.15,
           "GLN" : 146.15,
           "GLU" : 147.13,
           "GLX" : 146.64,
           "GLY" : 75.07,
           "HIS" : 155.16,
           "ILE" : 131.17,
           "LEU" : 131.17,
           "LYS" : 146.19,
           "MET" : 149.21,
           "PHE" : 165.19,
           "PRO" : 115.13,
           "SER" : 105.09,
           "THR" : 119.12,
           "TRP" : 204.23,
           "TYR" : 181.19,
           "VAL" : 117.15,
           "IVA" : 85.12,
           "STA" : 157.15,
           "ACE" : 43.30}
MW_H2O = 18.0

# Dictionary of van der Waal radii
VDW_DICT = {"N  ":1.650,
            "CA ":1.870,
            "C  ":1.760,
            "O  ":1.400,
            "AD1":1.700,
            "AD2":1.700,
            "AE1":1.700,
            "AE2":1.700,
            "CB ":1.870,
            "CD ":1.870,
            "CD1":1.760,
            "CD2":1.760,
            "CG ":1.870,
            "CG1":1.870,
            "CG2":1.870,
            "CE ":1.870,
            "CE1":1.760,
            "CE2":1.760,
            "CE3":1.760,
            "CH2":1.760,
            "CZ ":1.870,
            "CZ2":1.760,
            "CZ3":1.760,
            "ND ":1.650,
            "ND1":1.650,
            "ND2":1.650,
            "NE ":1.650,
            "NE1":1.650,
            "NE2":1.650,
            "NH1":1.650,
            "NH2":1.650,
            "NZ ":1.650,
            "OD ":1.400,
            "OD1":1.400,
            "OD2":1.400,
            "OT2":1.400,
            "OE ":1.400,
            "OE1":1.400,
            "OE2":1.400,
            "OH2":1.400,
            "OG ":1.400,
            "OG1":1.400,
            "OT1":1.400,
            "OXT":1.400,
            "OH ":1.400,
            "SD ":1.850,
            "SG ":1.850,
            "P  ":1.900}


# --------------------------------------------------------------------------- #
# Amino acid name data
# --------------------------------------------------------------------------- #

# List of three and one letter amino acid codes
_aa_index = [('ALA','A'),
             ('CYS','C'),
             ('ASP','D'),
             ('GLU','E'),
             ('PHE','F'),
             ('GLY','G'),
             ('HIS','H'),
             ('HSE','H'),
             ('HSD','H'),
             ('ILE','I'),
             ('LYS','K'),
             ('LEU','L'),
             ('MET','M'),
             ('MSE','M'),
             ('ASN','N'),
             ('PRO','P'),
             ('GLN','Q'),
             ('ARG','R'),
             ('SER','S'),
             ('THR','T'),
             ('VAL','V'),
             ('TRP','W'),
             ('TYR','Y')]

AA3_TO_AA1 = dict(_aa_index)
AA1_TO_AA3 = dict([(aa[1],aa[0]) for aa in _aa_index])

# --------------------------------------------------------------------------- #
# PDB record data
# --------------------------------------------------------------------------- #

# Deprecated records
DEPRECATED_RECORDS = ["TURN  ","HYDBND","SLTBRG"]

# Types of coordinate entries
COORD_RECORDS = ["ANISOU","ATOM  ","HETATM","END   ","ENDMDL","MODEL","TER   "]

# Problems with structure warranting user attention
ERROR_RECORDS = ["CAVEAT","OBSLTE"]


