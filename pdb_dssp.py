#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

"""
NOT FINISHED!!!"
"""

import sys

__author__ = "Michael J. Harms"
__date__ = "070604"
__version__ = "0.1"
__description__ = \
"""
pdb_dssp.py

Uses dssp to determine secondary structure assignment for each residue in a pdb
file.  This assignment is then placed in header HELIX/SHEET records.
"""

SS_RECORDS = ["HELIX ","SHEET ","TURN  "]
RECORD_ORDER = ["FORMUL","SSBOND","LINK  ","CISPEP","SITE  ","CRYST1"]

def runDssp(pdb):
    """
    Run DSSP on a pdb file.
    """

    pass

def parseDsspOutput(dssp_assign):
    """
    Parse DSSP output.
    """

    pass

def pdbDssp(pdb):
    
    dssp_assign = runDssp(pdb)
    ss_records = parseDsspOutput(dssp_assign)


def main():
    """
    Called if program is called from the command line.
    """

    print "This program is not finished!"
    sys.exit()

    import pdb_cmdline

    parser = pdb_cmdline.initializeParser(__description__,__version__)
    file_list, options = parser.parseCommandLine()

    for pdb_file in file_list:
        f = open(pdb_file,"r")
        pdb = f.readlines()
        f.close()

        pdb = pdbDssp(pdb)
   
    

if __name__ == "__main__":
    main()


