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
    ligands = dict([(l,[]) for l in ligands if l not in BORING_LIGANDS]).keys()

    return ";".join(ligands)

def main():
    """
    Function to call if run from command line.
    """

    from helper import cmdline

    
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="s",
                          long_flag="skip_boring",
                          action="store_true",
                          default=False,
                    help="Ignore boring ligands (from BORING_LIGANDS list)")

    file_list, options = cmdline.parseCommandLine()

    out = []
    for pdb_file in file_list:
        
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        pdb_id = pdb_file[:pdb_file.index(".pdb")]
        pdb_id = os.path.split(pdb_id)[-1] 
        
        ligand_out = pdbLigand(pdb,options.skip_boring)

        out.append("%s\t%s" % (pdb_id,ligand_out))

    out = ["%i\t%s\n" % (i,l) for i, l in enumerate(out)]

    print "".join(out)


# If run from command line...
if __name__ == "__main__":
    main()
        

