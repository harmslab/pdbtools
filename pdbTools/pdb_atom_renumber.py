#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_atom_renumber.py

Renumbers atom records in a pdb without touching residue numbers.
"""

__author__ = "Michael J. Harms"
__date__ = "070529"

import os, sys

def pdbAtomRenumber(pdb,renumber_het=True):
    """
    Renumber all atoms in pdb file, starting from 1.
    """

    entries_to_renumber = ["ATOM  ","TER   ","ANISOU"]
    if renumber_het == True:
        entries_to_renumber.append("HETATM")

    out = []
    counter = 1
    for line in pdb:
        # For and ATOM record, update residue number
        if line[0:6] in entries_to_renumber:
            out.append("%s%5s%s" % (line[0:6],counter,line[11:]))
            counter += 1
        else:
    
            # reset the counter for a new model
            if line[0:6] == "ENDMDL":
                counter = 1

            out.append(line)

    return out

def main():
    """
    Function to call if this is called from command line.
    """

    from helper import cmdline

    # Parse command line
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="r",
                      long_flag="renumber-het",
                      action="store_false",
                      default=True,
                      help="include hetatm entries in the renumbering.")

    file_list, options = cmdline.parseCommandLine()

    for pdb_file in file_list:
    
        # Read in the pdb file
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        out = pdbAtomRenumber(pdb,options.renumber_het)
        if len(file_list) == 1:
            print "".join(out)
        else:
            out_file = "%s_renum.pdb" % pdb_file[:-4]
            g = open(out_file,'w')
            g.writelines(out)
            g.close()

if __name__ == "__main__":
    main()

