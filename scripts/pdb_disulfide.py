#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ = \
"""
pdb_disulfide.py

Returns a list of cys in disulfide bonds in a protein using a simple distance
constraint.
"""

__author__ = "Michael J. Harms"
__date__ = "070726"

from pdbtools.helper import cmdline
from pdbtools import disulfide

def main():
    """
    Main function, if called from command line.
    """
    # Parse command line
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="c",
                      long_flag="cutoff",
                      action="store",
                      default=4,
                      help="cutoff for disulfide bonds",
                      nargs=1,
                      type=float)


    file_list, options = cmdline.parseCommandLine()

    for pdb_file in file_list:

        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        free_cys, disulfide_bonds = disulfide.pdbDisulfide(pdb,options.cutoff)

        out = ["<<%s>>\nFree Cys:        " % pdb_file]
        out.extend(["%s, " % c.strip() for c in free_cys])
        out.append("\nDisulfide Bonds: ")
        out.extend(["%s, " % d.strip() for d in disulfide_bonds])

        print "".join(out)

if __name__ == "__main__":
    main()
