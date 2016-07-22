#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# Copyright 2008, Marcin Cieslik
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ = \
"""
pdb_centerasu.py

Moves the assymmetric unit into the unit cell.
"""

__author__ = "Michael J. Harms and Marcin Cieslik"
__date__ = "080408"

from pdbtools.helper import cmdline
from pdbtools import centerasu

def main():
    """
    If called from command line...
    """

    # Parse command line
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="c",
                      long_flag="coord",
                      action="store_true",
                      default=False,
                      help="write out PDB file with re-centered ASU")
    file_list, options = cmdline.parseCommandLine()

    for pdb_file in file_list:

        # Load in pdb file
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        translation, pdb_out = centerasu.pdbCenterasu(pdb, options.coord)
        print "%s %s" % (pdb_file, tuple(translation.ravel()))
        # If the user wants re-centered PDB files, write them out
        if options.coord:
            out_file = "%s_asucenter.pdb" % pdb_file[:-4]
            g = open(out_file,'w')
            g.writelines(pdb_out)
            g.close()


# If called from command line:
if __name__ == "__main__":
    main()
