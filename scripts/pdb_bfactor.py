#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ = \
"""
pdb_bfactor.py:

Alters the bfactor column of a pdb file.  If a data file is supplied, it reads
residue number/value pairs and places the values in the b-factor column of a
pdb file.
"""
__author__ = "Michael J. Harms"
__date__ = "070706"

import sys
from pdbtools.helper import cmdline
from pdbtools import bfactor


def main():
    """
    Call if program called from command line.
    """

    from helper import cmdline

    # Parse command line
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="d",
                          long_flag="data_input",
                          action="store",
                          default=None,
                          help="FILE X Y; data from FILE columns X and Y",
                          nargs=3)
    cmdline.addOption(short_flag="a",
                          long_flag="abs_value",
                          default=False,
                          action="store_true",
                          help="Take absolute value of data in data file")
    cmdline.addOption(short_flag="s",
                      long_flag="set_value",
                      default=20.0,
                      action="store",
                      help="Single value to set the b-factors to.",
                      nargs=1,
                      type=float)

    file_list, options = cmdline.parseCommandLine()

    for pdb_file in file_list:
        # Read in file
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        # If a data file is specified, take the values from it and place in
        # data_dict.  Otherwise, set all residues to 20.0
        if options.data_input != None:

            # Finish processing command line options
            file = options.data_input[0]
            col1 = int(options.data_input[1])
            col2 = int(options.data_input[2])
            if options.abs_value:
                abs_value = True
            else:
                abs_value = False

            data_dict = bfactor.loadDataFile(file,col1,col2,abs_value)
        else:
            ca_list = [l for l in pdb if l[0:4] == "ATOM" and l[13:16] == "CA "]
            data_dict = dict([(ca[22:26].strip(),options.set_value)
                              for ca in ca_list])

        out = bfactor.pdbBfactor(pdb,data_dict)

        print "".join(out)

if __name__ == "__main__":
    main()
