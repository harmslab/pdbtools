#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ =\
"""
pdb_offset.py

Adds an offset to the residue numbers of a pdb without touching anything else.
"""

__author__ = "Michael J. Harms"
__date__ = "070131"


import os, sys

def pdbOffset(pdb_file,offset):
    """
    Adds an offset to the residue column of a pdb file without touching anything
    else.
    """

    # Read in the pdb file
    f = open(pdb_file,'r')
    pdb = f.readlines()
    f.close()

    out = []
    for line in pdb:
        # For and ATOM record, update residue number
        if line[0:6] == "ATOM  " or line[0:6] == "TER   ":
            num = offset + int(line[22:26])
            out.append("%s%4i%s" % (line[0:22],num,line[26:]))
        else:
            out.append(line) 

    return "".join(out)

def main():
    """
    Main function to run if script is called from the command line.
    """

    from helper import cmdline
   
    # Parse command line
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="o",
                          long_flag="offset",
                          action="store",
                          default=0,
                          help="offset; num atoms to offset",
                          nargs=1,
                          type=int)
    file_list, options = cmdline.parseCommandLine()
    offset = options.offset


    # Offset pdb files in file_list
    out_files = []
    for f in file_list:
        out_files.append(pdbOffset(f,offset))

    # If it is a single file, dump to command line.  If it is a directory,
    # write each file to r%s % file.  
    if len(out_files) == 1:
        print out_files[0]
    else:
        for i, f in enumerate(file_list):
            out_name = list(os.path.split(f))
            out_name[1] = out_name[1][:-4]
            out_name = os.path.join(out_name[0],"%s_offset.pdb" % out_name[1])
            g = open(out_name,'w')
            g.write(out_files[i])
            g.close()

if __name__ == "__main__":
    main()



