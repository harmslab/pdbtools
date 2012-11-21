#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_centermass.py

Calculates the center of mass of a protein (assuming all atoms have equal mass).
Returns either the center or the recentered coordinates of the pdb (if coord 
command line option is specified).
"""

__author__ = "Michael J. Harms"
__date__ = "061109"

import sys

def pdbCentermass(pdb,write_coord=False):
    """
    Calculates center of mass assuming the same mass for all atoms.  Either 
    returns recentered pdb or center of mass.
    """

    # Calculate the center of mass of the protein (assuming all atoms have the 
    # same mass).
    coord = [[float(line[30+i*8:38+i*8])for i in range(3)] for line in pdb
             if line[0:6] == "ATOM  "]
    num_atoms = len(coord)
    center = [sum([c[i] for c in coord])/num_atoms for i in range(3)]

    out = []
    # Re-center pdb file
    if write_coord:
        for line in pdb:
            if line[0:6] in ["ATOM  ","HETATM"]:
                line_out = [line[0:30],0.,0.,0.,line[54:]]
                for i in range(3):
                    line_out[i+1] = float(line[30+i*8:38+i*8]) - center[i]
                out.append("%s%8.3F%8.3F%8.3F%s" % tuple(line_out))
            else:
                out.append(line)

    center_out = ["%10.4F" % c for c in center]

    return center_out, out

def main():
    """
    If called from command line...
    """

    from helper import cmdline

    # Parse command line
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="c",
                      long_flag="coord",
                      action="store_true",
                      default=False,
                      help="write out re-centered coordinates")
    

    file_list, options = cmdline.parseCommandLine()
   
    for pdb_file in file_list:

        # Load in pdb file
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        center, pdb_out = pdbCentermass(pdb,options.coord)
        print "%s %s" % (pdb_file,"".join(center))

        # If the user wants re-centered coordinates, write them out
        if options.coord:
            out_file = "%s_center.pdb" % pdb_file[:-4]
            g = open(out_file,'w')
            g.writelines(pdb_out)
            g.close()


# If called from command line:
if __name__ == "__main__":
    main()
    

