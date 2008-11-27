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

import sys

try:
    from numpy import array, linalg, dot, mean, squeeze
except ImportError:
    print "numpy package is required for this program!"
    print "(http://numpy.scipy.org/)"
    sys.exit()

class PdbCenterAsuError(Exception):
    """
    General error class for this module.
    """

    pass


def pdbCenterasu(pdb, write_coord =False):
    """
    Moves the center of mass of the assymmetric unit into the unit cell.
    Does not search for symmetry-generated equivalent asymmetric units
    """

    # get coords from file
    coords = array([[float(line[30+i*8:38+i*8])for i in range(3)] for line in pdb
             if line[0:6] == "ATOM  " or line[0:6] == "HETATM"])

    # get conversion matrices from pdb header
    fmx = []
    for line in pdb:
        if line.startswith('SCALE'):
            data = map(float, line[6:].split()[:-1])
            fmx.append(data)
    
    if len(fmx) == 0:
        err = "No SCALE records found in pdb file!\n"    
        raise PdbCenterAsuError(err)

    fmx = array(fmx)
    omx = linalg.inv(fmx)

    # main procedure
    fcoords = dot(coords, fmx.transpose())            # fractionalize coords
    fcntr = mean(fcoords, axis =0)
    new_fcntr = fcntr % 1
    fvector = new_fcntr - fcntr
    new_fcoords = fcoords + fvector
    new_coords = dot(new_fcoords, omx.transpose())    # back to orthogonal
    nc = new_coords

    # modify coordinates
    pdb_out = []
    coord = 0
    for line in pdb:
        if line[0:6] not in ('ATOM  ', 'HETATM'):
            pdb_out.append(line)
            continue
        line_out = [line[0:30],nc[coord][0],nc[coord][1],nc[coord][2],line[54:]]
        pdb_out.append("%s%8.3F%8.3F%8.3F%s" % tuple(line_out))
        coord += 1

    return fvector, pdb_out


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
                      help="write out PDB file with re-centered ASU")
    file_list, options = cmdline.parseCommandLine()

    for pdb_file in file_list:

        # Load in pdb file
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        translation, pdb_out = pdbCenterasu(pdb, options.coord)
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


