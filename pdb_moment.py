#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_moment.py

Calculate the approximate dipole moment of a protein.  This program needs polar
hydrogens to do its calculation.  The pdb file must be in the slightly non-
standard (UHBD) format written out by pdb_addH.  Calculation is done assuming 
that Asp, Glu, Arg, Lys, and His are fully charged.  Charges are placed on the
CD, CG, CZ, NZ, and NE2 atoms respectively (see pdb_data/polar_param.txt).
"""
__author__ = "Michael J. Harms"
__date__ = ""

import os, sys
from helper import geometry
from pdb_data.polar import POLAR_CHARGE_DICT
from math import pi

class PdbMomentError(Exception):
    """
    General error class for this module.
    """

    pass

def calcMoment(charges,coord):
    """
    Actually do moment calculation on centered coordinates.
    """

    # unit conversion
    elementary_q = 1.602176487e-19 # elementary charge in C
    debye_conv = 3.33564e-30       # 1 Debye in C*m
    k = elementary_q*1e-10/debye_conv  

    num_atoms = len(charges)

    # Calculate the moment by sum(q(i)*r(i)) for x, y, and z
    moment = [sum([coord[j][i]*charges[j] for j in range(num_atoms)])
            for i in range(3)]

    moment = [m*k for m in moment]

    return moment


def pdbMoment(pdb):
    """
    Calculate the dipole moment of a pdb file using the charge set in 
    pdb_data.polar.
    """

    # Parse pdb file
    pdb = [l for l in pdb if l[0:4] == "ATOM"]

    # Grab fractional charge of each atom and coordinates from pdb file
    try:
        charges = [POLAR_CHARGE_DICT[l[17:21]][l[12:16]] for l in pdb]
    except KeyError:
        err = l
        err += "Residue/atom pair not recognized"
        raise PdbMomentError(err)
    
    coord = [[float(line[30+i*8:38+i*8])for i in range(3)] for line in pdb]
    
    num_atoms = len(charges)

    # Center coordinates on center of mass
    center = [sum([c[i] for c in coord])/num_atoms for i in range(3)]
    coord = [[c[i] - center[i] for i in range(3)] for c in coord]

    return calcMoment(charges,coord)


def main():
    """
    Function to call if run from command line.
    """

    from helper import cmdline

    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="a",
                          long_flag="addH",
                          action="store_true",
                          default=False,
                          help="Add hydrogens automatically if missing.")


    file_list, options = cmdline.parseCommandLine()

    out = []
    for pdb_file in file_list:
        
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        short_pdb = os.path.split(pdb_file)[-1][:-4]

        # See if the file is in the proper (wacky) UHBD format
        ca = [l for l in pdb if l[0:4] == "ATOM" and l[12:15] == "CA "]
        if len(ca) == 0:

            # If not in the correct format, try to add hydrogens
            print "%s does not appear to have polar hydrogens!" % pdb_file
            if options.addH:
                
                import pdb_addH

                print "Attempting to add hydrogens with pdb_addH.py"
                try:
                    pdb = pdb_addH.pdbAddH(pdb,pdb_id=pdb_file[:-4])
       
                    g = open("%sH.pdb" % pdb_file[:-4],"w")
                    g.writelines(pdb)
                    g.close() 

                except pdb_addH.PdbAddHError, (strerror):
                    print "Addition of hydrogens failed for %s" % file
                    print strerror
                    sys.exit()
            else:
                print "Please add hydrogens to the file using pdb_addH.py"
                sys.exit()

        # Calculate the dipole moment and write out
        try:
            moment = pdbMoment(pdb)
        except PdbMomentError, (strerror):
            print "Problem with %s" % pdb_file
            print strerror
            sys.exit()
        out.append("%30s" % short_pdb)
        out.append("%10.3F%10.3F%10.3F" % tuple(moment))
        out.append("%10.3F\n" % geometry.dist(moment))

    # Combine output, split only on line breaks, remove blank lines
    out = "".join(out)
    out = out.split("\n")
    out = [x for x in out if x.strip() != ""]

    # Add line numbers and header
    out = ["%10i%s\n" % (i,x) for i, x in enumerate(out)]
    out.insert(0,"%10s%30s%10s%10s%10s%10s\n" %
               (" ","pdb","x","y","z","length"))
    
    print "".join(out)


# If run from command line...
if __name__ == "__main__":

    main()
        

