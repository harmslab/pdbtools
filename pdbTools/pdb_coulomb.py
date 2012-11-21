#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_coulomb.py

Calculates the total coulomb energy (kcal/mol) of a protein structure in a pdb
file assuming that all groups titrate with model compound pKa values.
"""

__author__ = "Michael J. Harms"
__date__ = "070520"

import sys, os, copy
from math import sqrt, exp
from pdb_data.common import *

def readPDB(pdb_file):
    """
    Takes a pdb file and reads in the coordinates of each titratable group.  
    Assigns pka and charge state of each group.
    """


    # Open pdb_file and read each line into pdb (a list of lines)
    f = open(pdb_file,'r')
    pdb = f.readlines()
    f.close()

    # Grab only ATOM entries that are titratable
    pdb = [l for l in pdb if l[0:4] == "ATOM" and 
                             l[17:20] in TITR_ATOM.keys() and 
                             l[13:16] == TITR_ATOM[l[17:20]]]

    # Initalize lists to hold coordinates, pkas, and charge
    coord, pKa, charge = [], [], []

    # Go through each line in the pdb file
    for line in pdb:
        amino_acid = line[17:20]

        # Grab the xyz coordinates
        coord.append([float(line[30+8*i:38+8*i]) for i in range(3)])

        pKa.append(PKA_DICT[amino_acid])

        # Look up charge 
        charge.append(CHARGE_DICT[amino_acid])

    # Return the coordinates, pka, and charge
    return coord, pKa, charge

def hendersonHasselbach(pKa,charge,pH):
    """
    Calculate the fractional charge on a group with pKa and charge at some
    pH value.
    """

    return charge/(1 + 10**(charge*(pH-pKa)))
 

def pdbCoulomb(coord,pKa,charge,dielec_const,ionic_str,pH,temperature):
    """
    Calculates the energy of a structure given the coordinates of each 
    charged atom, their fractional charge, the dielectric constant, and the 
    ionic strength.
    """

    ionic_str = ionic_str/1000

    # Initialize variables
    kappa = 50.29*sqrt(ionic_str/(dielec_const*temperature))
    num_groups = len(coord)
    energy = 0.

    hh_chg = [hendersonHasselbach(pKa[i],charge[i],pH) 
              for i in range(num_groups)]

    # Calculate energy of interaction of every ij interaction (making sure not 
    # to double count; note we start j at i + 1).
    for i in range(num_groups):
        for j in range(i+1,num_groups):

            # Calculate distance between atom i and atom j
            r = sqrt(sum([(coord[i][k]-coord[j][k])**2 for k in range(3)]))

            # Add the energy of this interaction to the total energy
            energy += 332*hh_chg[i]*hh_chg[j]/(r*dielec_const)*exp(-kappa*r)

    # Return energy
    return energy


def main():
    """
    If called from command line, calculate energy and print to standard out.
    """

    from helper import cmdline
   
    # Parse command line
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="d",
                          long_flag="dielec_const",
                          action="store",
                          default=40.0,
                          help="dielectric constant",
                          nargs=1,
                          type=float)
    cmdline.addOption(short_flag="i",
                          long_flag="ionic_str",
                          default=100,
                          action="store",
                          help="ionic strength (mM)",
                          type=float)
    cmdline.addOption(short_flag="p",
                          long_flag="pH",
                          default=7.0,
                          action="store",
                          help="pH",
                          type=float)
    cmdline.addOption(short_flag="T",
                          long_flag="temperature",
                          default=298.0,
                          action="store",
                          help="temperature (K)",
                          type=float)
    cmdline.addOption(short_flag="t",
                          long_flag="titrate",
                          default=(None,None),
                          action="store",
                          help="titrate a variable",
                          type=str,
                          nargs=2)

    file_list, options = cmdline.parseCommandLine()
    
    # create dictionary of option values where values are in lists
    value_dict = dict([(k,[options.__dict__[k]])
                       for k in options.__dict__.keys() if k != "titrate"])

    # Deal with whether the user has specified a text file containing values
    # over which to titrate.
    if options.titrate != (None,None):
        available_options = value_dict.keys()

        # Make sure that the specified option can titrate and the file exists. 
        titration = options.titrate[0]
        data_file = options.titrate[1]
        if titration in available_options:
            if os.path.isfile(data_file):
                f = open(data_file,'r')
                titr_data = f.readlines()
                f.close()

                # Strip comments and blank lines, then re-join all lines
                titr_data = [l for l in titr_data
                             if l[0] != "#" and l.strip() != ""]
                titr_data = "".join(titr_data)
                
                # Parse file
                try:
                    titr = [float(x) for x in titr_data.split()]
                    value_dict[titration] = titr[:]

                # Do some basic error checking
                except ValueError:
                    print "Data file \"%s\" has mangled data!" % data_file
                    sys.exit()
                
                if len(titr) == 0:
                    print "Data file \"%s\" is empty!" % data_file
                    sys.exit()
            else:
                print "Data file \"%s\" does not exist!" % data_file
                sys.exit()
        else:
            print "\"%s\" cannot be titrated!" % titration
            print "Available titrations:"
            for option in available_options:
                print "\t%s" % option
            sys.exit()


    out = ["%10s%30s%10s%10s%10s%10s%10s\n" % 
           (" ","pdb","ep","ion_str","pH","T","dG")]
    counter = 0
    for pdb_file in file_list:

        # Read in coordinates
        coord, pKa, charge = readPDB(pdb_file)
        short_pdb = os.path.split(pdb_file)[-1][:-4]

        # Determine electrostatic energy, titrating over all relavent variables
        for temperature in value_dict["temperature"]:
            for pH in value_dict["pH"]:
                for ionic_str in value_dict["ionic_str"]:
                    for dielec in value_dict["dielec_const"]:
                        
                        energy = pdbCoulomb(coord,pKa,charge,dielec,ionic_str,
                                            pH,temperature)

    
                        # Write out results
                        out.append("%10i%30s%10.3F%10.3F%10.3F%10.3F%10.3F\n" % 
                                   (counter,short_pdb,dielec,ionic_str,
                                    pH,temperature,energy))
                        counter += 1

    print "".join(out)

if __name__ == "__main__":
    main()

