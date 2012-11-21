#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_contact.py

Generate a distance-based contact plot for a protein.
"""

__author__ = "Michael J. Harms"
__date__ = ""

import os, sys
from helper import geometry

def pdbContact(pdb,all=False):
    """
    Generate array of all atom-atom distances in a pdb file.  Defaults to only 
    look at CA atoms
    """

    # Remove non-CA atoms
    if not all:
        pdb = [l for l in pdb if l[0:4] == "ATOM" and l[13:16] == "CA "]
    else:
        pdb = [l for l in pdb if l[0:4] == "ATOM"]

    coord = [[float(l[30+8*i:38+8*i]) for i in range(3)] for l in pdb]

    dist = geometry.calcDistances(coord)

    return dist    

def main():
    """
    Function to call if run from command line.
    """

    from helper import cmdline

    
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="a",
                          long_flag="all",
                          action="store_true",
                          default=False,
                          help="Use all atoms (not just CA)")
    cmdline.addOption(short_flag="d",
                          long_flag="difference",
                          action="store",
                          default=None,
                          help="Take difference between two contact maps. " +\
                               "DIFFERENCE is a pdb file",
                          nargs=1,
                          type=str)


    file_list, options = cmdline.parseCommandLine()

    # If we are calculating a difference map, calculate this first.
    if options.difference != None:
        if not os.path.isfile(options.difference):
            print "\"%s\" is not a file!" % options.difference
            sys.exit()
        
        f = open(options.difference,'r')
        pdb = f.readlines()
        f.close()

        to_compare = pdbContact(pdb,options.all)


    for pdb_file in file_list:
        
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        dist = pdbContact(pdb,options.all)

        # Make sure the difference and main pdbs have the same number of atoms
        if options.difference != None:
            if len(dist) != len(to_compare):
                print "\"%s\" and \"%s\" do not have same number of atoms!" % \
                    (pdb_file,options.difference)
                sys.exit()

        # Write output
        out = []
        out.append("# Distance contact map for \"%s\".\n" % pdb_file)
        if options.difference != None:
            out.append("# Difference map to \"%s\".\n" % options.difference)
            out.append("#%9s%10s%10s%10s%10s%10s%10s\n" % 
                       ("atom1","atom2","d12","d12_r","dd12","mean","d/m"))
        else:
            out.append("#%9s%10s%10s\n" %  ("atom1","atom2","d12"))
                        
        for i in range(len(dist)):
            for j in range(len(dist)):
                out.append("%10i%10i%10.3F" % (i,j,dist[i][j]))

                if options.difference != None:
                    d = dist[i][j] - to_compare[i][j]
                    mean = (dist[i][j] + to_compare[i][j])/2.
                    if mean != 0:
                        d_m = d/mean
                    else:
                        d_m = 0

                    out.append("%10.3F%10.3F%10.3F%10.3F\n" % 
                               (to_compare[i][j],d,mean,d_m))
                else:
                    out.append("\n")
                               
            out.append("\n")
        out = "".join(out)
                
        short_pdb = os.path.split(pdb_file)[-1][:-4]
        g = open("%s_contact.txt" % short_pdb,'w')
        g.write(out)
        g.close()


# If run from command line...
if __name__ == "__main__":
    main()
        

