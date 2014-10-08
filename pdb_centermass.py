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

import sys, re
from pdb_data import common

def pdbCentermass(pdb,write_coord=False,include_hetatm=False,include_mass=True):
    """
    Calculates center of mass assuming the same mass for all atoms.  Either 
    returns recentered pdb or center of mass.
    """

    if include_hetatm:
        grab_lines = ["ATOM  ","HETATM"]
    else:
        grab_lines = ["ATOM  "] 

    # Create a regular expression to strip out charge information from atom type
    # Look for any number or + or -
    charge_pattern = re.compile("[0-9]|\+|\-")

    # Calculate the center of mass of the protein (assuming all atoms have the 
    # same mass).
    warn = False
    coord = []
    masses = []
    for l in pdb:
        
        # Skip non ATOM/HETATM lines
        if l[0:6] not in grab_lines:
            continue

        # Grab coordinates
        coord.append([float(l[30+i*8:38+i*8])for i in range(3)])

        # If we're ignoring atom mass, record all masses as 1
        if not include_mass:
            masses.append(1.0)
            continue

        # Otherwise, grab mass of each atom.  First try to grab atom type entry.
        # If it's missing, guess from the atom name column.
        atom_type = l[73:].strip()
        if atom_type == "":

            warn = True
            if l[12] == " ":
                atom_type = l[13]
            elif l[12] == "H":
                atom_type = "H"
            else:
                atom_type = l[12:14].strip()

        # Strip charge information from atom type
        atom_type = re.sub(charge_pattern,"",atom_type)

        # Now try to grab the mass 
        try: 
            masses.append(common.ATOM_WEIGHTS[atom_type])
        except:
            print "File contains atoms of unknown type (%s)" % atom_type
            print "Will assign them mass of carbon (12.011 g/mol)"
            print "To fix, edit ATOM_WEIGHTS dictionary in pdb_data/common.py"
            masses.append(12.011)

    num_atoms = len(coord)
    
    total_mass = sum(masses)
    weights = [m/total_mass for m in masses] 
    center = [sum([coord[i][j]*weights[i] for i in range(num_atoms)])
              for j in range(3)]

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

    if warn:
        print "Warning.  No element entries in file.  Attempting to extract"
        print "from the atom names.  Not always reliable..."

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
    cmdline.addOption(short_flag="m",
                      long_flag="mass",
                      action="store_false",
                      default=True,
                      help="weight atoms by their mass")
    cmdline.addOption(short_flag="a",
                      long_flag="hetatm",
                      action="store_true",
                      default=False,
                      help="include hetatms in center of mass calc")
    

    file_list, options = cmdline.parseCommandLine()
   
    for pdb_file in file_list:

        # Load in pdb file
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        center, pdb_out = pdbCentermass(pdb,options.coord,options.hetatm,
                                        options.mass)
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
    

