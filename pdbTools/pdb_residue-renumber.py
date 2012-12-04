#!/usr/bin/env python

# Copyright 2007, 2012, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_residue-renumber.py

A tool for renumbering the residues of a pdb file.  It allows the user to specify
a new starting residue number, then propagates that numbering down the file.  
The user can specify whether or not to modify hetatms, whether to remove gaps 
(that is, make the numbering continuous), and whether to modify a particular 
chain.
"""

__author__ = "Michael J. Harms"
__date__ = "121204"

import os, sys

def pdbResidueRenumber(pdb,start_res=None,renumber_het=True,remove_gaps=False,
                       chain=None):
    """
    Renumber the residues in a pdb file.
    """

    # Initialize output lists
    fmt = "%-8s%6s%6i%10s\n"
    align_out = []
    pdb_out = []

    # Renumber het atoms
    if renumber_het:
        to_renumber = ["ATOM  ","HETATM"]
    else:
        to_renumber = ["ATOM  "]

    seen_chains = {}
    for line in pdb:
        # For and ATOM record, update residue number
        if line[0:6] in to_renumber:

            # If this is the first residue in a chain, set up counter etc.
            if line[21] not in seen_chains.keys():

                # If we're only looking at a specific chain (and this isn't it)
                # don't modify th eline and move along.
                if chain != None and line[21] != chain:
                    pdb_out.append(line)
                    continue

                current_res = line[22:26]
                if start_res != None:
                    output_res = start_res
                    offset = (start_res - int(current_res))
                else:
                    output_res = int(current_res)
                    offset = 0

                align_out.append(fmt % (line[17:22],current_res,output_res,
                                        int(current_res)==output_res))

                seen_chains[line[21]] = [output_res,offset]

            # If this is a new residue, update the counter
            if line[22:26] != current_res:
                current_res = line[22:26]
                if remove_gaps:
                    seen_chains[line[21]][0] += 1
                    output_res = seen_chains[line[21]][0]
                else:
                    output_res = int(current_res) + seen_chains[line[21]][1]


                align_out.append(fmt % (line[17:22],current_res,output_res,
                                        int(current_res)==output_res))

            pdb_out.append("%s%4i%s" % (line[0:22],output_res,line[26:]))

        # If this is a terminii record, reset the renumbering tool
        elif line[0:3] == "TER":
            if chain != None and line[21] != chain:
                pdb_out.append(line)
            else:
                pdb_out.append("%s%4i%s" % (line[0:22],output_res,line[26:]))

        # Reset the chains for a new model entry
        elif line[0:5] == "MODEL":
            seen_chains = {}
            pdb_out.append(line)

        # Print every other record without modification
        else:
            pdb_out.append(line)

    # Insert REMARK about renumber
    try:
        remark_index = [l[0:6] for l in pdb].index("REMARK")
    except ValueError:
        remark_index = 0
    pdb_out.insert(remark_index,"%-80s\n" %
                   ("REMARK    PDB RENUMBERED USING pdb_residue-renumber.py"))
   
    return pdb_out, align_out

def main():
    """
    Function to call if called from the command line.
    """

    from helper import cmdline

    # Parse command line
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="a",
                      long_flag="align-out",
                      action="store_true",
                      default=False,
                      help="write out residue alignment file.")
    cmdline.addOption(short_flag="s",
                      long_flag="start-res",
                      action="store",
                      type=int,
                      default=None,
                      help="Make the first residue in the file have s")
    cmdline.addOption(short_flag="r",
                      long_flag="renumber-het",
                      action="store_true",
                      default=True,
                      help="include hetatm entries in the renumbering.")
    cmdline.addOption(short_flag="g",
                      long_flag="remove-gaps",
                      action="store_true",
                      default=False,
                      help="make the sequence numbering continuous (remove gaps)")
    cmdline.addOption(short_flag="c",
                      long_flag="chain",
                      action="store",
                      type=str,
                      default=None,
                      help="only renumber this chain")

    
    file_list, options = cmdline.parseCommandLine()

    for pdb_file in file_list:
        
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        pdb_out, alignment = pdbResidueRenumber(pdb,options.start_res,
                                                options.renumber_het,
                                                options.remove_gaps,
                                                options.chain)

        short_pdb = os.path.split(pdb_file)[-1][:-4]
        g = open("%s_res-renum.pdb" % short_pdb,"w")
        g.writelines(pdb_out)
        g.close()

        if options.align_out:
            print "Residue alignment written out."

            alignment.insert(0,"#%s" % pdb_out[0])
            g = open("%s_res-algn.out" % short_pdb,'w')
            g.writelines(alignment)
            g.close()

if __name__ == "__main__":
    main()


