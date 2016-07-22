#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_subset.py

Takes a subset of residues/chains from a pdb file.  All non-atom entries are 
left intact.
"""

__author__ = "Michael J. Harms"
__date__ = "080121"



class PdbSubsetError(Exception):
    """
    General error class for this module.
    """

    pass

def pdbSubset(pdb,chain,residues,
              atom_list=["ATOM  ","ANISOU","HETATM","TER   "]):
    """
    Grabs a subset of a pdb file, spitting out the pdb, the chains, and the 
    actual min and max in the pdb file.
    """
    
    if chain == "all":
        chains = ()
        chain_string = "all" 
    else:
        chains = tuple(chain.split(","))
        chain_string = "".join(chains)

    seen_chains = []
    out = []
    for line in pdb:

        # only look at records indicated by atom_list
        if line[0:6] not in atom_list:
            out.append(line)
            continue    
        
        # Grab only residues belonging to specified chains.  If chains is empty
        # take every chain
        if chains == () or line[21:22] in chains:

            if line[21:22] not in seen_chains:
                seen_chains.append(line[21:22])
   
            # Grab only residues above minimum and below maximum 
            if residues[0] == 0 and residues[1] == 0:
                out.append(line)
            elif int(line[22:26]) >= residues[0] and \
                        int(line[22:26]) <= residues[1]:
                out.append(line)
            else:
                continue
        else:
            continue

    if len([l for l in out if l[0:6] in atom_list]) == 0:
        err = "No residues in pdb meet subset criteria.\n"
        raise PdbSubsetError(err)

    seen_chains.sort()
    input_chains = list(chains)
    input_chains.sort()

    if len(input_chains) > 0:
        if seen_chains != input_chains:
            missing_chains = ",".join([c for c in input_chains
                                       if c not in seen_chains])
            err = "pdb file did not have chain(s): %s.\n" % missing_chains
            raise PdbSubsetError(err)
        

    # Determine the actual min and max in the pdb file
    residues = [0,0]
    all_residues = [int(l[22:26]) for l in out if l[0:6] in atom_list]
    residues[0] = min(all_residues)
    residues[1] = max(all_residues)                                 

    return out, chain_string, residues


def main():
    """
    If called from command line.
    """

    import sys
    from helper import cmdline
   
    # Parse command line
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="c",
                          long_flag="chain",
                          action="store",
                          default="all",
                          help="chain to select (separte mulitiple by commas)",
                          nargs=1,
                          type=str)
    cmdline.addOption(short_flag="r",
                          long_flag="residues",
                          default=[0,0],
                          action="store",
                          help="residues to select",
                          nargs=2,
                          type=int)

    file_list, options = cmdline.parseCommandLine()

    for pdb_file in file_list:

        # Attempt to strip .pdb extension
        try:
            pdb_id = pdb_file[:pdb_file.index(".pdb")]
        except IndexError:
            pdb_id = pdb_file

        # Read file
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        # Pop out subset
        pdb, chain, residues = pdbSubset(pdb,options.chain,options.residues)
      
        # Write to file 
        out_file = "%s_%s_%i-%i.pdb" % (pdb_id,chain,residues[0],residues[1])
        g = open(out_file,'w')
        g.writelines(pdb)
        g.close()



if __name__ == "__main__":
    main()

