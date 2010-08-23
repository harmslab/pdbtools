#!/usr/bin/env python

# Copyright 2009, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_oligomer.py

Grabs the name of every molecule in the file, then looks at the biological 
assembly and decides whether or not the asymmetric unit has the relevant 
assembly.  If it does not, the program spits out false.  (Some day, it will
spit out the relevant assembly...)

"""
__author__ = "Michael J. Harms"
__date__ = "090305"

import sys, os

def parseBioMat(pdb):
    """
    Parse a BIOMT pdb entry.  This entry contains the matrix required to 
    transform the atoms in the pdb file into the relevant biological 
    assembly.

    The output of this function is a list of two-field tuples (chains,
    matricies).  The matricies entry is a nested tuple of 3x4 matricies.
    """

    # Grab biomolec entries
    biomolec_entries = [l for l in pdb if l[0:10] == "REMARK 350"]
    biomolec = [i for i, l in enumerate(biomolec_entries)
                if l[11:22] == "BIOMOLECULE"]

    biomolec_output = []
    for i in biomolec:
        
        j = i + 1
        while biomolec_entries[j][34:40] != "CHAINS":
            j += 1

        # Find chains associated with transformation within this biomolecule
        chains = []
        while biomolec_entries[j][34:40] == "CHAINS":
            chains.extend(biomolec_entries[j][42:].split(","))
            j += 1
        chains = tuple([c.strip() for c in chains])

        # Create a list of all transformation matricies for this molecule. 
        # This is rather hacked, but functional...
        matrix_list =[]
        current_matrix = 0
        tmp_matrix = [[],[],[]]

        # Go through all entries in the matrix.  Note that "j" is set above
        while j < len(biomolec_entries) and \
                                    biomolec_entries[j][13:18] == "BIOMT":

            # If we are on the next matrix, record the old matrix and reset
            # the temproary matrix
            if (int(biomolec_entries[j][22]) - 1) != current_matrix:
                matrix_list.append(tuple([tuple(m) for m in tmp_matrix]))
                tmp_matrix = [[],[],[]]
                current_matrix += 1
            
            # Place this line of the matrix in the temporary matrix 
            line = biomolec_entries[j]
            mat_line = int(line[18]) - 1
            tmp_matrix[mat_line] = [float(e) for e in line[23:].split()]

            j += 1
            
        matrix_list.append(tuple([tuple(m) for m in tmp_matrix]))

        biomolec_output.append((chains,tuple(matrix_list))) 

    return tuple(biomolec_output)            
  
 
def findAllChains(pdb):
    """
    Create a dictionay that keys chains to their molecule type.  If no molecule
    is assigned to the chain in the header (e.g. for HETATM entries), the 
    residue type for the chain is used as the name.  If a HETATM chain is 
    heterogeneous (having CA and HOH, for example), one or the other name will
    be used.
    """
    
    # Parse COMPND part of pdb header 
    compounds = [l for l in pdb if l[0:6] == "COMPND"]
    unique_comp = [(i,l[18:].strip()) for i, l in enumerate(compounds)
                   if l[11:16] == "CHAIN"] 

    # Grab the molecule name from the COMPND MOLECULE entry
    molec_name = []
    j = 1
    for c in unique_comp:
        while compounds[c[0]-j][11:19] != "MOLECULE":
            j += 1
        molec_name.append(compounds[c[0]-j][21:].strip()[:-1])
 
    # Create a dictionary of chains to compounds
    chain_dict = {}
    for i in range(len(molec_name)):
        tmp_chains = [c.strip() for c in unique_comp[i][1][:-1].split(",")]
        chain_dict.update([(c,molec_name[i]) for c in tmp_chains])

    # grab chains from HETATM entries (if they were not listed in the COMPND
    # section)
    known_chains = chain_dict.keys()
    hetatm = [l for l in pdb if l[0:6] == "HETATM"]
    hetatm_chains = dict([(l[21],l[17:20]) for l in hetatm])
    for k in hetatm_chains.keys():
        if k not in known_chains:
            chain_dict[k] = hetatm_chains[k]
    
    return chain_dict 


def pdbOligomer(pdb,collapse_repeat=True):
    """
    Create a report about the oligomerization state of a pdb file.  If 
    collapse_repeat is specified, sequential repeating entries will be 
    reported only once.
    """

    chain_dict = findAllChains(pdb)

    # Determine if a transformation must be applied to the pdb, or if we have
    # a WYSIWYG oligomerization state.
    out = []
    wysiwyg = True
    biomat = parseBioMat(pdb)

    # If no biological matrix is present, create dummy output that will flag
    # the system as non wysiwg, but will allow chains, etc. to be parsed for
    # output
    if biomat == ():
        keys = chain_dict.keys()
        keys.sort()
        biomat = [[keys,(False,False)]]

    for biomolec in biomat:
      
        # If there is more than one matrix, there must be a real transformation
        if len(biomolec[1]) != 1:
            wysiwyg = False
        else:

            # If the transformation is not the identity matrix...
            if biomolec[1][0] != ((1.0,0.0,0.0,0.0),
                                  (0.0,1.0,0.0,0.0),
                                  (0.0,0.0,1.0,0.0)):
                wysiwyg = False

        # Createa a list of molecule names.  If collapse_repeat is specified,
        # a space (" ") is inserted for repeated molecule names.
        if collapse_repeat:
            molec_names = []
            tmp_molec_names = [chain_dict[c] for c in biomolec[0]]
            molec_names.append(tmp_molec_names[0])
            for i in range(1,len(tmp_molec_names)):
                if tmp_molec_names[i] != tmp_molec_names[i-1]:
                    molec_names.append(tmp_molec_names[i])
                else:
                    molec_names.append(" ")
        else:
            molec_names = [chain_dict[c] for c in biomolec[0]]
        
        # Create pretty output    
        out.extend(["%s;" % m for m in molec_names])
        out.append("\t")
        out.extend(["%s;" % c for c in biomolec[0]])
        out.append("\t%s\t|" % wysiwyg)
 
    return "".join(out)

def main():
    """
    Function to call if run from command line.
    """

    from helper import cmdline

    
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="c",
                          long_flag="collapse_repeat",
                          action="store_true",
                          default=True,
                          help="Collapse repeated molecule names")


    file_list, options = cmdline.parseCommandLine()

    out = []
    for pdb_file in file_list:
        
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        olig_output = pdbOligomer(pdb,options.collapse_repeat)
        
        pdb_id = pdb_file[:pdb_file.index(".pdb")]
        pdb_id = os.path.split(pdb_id)[-1] 

        out.append("%s\t%s" % (pdb_id,olig_output))

    out = ["%i\t%s\n" % (i,l) for i, l in enumerate(out)]
    out.insert(0,"\tpdb\tcompounds\tchains\tWYSIWYG\n")

    print "".join(out)


# If run from command line...
if __name__ == "__main__":
    main()
        

