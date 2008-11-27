# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdbObject.py

A heirarchical object to represent structures stored in pdb files.  The
Structure class holds instances of the Chain class which holds instances of the
Fragments class.  The purpose of this module is to 1) find protein fragments,
regardless of numbering/chain scheme.  2) Create dictionaries that can flip the
protein sequence between "raw" (e.g. what is in the original pdb file) and
"fixed" (e.g. continuous numbering from 1 with no chain identifiers).
"""
__author__ = "Micahel J. Harms"
__date__ = "070725"

from geometry import dist
from math import sqrt

class PdbContainerError(Exception):
    """
    General error class for this module.
    """

    pass


class Fragment:
    """
    Holds a continuous fragment of amino acids and whether or not the fragment
    should have termini.
    """
    
    def __init__(self,atom_lines,n_terminus=False,c_terminus=False):
        """
        Initialize class.
        """
        
        self.atom_lines = atom_lines
        self.n_terminus = n_terminus
        self.c_terminus = c_terminus

    def _renumberAtoms(self,renumber_dict):
        """
        Renumber atoms in fragment according to renumber dictionary.
        """
        
        for i in range(len(self.atom_lines)):
            key = self.atom_lines[i][21:26]
            self.atom_lines[i] = "%s%s%s" % (self.atom_lines[i][:21],
                                             renumber_dict[key],
                                             self.atom_lines[i][26:])
        
    def write(self,filename=None):
        """
        Write out the fragment.
        """
        
        if filename != None:
            f = open(filename,'w')
            f.writelines(self.atom_lines)
            f.close()
        else:
            return "".join(self.atom_lines)

class Chain:
    """
    Hold a single molecule chain (as a set of fragments).
    """
    
    def __init__(self,chain_id,atom_lines,seq_lines,frag_dist_cutoff=5.0):
        """
        Hold set of fragments that represents the chain of atoms held in
        atom_lines.
        """
        
        self.chain_id = chain_id
        self.atom_lines = atom_lines
        self.seq_lines = seq_lines
        
        # Use comparison of atom_lines and seq_lines to determine whether or not
        # to add N- and C-termini to chain.  If seq_lines is empty, do not add
        # termini.
        if len(seq_lines) != 0:
            self.findTermini()
        else:
            self.n_terminus = False
            self.c_terminus = False
        self.findFragments(frag_dist_cutoff)        
        
        num_breaks = len(self.frag_index)
        # If there is no break in the chain, create a single fragment
        if num_breaks == 0:
            self.fragments = [Fragment(self.atom_lines,self.n_terminus,
                                       self.c_terminus)]
            
        # Otherwise, create correct number of fragments
        else:
            # Create n-terminal fragments
            self.fragments = [Fragment(self.atom_lines[:self.frag_index[0]],
                                       n_terminus=self.n_terminus)]
            
            # Create internal fragments
            for i in range(0,num_breaks-1):
                index1 = self.frag_index[i]
                index2 = self.frag_index[i+1]
                self.fragments.append(Fragment(self.atom_lines[index1:index2]))
            
            # Create c-terminal fragment
            self.fragments.append(Fragment(self.atom_lines[self.frag_index[-1]:],
                                           c_terminus=self.c_terminus))
    

    def findTermini(self,num_compare=10):
        """
        Determine whether or not to add n and c termini to this chain based on
        comparision of sequence in SEQRES entries and sequence in structure. I
        do not want to attempt a full alignment to determine if the termini are
        represented.  I'll just check the first and last $num_compare residues
        of gene_seq and struct_seq to make sure that they are the same to add
        termini.  This will fail in cases where:
          1) long sequences are exactly repeated throughout protein and this seq
             happens to align with start (not terribly likely)
          2) There is a gap in the structure within the first $num_compare
             residues.
        """

        # Default is not to add termini
        self.n_terminus = False
        self.c_terminus = False    

        # Determine the gene sequence using SEQRES entires.
        gene_seq = []
        for l in self.seq_lines:
            gene_seq.extend(l.split())
        
        # Determine the sequence in the structure using CA atoms
        struct_seq = [l[17:20] for l in self.atom_lines if l[13:16] == "CA "]
    
        # Determine if termini should be added.
        if struct_seq[:num_compare] == gene_seq[:num_compare]:
            self.n_terminus = True
        if struct_seq[num_compare:] == gene_seq[num_compare:]:
            self.c_terminus = True

    def findFragments(self,dist_cutoff=5.0):
        """
        Take lines of pdb in self.atom_lines and find breaks in chain.
        """
        
        # For every residue in the pdb file, grab the CA and N lines.  This
        # assumes that there are no duplicate residues and that all residues
        # have CA and N.  If there is a duplicate residue or missing atom, an
        # error is raised.

        residue_list = []
        for line in self.atom_lines:
            if line[21:26] not in residue_list:
                residue_list.append(line[21:26])

        ca_list = []
        n_list = []        
        for resid in residue_list:
            resid_atoms = [l for l in self.atom_lines if l[21:26] == resid]

            ca = [l for l in resid_atoms if l[13:16] == "CA "]
            n = [l for l in resid_atoms if l[13:16] == "N  "]
            if len(ca) > 1 or len(n) > 1:
                err = "Residue \"%s\" is duplicated!" % resid
                raise PdbContainerError(err)
            elif len(ca) == 0 or len(n) == 0:
                err = "Residue \"%s\" has missing CA or N atoms!" % resid
                raise PdbContainerError(err)

            ca_list.append(ca[0])
            n_list.append(n[0])
   
        # Grab coordinates of CA atoms and index of residues indexes of 
        # amide nitrogens.
        ca_coord = [[float(l[30+8*i:38+8*i]) for i in range(3)] for l in ca_list]
        res_index = [self.atom_lines.index(n) for n in n_list]
            
        # Check the distance between every ca carbon and the next.  If these
        # are further apart than dist_cutoff, place in different fragments.
        self.frag_index = []
        for i in range(1,len(ca_list)):
            if dist(ca_coord[i-1],ca_coord[i]) > dist_cutoff:
                self.frag_index.append(res_index[i])

    def _renumberAtoms(self,renumber_dict):
        """
        Renumber fragments in fragment list according to renumber dictionary.
        Use renumbered fragments to renumber self.atom_list.
        """

        self.atom_lines = []
        for fragment in self.fragments:
            fragment._renumberAtoms(renumber_dict)
            self.atom_lines.extend(fragment.atom_lines)
        

    def write(self,filename=None):
        """
        Write out the Chain.
        """
        
        out = []
        for frag in self.fragments:
            out.append(frag.write())

        if filename != None:
            f = open(filename,'w')
            f.writelines(out)
            f.close()
        else:
            return "".join(out)
        

class Structure:
    """
    Class to hold a pdb file (broken into chains, then into fragments).
    """
    
    def __init__(self,pdb_id,seq_lines,atom_lines,frag_dist_cutoff=5.0):
        """
        Read pdb file and find chains in the pdb file.
        """
        
        self.pdb_id = pdb_id
        self.seq_lines = seq_lines
        self.atom_lines = atom_lines
  
        # Create self.chains (list containing Chain objects for each chain in
        # the pdb file).
        self.findChains(frag_dist_cutoff)
        
        # Create dictionaries to convert between raw numbering scheme and
        # fixed numbering scheme
        self.residue_list = []
        for line in self.atom_lines:
            if line[21:26] not in self.residue_list:
                self.residue_list.append(line[21:26])
        
        self.raw2fix = dict([(r,"%5i" % (i+1))
                             for i, r in enumerate(self.residue_list)])
        self.fix2raw = dict([("%5i" % (i+1),r)
                             for i, r in enumerate(self.residue_list)])
        
        self.current_numbering = "raw"

   
    def findChains(self,frag_dist_cutoff=5.0):
        """
        Find the indiviual chains in a pdb file given all atom entries in the
        file.
        """
        
        # Create list of chains and dictionary of atoms in chain
        self.chain_list = []
        self.chain_atoms = {}

        for line in self.atom_lines:
            chain = line[21]
            if chain not in self.chain_list:
                self.chain_list.append(chain)
                self.chain_atoms.update([(chain,[])])
            self.chain_atoms[chain].append(line)    

                
        # Populate self.chain_seq (from SEQRES) keyed to chain_id
        self.chain_seq = dict([(k,[]) for k in self.chain_list])        
        for line in self.seq_lines:
            
            # This try/except statement skips SEQRES entries that do not have
            # atom entries in the file.
            try:
                self.chain_seq[line[11]].append(line[19:70])
            except KeyError:
                pass           

        # Create self.chains, which holds Chain object for each chain
        self.chains = []
        for chain_id in self.chain_list:
            self.chains.append(Chain(chain_id,self.chain_atoms[chain_id],
                                     self.chain_seq[chain_id],
                                     frag_dist_cutoff))

    def renumberAtoms(self):
        """
        Renumber all atoms.  Toggles between "fixed" and "raw" states.  The
        "fixed" state has all residues in the protein renumbered from 1, and no
        chain identifiers.  The "raw" state has whatever numbering was
        originally assigned in the pdb file.
        """ 
        
        # Toggle status between "fixed" and "raw" numbering.
        if self.current_numbering == "fixed":
            renumber_dict = self.fix2raw
            self.current_numbering = "raw"
        elif self.current_numbering == "raw":
            renumber_dict = self.raw2fix
            self.current_numbering = "fixed"
        
        # Renumber all atoms in chains
        self.atom_lines = []
        for chain in self.chains:
            chain._renumberAtoms(renumber_dict)
            self.atom_lines.extend(chain.atom_lines)


    def write(self,filename=None):
        """
        Write out the entire structure.
        """
        
        out = []
        for chain in self.chains:
            out.append(chain.write())
            if self.current_numbering == "raw": out.append("%-79s\n" % "TER")
        out.append("%-80s" % "END")

        if filename != None:
            f = open(filename,'w')
            f.writelines(out)
            f.close()
        else:
            return "".join(out)

    def dumpStructures(self):
        """
        Dump list of all fragments in structure with termini information.
        """
        
        structure_list = [] 
        for chain in self.chains:
            for frag in chain.fragments:
                structure_list.append([frag.atom_lines,frag.n_terminus,
                                       frag.c_terminus])
        
        return structure_list

    def dumpNumberConversion(self,output_file):
        """
        Write dictionaries used to convert number schemes to a file.
        """
        
        to_write = ["\"%s\":\"%s\"\n" % (k,self.raw2fix[k]) for k in self.raw2fix]
        to_write.sort()
        to_write.insert(0,"# raw:fixed\n")

        f = open(output_file,'w')
        f.writelines(to_write)
        f.close()
    
    def loadNumberConversion(self,input_file,current_numbering):
        """
        Take file written by dumpNumberConversion and load into raw2fix and
        fix2raw.  current_numbering must be specified to indicate the current
        state of the pdb relative to the loaded dictionaries.*** NOTE!  THIS
        WILL WIPE OUT THE AUTOMATICALLY GENERATED CONVERSION DICTIONARIES! ***
        """
        
        f = open(input_file,'r')
        input = f.readlines()
        f.close()
      
        # Parse the file; note lines beginning with # are treated as comments 
        input = [l for l in input if l[0] != "#"]
        input = [(l.split(":")[0],l.split(":")[1].strip("\n")) for l in input]
        input = [(x[0].strip("\""),x[1].strip("\"")) for x in input]

        # Load numbering into instance dictionaries. 
        self.current_numbering = current_numbering
        self.raw2fix = dict(input) 
        self.fix2raw = dict([(x[1],x[0]) for x in input])


def main():
    """
    Main function (for testing purposes).  Program is not meant to be called.
    """

    import sys
    
    input_pdb = sys.argv[1]

    f = open(input_pdb)
    pdb_lines = f.readlines()
    f.close()

    atom_lines = [l for l in pdb_lines if l[0:4] == "ATOM"]
    seq_lines = [l for l in pdb_lines if l[0:6] == "SEQRES"]

    x = Structure(input_pdb[:-4],seq_lines,atom_lines)
    for c in x.chains:
        for i, f in enumerate(c.fragments):
            f.write("%s_%s_%s.pdb" % (x.pdb_id,c.chain_id,i))

    x.renumberAtoms()
    x.write("%s_fixed.pdb" % x.pdb_id)

    x.renumberAtoms()
    x.write("%s_recovered.pdb" % x.pdb_id)

    x.dumpNumberConversion("junk.conversion")

if __name__ == "__main__":
    main()

