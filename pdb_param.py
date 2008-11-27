#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_param.py:

Calculates some molecular parameters given the sequence of the protein in a pdb
file.
    0) The molecular weight of the protein
    1) The pI of the protein assuming model compound pKa values
    2) The fraction titratable, acidic, basic
    3) The charge on the molecule assuming all ASP, GLU, HIS, LYS, and ARG are
       ionized
"""

__author__ = "Michael J. Harms"
__date__ = "080125"


import os, sys
import pdb_seq
from pdb_data.common import *

class PdbParamError(Exception):
    """
    General error class for this module.
    """
    pass

def calcChargeState(count_dict,pH):
    """
    Calculates the charge state of a protein at a particular pH.
    """
    
    charge_state = 0.
    for aa in count_dict.keys():
        pKa = PKA_DICT[aa]
        charge = CHARGE_DICT[aa]
        charge_state += count_dict[aa]*charge/(1 + 10**(charge*(pH-pKa)))

    return charge_state


def pdbPi(aa_list,initial_pH_step=2.0,cutoff=0.001):
    """
    Calculate the PI of a sequence of amino acids.
    """

    # Count the number of titratable amino acids in a sequence
    count_dict = dict([(k,0) for k in PKA_DICT.keys()])
    for k in count_dict.keys():
        count_dict[k] = len([aa for aa in aa_list if aa == k])
    count_dict["NTERM"] = 1
    count_dict["CTERM"] = 1
 
    # Determine the pI using a simple convergence algorithm
    #
    #  |<---------------no
    #  |                 |
    #  |             ------------
    # pH += step --> |Q(pH) < 0? | <-----------------|
    #                ------------                    |
    #                    |                           |
    #                   yes --> step *= 0.5 --> (pH -= step)
    
    step = initial_pH_step
    pH = 0
    while step > cutoff:
        if calcChargeState(count_dict,pH) > 0:
            pH += step
        else:
            step *= 0.5
            pH -= step
   
    return pH

def calcMW(aa_list):
    """
    Calculate the molecular weight of a sequence.
    """
    
    # Convert list to molecular weight
    try:
        mw = sum([MW_DICT[aa] for aa in aa_list]) 
    except KeyError:
        err = "Sequence contains non-standard amino acids!"
        raise PdbParamError(err)

    # Subtract the molecular for every bond in the protein
    mw = mw - (len(aa_list) - 1)*MW_H2O

    return mw


def pdbParam(pdb,chain="all",use_atoms=False):
    """
    Calculate some general electrostatic properties of a protein.
    """

    # Use pdb_seq to grab the sequence of the protein(s)
    chain_dict, seq_type = pdb_seq.pdbSeq(pdb,use_atoms)
    
    # Convert chain dictionary to list of amino acids
    if chain == "all":
        aa_list = []
        for c in chain_dict.keys():
            aa_list.extend(chain_dict[c])
    else:
        if chain in chain_dict.keys():
            aa_list = chain_dict[chain]
        else:
            err = "Chain \"%s\" is not in pdb file!" % chain
            raise PdbParamError(err)
    
    # Count number of each type of group 
    count_dict = dict([(k,0) for k in MW_DICT.keys()])
    for k in count_dict.keys():
        count_dict[k] = len([aa for aa in aa_list if aa == k]) 
   
    # Calculate molecular weight and pI
    mw = calcMW(aa_list) 
    pI = pdbPi(aa_list)
   

    return count_dict, mw, pI, seq_type

    
def main():
    """
    If called from the command line, execute this.
    """

    from helper import cmdline
   
    # Parse command line
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="a",
                      long_flag="atomseq",
                      action="store_true",
                      default=False,
                      help="use ATOM sequence, not SEQRES")
    cmdline.addOption(short_flag="f",
                      long_flag="freq",
                      action="store_true",
                      default=False,
                      help="write out amino acid frequencies")
    cmdline.addOption(short_flag="c",
                      long_flag="chain",
                      action="store",
                      default="all",
                      help="chain to analyze",
                      nargs=1,
                      type=str)
    
    file_list, options = cmdline.parseCommandLine()

    out = []
    for pdb_index, pdb_file in enumerate(file_list):
        
        # Read in pdb file
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        # Calculate simple properties about protein charge
        try:
            count_dict, mw, pI, seq_type = pdbParam(pdb,chain=options.chain,
                                                    use_atoms=options.atomseq)
        except PdbParamError, (strerror):
            print "Error processing file \"%s\":" % pdb_file
            print "%s" % (strerror)
            sys.exit()

        aa_list = count_dict.keys()
        aa_list.sort()
        
        # Calculate fraction ionizable
        total = float(sum(count_dict.values()))
        titr_aa_list = [aa for aa in aa_list if aa in PKA_DICT.keys()]
        titr_total = float(sum([count_dict[aa] for aa in titr_aa_list]))

        fx_titr = titr_total/total

        # Calculate the apparent charge
        acid = sum([count_dict[aa] for aa in ["ASP","GLU"]])
        base = sum([count_dict[aa] for aa in ["HIS","LYS","ARG"]])
        app_charge = base - acid

        # Print to output in pretty way
        short_pdb = os.path.split(pdb_file)[-1][:-4]
        out.append("%30s%10s%10i%10.2F%10.2F%10i\n" % 
                   (short_pdb,seq_type.strip(),mw,pI,fx_titr,app_charge))
      
        # Write out amino acid frequencies if requested
        if options.freq:
            
            total_freq = dict([(aa,count_dict[aa]/total) for aa in aa_list])
            titr_freq = dict([(aa,count_dict[aa]/titr_total) 
                              for aa in titr_aa_list])

            freq_out = [5*"%10s" % (" ","aacid","counts","fx_total","fx_titr")]
            freq_out.append("\n")
            for aa_index, aa in enumerate(aa_list):
                if aa in titr_aa_list:
                    freq_out.append("%10i%10s%10i%10.2F%10.2F\n" % 
                        (aa_index,aa,count_dict[aa],100*total_freq[aa],
                        100*titr_freq[aa]))
                else:
                    freq_out.append("%10i%10s%10i%10.2F%10s\n" % 
                        (aa_index,aa,count_dict[aa],100*total_freq[aa],"NA"))

            g = open("%s_freq.txt" % (pdb_file[:-4]),"w")
            g.writelines(freq_out)
            g.close()

    out = ["%10i%s" % (i,l) for i, l in enumerate(out)]      
    out.insert(0,"%10s%30s%10s%10s%10s%10s%10s\n" % 
               (" ","pdb","seq","mw","pI","fx_titr","charge"))
    print "".join(out) 
   

if __name__ == "__main__":
    main()

