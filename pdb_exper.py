#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_experiment.py

Extracts data about experimental procedure used to generate structure.
"""
__author__ = "Michael J. Harms"
__date__ = "080128"

import os, sys

def pdbExperiment(pdb):
    """
    Grab experimental parameters out of pdb file.
    """

    # Extract protein name and organism
    header = [l for l in pdb if l[0:6] == "HEADER"][0][6:50].strip()
    header = "\'%s\'" % header

    organism = [l for l in pdb if l[11:30] == "ORGANISM_SCIENTIFIC"]
    if len(organism) == 0:
        organism = "NA"
    else:
        organism = organism[0][31:].strip()
        organism = organism.strip(";")
        organism = "\'%s\'" % organism

    exp_type = "".join([l for l in pdb if l[0:6] == "EXPDTA"])
    if exp_type[10:13] == "NMR":
        exp_type = "NMR"
        resolution = "NA"
        r_value = "NA"
        r_free = "NA"
    elif exp_type[10:15] == "X-RAY":
        exp_type = "XRAY"

        remarks = [l for l in pdb if l[0:6] == "REMARK"]

        # Extract resolution
        try:
            resolution = [l for l in remarks
                          if l[13:34] == "RESOLUTION RANGE HIGH"][0]
            resolution = resolution.split(":")[1]
            resolution = "%10.2F" % float(resolution)
        except (IndexError,ValueError):
            resolution = "%10s" % "NA"


        fit_hash = [l[12:42] for l in remarks]
        try:
            fit_start = fit_hash.index("FIT TO DATA USED IN REFINEMENT")
        except ValueError:
            fit_start = fit_hash.index("USING DATA ABOVE SIGMA CUTOFF.")
        fit_end = fit_hash[fit_start:].index("                              ")
        fit_remarks = remarks[fit_start:fit_start+fit_end]
       
        # Extract R-value (note that we take the *last* possible R-value in case
        # the first R-value is working + test set) 
        try:
            r_value = [l for l in fit_remarks
                       if l[13:25] == "R VALUE     "][-1]
            r_value = r_value.split(":")[1]
            r_value = "%10.1F" % (100*float(r_value))
        except (IndexError,ValueError):
            r_value = "%10s" % "NA"
       
        # Extract R-free 
        try:
            r_free = [l for l in fit_remarks
                      if l[13:30] == "FREE R VALUE     "][0]
            r_free = r_free.split(":")[1]
            r_free = "%10.1F" % (100*float(r_free))
        except (IndexError,ValueError):
            r_free = "%10s" % "NA"

    exp_data = "%48s%40s%5s%s%s%s" % (header,organism,exp_type,resolution,r_value,r_free)

    return exp_data.lower()

def main():
    """
    Function to call if run from command line.
    """

    from helper import cmdline

    
    cmdline.initializeParser(__description__,__date__)

    file_list, options = cmdline.parseCommandLine()

    out = []
    for pdb_file in file_list:
        
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()

        exp_data = pdbExperiment(pdb)
        
        pdb_id = pdb_file[:pdb_file.index(".pdb")]
        pdb_id = os.path.split(pdb_id)[-1] 

        out.append("%30s%s" % (pdb_id,exp_data))

    out = ["%10i%s\n" % (i,x) for i, x in enumerate(out)]
    out.insert(0,"%10s%30s%48s%40s%5s%10s%10s%10s\n" % (" ","pdb","protein",
               "organism","exp","res","r_value","r_free"))
    
    print "".join(out)


# If run from command line...
if __name__ == "__main__":
    main()
        

