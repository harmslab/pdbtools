#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdb_bfactor.py:

Alters the bfactor column of a pdb file.  If a data file is supplied, it reads
residue number/value pairs and places the values in the b-factor column of a
pdb file.
"""
__author__ = "Michael J. Harms"
__date__ = "070706"

import sys

def loadDataFile(data_file,col1=0,col2=1,abs_value=False):
    """
    Parses input from a two column data file, creating a dictionary of
    column[0]:column[1] key/value pairs.  Returns dictionary.
    """
    
    # Read in file
    f = open(data_file,'r')
    data = f.readlines()
    f.close()

    # Strip out blank lines and comments
    data = [d for d in data if d[0] != "#" and d.strip() != ""]

    data_dict = {}
    for record in data:
        try:
            field = record.split()
            key = field[col1]
            data_dict.update([(key,float(field[col2]))])

            if abs_value:
                data_dict[key] = abs(data_dict[key])
            
        except IndexError:
            sys.stderr.write("Mangled data, skipping line:\n%s" % record)
            continue
        except ValueError:
            sys.stderr.write("Mangled data, skipping line:\n%s" % record)
            continue
        
    return data_dict


def pdbBfactor(pdb,data_dict):
    """
    Goes through pdb line by line.  If residue is in dictionary data_dict,
    the b-factor is replaced in output_file by the dictionary value.  If the
    residue is not in the dictionary, the b-factor is given value 0.0.
    Returns void.
    """

    out = []
    for line in pdb:
        if line[0:6] == "ATOM  ":
            resnum = line[23:26].strip()
            if resnum in data_dict.keys():
                out.append("%s%6.2F%s" % (line[:60],data_dict[resnum],
                                          line[66:]))
            else:
                out.append("%s%6s%s" % (line[:60],"NA",line[66:]))
        
        elif line[0:6] == "HETATM":
            out.append("%s%6s%s" % (line[:60],"NA",line[66:]))
        
        else:
            out.append(line)
    
    return out


def main():
    """
    Call if program called from command line.
    """

    from helper import cmdline
   
    # Parse command line
    cmdline.initializeParser(__description__,__date__)
    cmdline.addOption(short_flag="d",
                          long_flag="data_input",
                          action="store",
                          default=None,
                          help="FILE X Y; data from FILE columns X and Y",
                          nargs=3)
    cmdline.addOption(short_flag="a",
                          long_flag="abs_value",
                          default=False,
                          action="store_true",
                          help="Take absolute value of data in data file")
    cmdline.addOption(short_flag="s",
                      long_flag="set_value",
                      default=20.0,
                      action="store",
                      help="Single value to set the b-factors to.",
                      nargs=1,
                      type=float)

    file_list, options = cmdline.parseCommandLine()

    for pdb_file in file_list:
        # Read in file
        f = open(pdb_file,'r')
        pdb = f.readlines()
        f.close()
    
        # If a data file is specified, take the values from it and place in 
        # data_dict.  Otherwise, set all residues to 20.0
        if options.data_input != None:
       
            # Finish processing command line options
            file = options.data_input[0]
            col1 = int(options.data_input[1])
            col2 = int(options.data_input[2])
            if options.abs_value:
                abs_value = True
            else:
                abs_value = False

            data_dict = loadDataFile(file,col1,col2,abs_value)
        else:
            ca_list = [l for l in pdb if l[0:4] == "ATOM" and l[13:16] == "CA "]
            data_dict = dict([(ca[22:26].strip(),options.set_value)
                              for ca in ca_list])

        out = pdbBfactor(pdb,data_dict)

        print "".join(out)

if __name__ == "__main__":
    main()


