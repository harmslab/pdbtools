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
