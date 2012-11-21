# Simple Makefile for satk programs.
#
# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  
#

F77 = /usr/bin/g77
F77_FLAGS = -w -O2
BINARIES = crdfrmt.out mkwij.out  prescan.out  satkelni.out  statacc.out 


all: $(BINARIES)

.f.out:
	${F77} ${F77_FLAGS} $*.f -o $*.out

clean:
	rm *.out -f
