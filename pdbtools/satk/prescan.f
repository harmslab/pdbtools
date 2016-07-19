C      PROGRAM PRESCAN (INFILE1, OUTFILE1) 
C...   Program call commented out in i686 Linux port (MJH, 5/3/05)
C***********************************************************************
C...
C...   prescan.f is a program that identifies potential charge bearing
C...   atoms in the lists of coordinates of proteins and nucleic acids
C...   produced by crdfrmt.f. A new file is created that contains the
C...   coordinates for the charge bearing atoms.
C...
C...   The DATA statement is based on the appendix B and appendix F of
C...   the Protein Data Bank Atomic Coordinate and Bibliographic Entry
C...   Format Description published by Brookhaven in January 1985. It
C...   identifies explicitly all of the 20 amino acids, the 5 most
C...   frequently found nucleotides and a whole bunch of chemical-
C...   modified ones. The statement is already 19 lines long, which is
C...   the longest continuation on most Fortran compilers. If new data
C...   is entered the appropiate changes must accompany. Changes in the
C...   values of the upper limits of the DO 600, DO 610 and DO 620 must
C...   also be entered as necessary. The Fe and propionate residues of
C...   the heme group are handled explicitly. All other heteroatoms
C...   are copied textually unto the new file.
C...
C...   Bertrand Garcia-Moreno E.
C...
C...   Department of Biology         
C...   The Johns Hopkins University
C...   
C...   HP 9000-835SRX
C...
C...   Fortran 77
C...
C*********************************************************************

C        Copyright 1988, Bertrand Garcia-Moreno E.
C        This program is distributed under General Public License v. 3.  See the
C        file COPYING for a copy of the license.  
      
      CHARACTER ATNAM*5, RESNAM*5, HEADER*80, RLST*3, ALST*4, RLST1*3,
     *FILNAM1*20, FILNAM2*20, INFILE1*20, OUTFILE1*20
      INTEGER ATNUM, RESNUM, FLAG1, FLAG2
      REAL X, Y, Z
      DIMENSION RLST(30),ALST(25),RLST1(13)

      DATA RLST(1),RLST(2),RLST(3),RLST(4),RLST(5),RLST(6),RLST(7),
     *RLST(8),RLST(9),RLST(10),RLST(11),RLST(12),RLST(13),RLST(14),
     *RLST(15),RLST(16),RLST(17),RLST(18),RLST(19),RLST(20),RLST(21),
     *RLST(22),RLST(23),RLST(24),RLST(25),RLST(26),RLST(27),RLST(28),
     *RLST(29),RLST(30),RLST1(1),RLST1(2),RLST1(3),RLST1(4),RLST1(5),
     *RLST1(6),RLST1(7),RLST1(8),RLST1(9),RLST1(10),RLST1(11),
     *RLST1(12),RLST1(13),ALST(1),ALST(2),ALST(3),
     *ALST(4),ALST(5),ALST(6),ALST(7),ALST(8),ALST(9),ALST(10),ALST(11),
     *ALST(12),ALST(13),ALST(14),ALST(15),ALST(16),ALST(17),ALST(18),
     *ALST(19),ALST(20),ALST(21),ALST(22),ALST(23),ALST(24),ALST(25)
     */'GLU','ASP','ARG','LYS','HIS','TYR','CYS','HEM','ASX',
     *'GLX','UNK','  A',
     *'  C','  G','  T','  U','  I','1MA','5MC','OMC','1MG','2MG',
     *'M2G','7MG','OMG',' YG',' +U','H2U','5MU','PSU','ALA','ASN',
     *'GLN','GLY','ILE','LEU','MET','PHE','PRO','SER','THR',
     *'TRP','VAL',' NH2',' NH1',
     *' NZ ',' ND1',' NE2',' AD1',' AE2',' OE1',' OE2',' AE1',' AE2',
     *' OD1',' OD2',' AD1',' AD2',' OH ',' SG ',' O1A',' O2A',' O1D',
     *' O2D',' OXT','FE ',' O1P',' O2P'/

      FLAG1=0
      FLAG2=0

C...  Following 2 lines of code added by Mike Harms, 12/09/2004 
      CALL GETARG(1,INFILE1)
      CALL GETARG(2,OUTFILE1)

      OPEN (UNIT=10, FILE=INFILE1, STATUS='OLD', IOSTAT=IOS)
      OPEN (UNIT=20, FILE=OUTFILE1, STATUS='NEW', IOSTAT=IOS)

      DO 10 I=1,5
      READ(10,'(A80)') HEADER
      WRITE(20,'(A80)') HEADER
10    CONTINUE

C***********************************************************************
C...
C...  Charged atoms are found and written unto TAPE20. 
C...
C**********************************************************************

501   CONTINUE 
      READ(10,9000,END=9999) RESNUM, RESNAM, ATNAM, X, Y, Z, ATNUM
      IF((RESNAM.EQ.'     ').OR.(ATNAM.EQ.'     ')) GO TO 500

C***********************************************************************
C...
C...    The amino termini of polypeptides is dealt with
C...
C***********************************************************************

        IF((RESNUM.EQ.1).AND.(ATNAM(1:4).EQ.' N  ')) THEN
         WRITE(20,9000) RESNUM, RESNAM, ATNAM, X, Y, Z, ATNUM
         FLAG1=0

C***********************************************************************
C...
C...    The 5' end of polynucleotides are dealt with
C...
C***********************************************************************

        ELSE IF((RESNUM.EQ.1).AND.(ATNAM(1:4).EQ.' OXT')) THEN
         WRITE(20,9000) RESNUM, RESNAM, ATNAM, X, Y, Z, ATNUM

C***********************************************************************
C...
C...    The carboxyl termini of polypeptides is deal with
C...
C***********************************************************************

        ELSE IF((ATNAM(1:4).EQ.' OXT').AND.(FLAG1.EQ.0)) THEN
         WRITE(20,9000) RESNUM, RESNAM, ATNAM, X, Y, Z, ATNUM
         BACKSPACE 10
           DO 200 I=1,200
           READ(10,9000,END=9999) RESNUM, RESNAM, ATNAM, X, Y, Z, ATNUM
           BACKSPACE 10
           BACKSPACE 10
           IF((RESNAM.EQ.'     ').OR.(ATNAM.EQ.'     ')) GO TO 500
             IF (ATNAM(1:4).EQ.' O  ') THEN
              WRITE(20,9000) RESNUM, RESNAM, ATNAM, X, Y, Z, ATNUM
              FLAG1=1
              FLAG2=1
              GO TO 500
             ENDIF
200        CONTINUE
        ELSE IF((ATNAM(1:4).EQ.' OXT').AND.(FLAG2.EQ.1)) THEN
         FLAG2=0 

C***********************************************************************
C...
C...   Everything else is dealt with. If the residue type (aminoacid
C...   or nucleotide) is listed in the DATA statement, the record is
C...   'caught' in the DO 600 loop. If the atom type is one of the
C...   ones listed in the DATA statement, the record is 'caught' in
C...   the DO 610 loop. If the atom and residue in question is not
C...   one of the charged atoms listed in the DATA statement, then
C...   control of the program passes on to the DO 620, where if the
C...   residue type is not identified as one of the uncharged residues
C...   listed in the DATA statement, it is entered into the new file.
C...
C***********************************************************************

        ELSE IF (FLAG2.EQ.0) THEN
          DO 600 III=1,30
            IF(RESNAM(1:3).EQ.RLST(III)) THEN
             DO 610 JJJ=1,25
               IF(ATNAM(1:4).EQ.ALST(JJJ))THEN
                WRITE(20,9000) RESNUM, RESNAM, ATNAM, X, Y, Z, ATNUM
                GO TO 500
               ENDIF
610          CONTINUE
             GO TO 500
            ENDIF
600       CONTINUE
          DO 620 III=1,13
          IF(RESNAM(1:3).EQ.RLST1(III)) GO TO 500
620       CONTINUE
          WRITE(20,9000) RESNUM, RESNAM, ATNAM, X, Y, Z, ATNUM
        ENDIF
500   GO TO 501

C***********************************************************************
C...
C...   T H E     E N D
C...
C***********************************************************************

9000  FORMAT(5X,I4,3X,A5,3X,A5,2X,3(2X,F8.3),6X,I5)
9999   CONTINUE

       CLOSE(10)
       CLOSE(20)

       END
