C      PROGRAM CRDFRMT (INFILE1,OUTFILE1,OUTFILE2)
C...   Program call commented out in i686 Linux port (MJH, 5/3/05)
C***********************************************************************
C...
C...     This program is used to change the format of the coordinates
C...     of macromolecules. The files read from the Brookhaven Protein
C...     Data Bank are stripped of comments and remarks, connectivi-
C...     ties, etc, etc. In fact the only records that are retained 
C...     include a header identifying the structure, the ATOM records
C...     containing the X, Y, Z coordinates for the atoms, and the 
C...     HETATOM records containing the coordinates for all the hete-
C...     roatoms seen in the crystallographic structure. Water molecules
C...     are treated as special cases. The coordinates of these
C...     molecules are written unto TAPE30. The reformatted coordinates
C...     of the macromolecules are written unto TAPE20.
C...     When going from .pdb format to the .cor format, 5 lines of 
C...     header are read, followed by the ATOM information.
C...    
C...
C...     Bertrand Garcia-Moreno E.
C...     Department of Biology
C...     The Johns Hopkins University
C...     Baltimore, MD 21218        (301) 338-7239
C...
C...     Fortran 77, HP 9000 835-SRX
C...
C***********************************************************************

C        Copyright 1988, Bertrand Garcia-Moreno E.
C        This program is distributed under General Public License v. 3.  See the
C        file COPYING for a copy of the license.  

       CHARACTER HEADER*80, LABEL*6, ATNAM*5, TATNAM*5, RESNAM*5,
     * TRESNAM*5, TLABEL*6, INFILE1*20, OUTFILE1*20, OUTFILE2*20
       INTEGER ATNUM, TATNUM, RESNUM, TRESNUM, DUMYRES
       REAL X, TX, Y, TY, Z, TZ 

C...   Following 3 commands added my Mike Harms, 12/09/2004      
       CALL GETARG(1,INFILE1)
       CALL GETARG(2,OUTFILE1)
       CALL GETARG(3,OUTFILE2)

       OPEN(UNIT=10,FILE=INFILE1,STATUS='OLD')
       OPEN(UNIT=20,FILE=OUTFILE1,STATUS='NEW')
       IF(OUTFILE2(1:1).NE.' ') OPEN(UNIT=30,FILE=OUTFILE2,STATUS='NEW')

       DO 100 I=1,4
       READ(10,'(A80)') HEADER
       WRITE(20,'(A80)') HEADER
       IF(OUTFILE2(1:1).NE.' ') WRITE(30,'(A80)') HEADER
100    CONTINUE

       WRITE(20,'(A1)') ' '
       IF(OUTFILE2(1:1).NE.' ') WRITE(30,'(A1)') ' '


C***********************************************************************
C...
C...   The DO 500 loop that follows covers the entire length of the 
C...   program. The Protein Data Bank Atomic Coordinate and Bilbiogra-
C...   phic Entry Format Description published by PDB on January
C...   1985 was used to design the code inside the DO 500 loop.
C...
C***********************************************************************

       IF (OUTFILE2(1:1).NE.' ') THEN 

500    CONTINUE          
       READ(10,'(A6)',END=9999) LABEL
         IF ((LABEL.EQ.'ATOM  ').OR.(LABEL.EQ.'HETATM')) THEN
            BACKSPACE 10
            READ(10,'(A6,I5,1X,A5,A5,I4,4X,3F8.3)',END=9999) LABEL,
     *      ATNUM, ATNAM, RESNAM, RESNUM, X, Y, Z
            WRITE(20,'(5X,I4,3X,A5,3X,A5,2X,3(2X,F8.3),6X,I5)')
     *      RESNUM, RESNAM, ATNAM, X, Y, Z, ATNUM
              DO 510 II=1,100000
              IF(II.EQ.99999) PRINT *, 'Beware of limits in DO 510'
              READ(10,'(A6,I5,1X,A5,A5,I4,4X,3F8.3)',END=9999) LABEL,
     *        TATNUM, TATNAM, TRESNAM, TRESNUM, TX, TY, TZ
              IF((TATNAM(2:2).EQ.'D').OR.(TATNAM(2:2).EQ.'H')) GO TO 510
                IF(LABEL.EQ.'ATOM  ') THEN
                IF(TLABEL.EQ.'HETATM') WRITE (20,'(A1)') ' '
                TLABEL='ATOM  '
                  IF(TRESNUM.NE.RESNUM) THEN
                    WRITE(20,'(A1)')' '
                    WRITE(20,'(5X,I4,3X,A5,3X,A5,2X,3(2X,F8.3),6X,I5)')
     *              TRESNUM, TRESNAM, TATNAM, TX, TY, TZ, TATNUM
                    RESNUM=TRESNUM
                    RESNAM=TRESNAM
                     ATNAM=TATNAM
                     ATNUM=TATNUM
                         X=TX
                         Y=TY
                         Z=TZ
                    GO TO 510
                  ELSE IF (TRESNUM.EQ.RESNUM) THEN
                    WRITE(20,'(5X,I4,3X,A5,3X,A5,2X,3(2X,F8.3),6X,I5)')
     *              TRESNUM, TRESNAM, TATNAM, TX, TY, TZ, TATNUM 
                    RESNUM=TRESNUM
                    RESNAM=TRESNAM
                     ATNAM=TATNAM
                     ATNUM=TATNUM
                         X=TX
                         Y=TY
                         Z=TZ
                    GO TO 510
                  ENDIF
                ELSE IF (LABEL.EQ.'TER   ') THEN
                  GO TO 510
                ELSE IF (LABEL.EQ.'HETATM') THEN
  

C**********************************************************************
C...                                                              
C...    The coordinates of all the heteroatoms are found and appended
C...    to the end of the data file or to the end of the listing of
C...    atoms of individual chains or subunits. Waters are treated
C...    separately and written unto TAPE30
C...                                                         
C**********************************************************************

          IF((RESNUM.GE.100).AND.(RESNUM.LT.200)) DUMYRES=200       
          IF((RESNUM.GE.200).AND.(RESNUM.LT.300)) DUMYRES=300      
          IF((RESNUM.GE.300).AND.(RESNUM.LT.400)) DUMYRES=400     
          IF((RESNUM.GE.400).AND.(RESNUM.LT.500)) DUMYRES=500    
          IF((RESNUM.GE.500).AND.(RESNUM.LT.600)) DUMYRES=600   
          IF((RESNUM.GE.600).AND.(RESNUM.LT.700)) DUMYRES=700  
          IF((RESNUM.GE.700).AND.(RESNUM.LT.800)) DUMYRES=800 
          IF((RESNUM.GE.800).AND.(RESNUM.LT.900)) DUMYRES=900
          IF((RESNUM.GE.900).AND.(RESNUM.LT.1000)) DUMYRES=1000  
215         IF(TRESNAM.EQ.'HOH') THEN                       
              GO TO 220                                    
            ELSE IF (TRESNAM.EQ.'H2O') THEN               
              GO TO 220                                  
            ELSE IF (TRESNAM.EQ.'WAT') THEN             
              GO TO 220                                
            ELSE IF (TRESNAM.EQ.'OH2') THEN           
              GO TO 220                              
            ELSE                                  
              GO TO 230 
            ENDIF                               

C***********************************************************************
C...                                                            
C...    The statements that follow deal with the water molecules.
C...                                                       
C***********************************************************************

220    CONTINUE
       L=1
       WRITE(30,'(5X,I4,3X,A5,3X,A5,2X,3(2X,F8.3),6X,I5)') L,
     * TRESNAM, TATNAM, TX, TY, TZ, L
         DO 520 L=2,10000
         IF(L.EQ.9999) PRINT *, 'Beware of DO 520 loop'
         READ(10,'(A6,I5,1X,A5,A5,I4,4X,3F8.3)',END=9999) LABEL,
     *   ATNUM, ATNAM, RESNAM, RESNUM, X, Y, Z 
         IF(TRESNAM.NE.RESNAM)THEN
           BACKSPACE 10
           GO TO 510
         ELSE
           WRITE(30,'(A1)') ' '
           WRITE(30,'(5X,I4,3X,A5,3X,A5,2X,3(2X,F8.3),6X,I5)') L,
     *     RESNAM, ATNAM, X, Y, Z, L
         ENDIF
520      CONTINUE
       GO TO 510
 
C**********************************************************************
C...
C...     The following lines deal with heteroatoms.
C...
C**********************************************************************

230   CONTINUE
      WRITE(20,'(A1)') ' '
      WRITE(20,'(5X,I4,3X,A5,3X,A5,2X,3(2X,F8.3),6X,I5)') DUMYRES,
     *TRESNAM, TATNAM, TX, TY, TZ, TATNUM
        DO 250 L=1,200
        IF (L.EQ.199) PRINT *, 'Beware of DO 250 loop'
        READ(10,'(A6,I5,1X,A5,A5,I4,4X,3F8.3)',END=9999) LABEL,
     *  ATNUM, ATNAM, RESNAM, RESNUM, X, Y, Z 
          IF(RESNAM.NE.TRESNAM) THEN
          BACKSPACE 10
          TLABEL='HETATM'  
          GO TO 510
          ELSE    
          WRITE(20,'(5X,I4,3X,A5,3X,A5,2X,3(2X,F8.3),6X,I5)')
     *    DUMYRES, RESNAM, ATNAM, X, Y, Z, ATNUM
          ENDIF
250     CONTINUE
      WRITE (20,'(A1)') ' '
      GO TO 510

C***********************************************************************
C...
C...   Next section does .cor to .pdb conversion.
C...
C...
C**********************************************************************

       ENDIF
510    CONTINUE
       ENDIF
       GO TO 500

       ELSE IF(OUTFILE2(1:1).EQ.' ') THEN
       ICOUNT = 1
600    CONTINUE

       READ(10,'(5X,I4,3X,A5,3X,A5,2X,3(2X,F8.3),6X,I5)',END=9999)
     *      RESNUM, RESNAM, ATNAM, X, Y, Z, ATNUM
       IF(RESNAM.EQ.'     '.AND.ATNAM.EQ.'     ') THEN
         GO TO 600
       ELSE
         WRITE(20,'("ATOM  ",I5,1X,A5,A5,I4,4X,3F8.3)') ICOUNT, ATNAM,
     *         RESNAM, RESNUM, X, Y, Z
         ICOUNT = ICOUNT + 1
       ENDIF
       GO TO 600
      
       ENDIF

9999   CONTINUE
       CLOSE (10)
       CLOSE (20)
       IF(OUTFILE2(1:1).EQ.' ')CLOSE (30)
      
       END  

