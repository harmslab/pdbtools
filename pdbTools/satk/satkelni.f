CMJH  PROGRAM SATKEL(INFILE1,INFILE3,OUTFILE1,OUTFILE2,OUTFILE3,
CMJH *               OUTFILE4,OUTFILE5)
C**************************************************************** 
C...
CMJH    NOTES ON LINUX (x86) PORT (MJH, 5/19/2005)
CMJH    This code as originaally written on an HP 9000-825.  At some
CMJH    point prior to May, 2005 it was recompiled for SGI using
CMJH    the MIPSpro 7.2.1 F77 compiler.  This version has been
CMJH    slightly altered to allow it to compile in a Linux x86
CMJH    environment using gcc 3.4.3.  A comparison of the output from
CMJH    the SGI and linux binaries is in the test/ folder.  No
CMJH    original code has been deleted, but rather commented out with
CMJH    the 'CMJH' comment line.  Specifically, four READ statements
CMJH    had to be altered (see commented lines).  I also had to add
CMJH    three junk variables (MJHDEB(2), DEBRIS, and DEBRIS1).  
CMJH
CMJH    Michael Harms, May 19, 2005
C...                                                             
C...  Elect calculates the proton and ion binding behavior of a 
C...  protein. Fed with a table of work factors calculated with Tanford-
C...  Kirkwood theory, or with an otherwise suitable dielectric model,
C...  the coordinates of a macromolecule, data such as the pKint of
C...  the titratable sites in the macromolecule, their normalized
C...  accessibility to solvent, and the valence charge of the titra-
C...  table site, this program will calculate the extent of binding at
C...  any one site as a function of pH and ionic strength, the total
C...  work necessary to load, fully or partially all sites, the pK1/2
C...  of binding sites, their titration profiles, etc, etc, etc. 
C...
C...  A pKint value of 0 identifies that atom as titratable ion. A
C...  pKint value of -1 identifies that atom as a site that is always
C...  fully charged.
C...                                                               
C...     Unit Number        File Extension
C...         60                  .sin
C...         10                  .wij
C...         20                  .pot
C...         30                  .plt
C...         40                  .ijd
C...         50                  output
C...         70                  .pov
C...        
C...     This program is a version of the original ELECT program written
C...  by S. J. Shire.                       
C...                                                             
C...  Bertrand Garcia-Moreno E.
C...                                                       
C...  Fortran 77 on an HP 9000-835 SRX
C...                                                      
C***********************************************************************

C        Copyright 1988, Bertrand Garcia-Moreno E.
C        This program is distributed under General Public License v. 3.  See the
C        file COPYING for a copy of the license.  

      IMPLICIT NONE
      DIMENSION SRION(20),SSRION(20),X(400),Y(400),Z(400)
      COMMON  /A/ WIJ(20,800),RIJ(400,400),PKINT(400),CHARGE(400),
     *        BURY(400),RESNAM(400),RESNUM(400),M,MM,PHL,PHH,PHIN,
     *        IONST,TEMP,IJNUM(40),IJRSNAM(40),IJATNAM(40),DIJNUM,
     *        DIJLIM,ATNAM(400),GRPHNUM,GNAM(40),GNUM(40),ANSWER1,
     *        ANSWER2,ANSWER3,ANSWER4,IRESNAM(40),IRESNUM(40),IRESI,
     *        SCONC(20),SALTNAM(2),SALTCONC(2,20),ANFCTR(2),CAFCTR(2),
     *        CACHAR(2),ANCHAR(2),POTMDL,SATERM,SAMOD,NEWD,IONC,IS,
     *        GCHN(40),MJHDEB(2)
      CHARACTER HEADER*80,IJATNAM*5,IJRSNAM*5,GNAM*5,
     *        RESNAM*5,ATNAM*5,INFILE1*20,INFILE2*20,INFILE3*20,
     *        OUTFILE1*20,OUTFILE2*20,OUTFILE3*20,OUTFILE4*20,
     *        ANSWER1*1,ANSWER2*1,ANSWER3*1,ANSWER4*1,IRESNAM*5,
     *        SALTNAM*7, SATERM*1,GCHN*1,OUTFILE5*20,DEBRIS*80,
     *        DEBRIS1*80,MJHDEB*1
      REAL PHL,PHH,PHIN,SRION,DIJLIM,B,BB,X,Y,Z,PKINT,CHARGE,BURY,NEWD,
     *        TEMP,IONC,IONST,SCONC,SALTCONC,WIJ,SSRION,
     *        RIJ,ANSTRN,CASTRN
      INTEGER M,GRPHNUM,DIJNUM,N,IJNUM,GNUM,MM,RESNUM,IRESNUM,ANFCTR,
     *        CAFCTR,POTMDL,SAMOD,CACHAR,ANCHAR,SALTNUM,CONCNUM,IS,NN,J,
     *        IRESI,II,I
      DIJNUM=0
      GRPHNUM=0

C***********************************************************************
C...                                                                
C...  The header is read and written.
C...                                                              
C***********************************************************************

CMJH  CALL statements required because PROGRAM statement commented out.  

      CALL GETARG(1,INFILE1)
      CALL GETARG(2,INFILE3)
      CALL GETARG(3,OUTFILE1)
      CALL GETARG(4,OUTFILE2)
      CALL GETARG(5,OUTFILE3)
      CALL GETARG(6,OUTFILE4)
      CALL GETARG(7,OUTFILE5) 

      OPEN (UNIT=90,FILE=INFILE1,STATUS='OLD')
      READ(90,'(A1)') ANSWER1
        IF((ANSWER1.EQ.'y').OR.(ANSWER1.EQ.'Y'))THEN
        ANSWER1='Y'
        READ(90,'(A1)') ANSWER4
          IF((ANSWER4.EQ.'n').OR.(ANSWER4.EQ.'N'))THEN
            ANSWER4='N'
            READ(90,*) IRESI
              DO 95 I=1,IRESI
              READ(90,*) IRESNAM(I),IRESNUM(I)
95            CONTINUE
          ELSE
            ANSWER4='Y'
          ENDIF
        ENDIF
      READ(90,'(A1)') ANSWER2
      IF(ANSWER2.EQ.'y'.OR.ANSWER2.EQ.'Y') ANSWER2='Y'
      READ(90,'(A1)') ANSWER3
      IF(ANSWER3.EQ.'y'.OR.ANSWER3.EQ.'Y') ANSWER3='Y'
      READ(90,*) M, PHH, PHL, PHIN 
      NN = (PHH - PHL)/PHIN + 1
      IF(M.GT.400) STOP 'There are more than 400 charged atoms!!'
  
      READ(90,'(I1)') POTMDL
      IF (POTMDL.EQ.1) THEN
        READ(90,'(A20)') INFILE2
        OPEN(UNIT=10,FILE=INFILE2,STATUS='OLD')
        READ(90,'(I2)') CONCNUM
          DO 100 I=1,CONCNUM
          READ(90,*) SSRION(I)
          IF(I.GT.20) STOP 'Raise the upper index of SRION array!'
100       CONTINUE
      ELSE IF (POTMDL.EQ.2.OR.POTMDL.EQ.3) THEN
        IF(POTMDL.EQ.2) THEN
          READ(90,*) TEMP
        ELSE IF(POTMDL.EQ.3) THEN 
          READ(90,*) TEMP
          READ(90,*) NEWD
          READ(90,*) NEWD
        ENDIF
        READ(90,*) SALTNUM 
          DO 205 I = 1, 2
          READ(90,'(A7)') SALTNAM(I)
          READ(90,*) CACHAR(I)
          READ(90,*) ANCHAR(I)
205       CONTINUE  
        READ(90,*) CONCNUM 
          DO 201 I = 1, CONCNUM
          SRION(I) = 0.0 
          SCONC(I) = 0.0
            DO 210 J = 1, 2
            READ(90,*) SALTCONC(J,I)
            IF(ABS(ANCHAR(J)).EQ.ABS(CACHAR(J))) THEN
              ANFCTR(J) = 1
              CAFCTR(J) = 1
            ELSE IF(ABS(ANCHAR(J)).EQ.1.AND.ABS(CACHAR(J)).EQ.2) THEN
              ANFCTR(J) = 2
              CAFCTR(J) = 1
            ELSE IF(ABS(ANCHAR(J)).EQ.1.AND.ABS(CACHAR(J)).EQ.3) THEN
              ANFCTR(J) = 3
              CAFCTR(J) = 1
            ELSE IF(ABS(ANCHAR(J)).EQ.2.AND.ABS(CACHAR(J)).EQ.3) THEN
              ANFCTR(J) = 3
              CAFCTR(J) = 2
            ELSE
              PRINT *,"I do not know this type of salt!!!"
              STOP
            ENDIF
            SCONC(I) = SALTCONC(J,I) + SCONC(I)
              CASTRN = SALTCONC(J,I) * CAFCTR(J) * (CACHAR(J)**2)
              ANSTRN = SALTCONC(J,I) * ANFCTR(J) * (ANCHAR(J)**2)
            SRION(I) = CASTRN + ANSTRN + SRION(I)
210         CONTINUE
  
          SRION(I) = SRION(I) / 2
201       CONTINUE
      ENDIF  
     
      READ(90,'(A1)') SATERM
        IF (SATERM.EQ.'y'.OR.SATERM.EQ.'Y') THEN
          SATERM = 'Y'
          READ(90,'(I1)') SAMOD
        ENDIF

      READ(90,*) GRPHNUM
      IF(GRPHNUM.GT.40) STOP 'Raise index of GNAM and GNUM!'
        IF(GRPHNUM.NE.0) THEN
          DO 105 I=1,GRPHNUM
          READ(90,*) GNAM(I), GNUM(I), GCHN(I)
105       CONTINUE
        ENDIF
      READ(90,*) DIJNUM
      IF(DIJNUM.GT.40) STOP 'Raise index of IJRSNAM and IJNUM!!!'
        IF(DIJNUM.NE.0)THEN
          READ(90,*) DIJLIM
            DO 110 I=1,DIJNUM
            READ(90,*) IJRSNAM(I), IJNUM(I)
110         CONTINUE
        ENDIF

      OPEN(UNIT=60,FILE=INFILE3,STATUS='OLD')
      OPEN(UNIT=20,FILE=OUTFILE1,STATUS='NEW')
      OPEN(UNIT=30,FILE=OUTFILE2,STATUS='NEW')
      OPEN (UNIT=70,FILE=OUTFILE5,STATUS='NEW')
      IF(DIJNUM.NE.0)OPEN(UNIT=40,FILE=OUTFILE3,STATUS='NEW')
      OPEN(UNIT=50,FILE=OUTFILE4,STATUS='NEW')

      DO 500 I=1,5                                                   
      READ(60,'(A80)') HEADER                                         
      WRITE(20,'(A80)') HEADER                                        
      WRITE(30,'(A80)') HEADER   
      IF(DIJNUM.NE.0)WRITE(40,'(A80)') HEADER     
      WRITE(50,'(A80)') HEADER
500   CONTINUE                                                  

      DO 130 I=1,5
        IF (POTMDL.EQ.1) THEN
          READ(10,'(A80)') HEADER
            DO 140 II=20,50,10
            IF((II.EQ.40).AND.(DIJNUM.EQ.0)) GO TO 140
            WRITE(II,'(A80)') HEADER
140         CONTINUE
        ELSE IF (POTMDL.EQ.2) THEN
          DO 131 II=20,50,10
            IF(II.EQ.40.AND.DIJNUM.EQ.0) GO TO 131
            IF (I.EQ.1) THEN
              WRITE(II,'("Distance-Dependent delectric was used")')
            ELSE IF (I.EQ.2) THEN
              WRITE(II,'("  T = ",F7.2)') TEMP
            ELSE
              WRITE(II,'(" ")') 
            ENDIF 
131       CONTINUE 
        ELSE IF (POTMDL.EQ.3) THEN
          DO 132 II=20,50,10 
            IF(II.EQ.40.AND.DIJNUM.EQ.0) GO TO 132
            IF (I.EQ.1) THEN
              WRITE(II,'("Constant dielectric was used")')
            ELSE IF (I.EQ.2) THEN
              WRITE(II,'("  T = ",F7.2)') TEMP
            ELSE IF (I.EQ.3) THEN 
              WRITE(II,'("  D = ",F7.2)') NEWD 
            ELSE
              WRITE(II,'(" ")') 
            ENDIF
132       CONTINUE
        ENDIF
130   CONTINUE

      IF (POTMDL.EQ.1) THEN
CMJH        READ(10,'("TEMP = ",F6.2)') TEMP
        READ(10,'(A7,F6.2)') DEBRIS, TEMP
        READ(10,'(A80)') HEADER
      ENDIF

      DO 142 II=20,50,10
      IF((II.EQ.40).AND.(DIJNUM.EQ.0)) GO TO 142
      IF (POTMDL.EQ.1) THEN
        WRITE(II,'("TEMP = ",F6.2)') TEMP
        WRITE(II,'(A80)') HEADER
      ELSE
        WRITE(II,'(" ")')
        WRITE(II,'(" ")')
      ENDIF
142   CONTINUE

      IF (POTMDL.EQ.1) THEN
        READ(10,'(A3,F5.2,A10)') DEBRIS, B, DEBRIS1
CMJH        READ(10,'("B= ",F5.2," angstroms")') B
        READ(10,'(A80)') HEADER
      ENDIF

      DO 141 II=20,50,10
      IF((II.EQ.40).AND.(DIJNUM.EQ.0)) GO TO 141
      IF (POTMDL.EQ.1) THEN 
        WRITE(II,'("B= ",F5.2," angstroms")') B
        WRITE(II,'(A80)') HEADER
      ELSE
        WRITE(II,'(" ")') 
        WRITE(II,'(" ")') 
      ENDIF
141   CONTINUE
 
      WRITE(20,'(A20)') INFILE3
      WRITE(50,'(///"IONIC. STRN.  IONIC CONC. SALT 1    CONC.    TYPE "
     *          "SALT 2    CONC.   TYPE  ")')
                  
      IF (POTMDL.EQ.1) THEN
        READ  (10,'(A35,I3,A21)'), DEBRIS, N, DEBRIS1
CMJH        READ  (10,'("The Wij factors were calculated at ",I3," ionic",
CMJH     *               " strengths (M):")') N
        DO 120 I = 1, N
        READ (10,'(F12.8,1X,F12.8,2(1X,A7,1X,F10.7,1X,I1,A1,I1))')
     *        SRION(I), SCONC(I), SALTNAM(1), SALTCONC(1,I), ANFCTR(1),
     *        MJHDEB(1),CAFCTR(1), SALTNAM(2), SALTCONC(2,I), ANFCTR(2), 
     *        MJHDEB(2), CAFCTR(2)

CMJH        READ (10,'(F12.8,1X,F12.8,2(1X,A7,1X,F10.7,1X,I1,":",I1))')
CMJH     *        SRION(I), SCONC(I), SALTNAM(1), SALTCONC(1,I),  ANFCTR(1),
CMJH     *        CAFCTR(1), SALTNAM(2), SALTCONC(2,I), ANFCTR(2), CAFCTR(2)

          DO 121 II = 1, CONCNUM
          IF (ABS(SRION(I)-SSRION(II)).LT..000000001) THEN
            WRITE (50,'(F12.8,1X,F12.8,2(1X,A7,1X,F10.7,1X,I1,":",I1))')
     *            SRION(I),SCONC(I),SALTNAM(1),SALTCONC(1,I), ANFCTR(1),
     *            CAFCTR(1),SALTNAM(2),SALTCONC(2,I),ANFCTR(2),CAFCTR(2)
          ENDIF
121       CONTINUE
120     CONTINUE
      ELSE
        DO 122 I = 1, CONCNUM
          WRITE (50,'(F12.8,1X,F12.8,2(1X,A7,1X,F10.7,1X,I1,":",I1))')
     *          SRION(I),SCONC(I),SALTNAM(1),SALTCONC(1,I), ANFCTR(1),
     *          CAFCTR(1),SALTNAM(2),SALTCONC(2,I),ANFCTR(2),CAFCTR(2)
122     CONTINUE
      ENDIF

      WRITE(30,'("Number of concentrations studied: ",I3)') CONCNUM
      WRITE(30,'("Number of residues included in this file: ",I2)')
     *   GRPHNUM
      WRITE(30,'("The number of pH values is: ",I3,2X,F5.2,2X,F5.2)')
     *      NN, PHL, PHIN

      IF(GRPHNUM.EQ.0) THEN   
        WRITE(30,'(////////////)')
      ELSE
        WRITE(30,9010)
      ENDIF 

      IF(DIJNUM.NE.0)THEN
        WRITE(40,'("Number of residues included in this file: ",I2,//)')
     *       DIJNUM
        WRITE(40,9010)
      ENDIF

      WRITE(50,9010)
9010  FORMAT(///////1X,'NUM',2X,'RESIDUE',1X,'ATOM',3X,'CHARGE',4X,
     *'BURY',5X,'PKINT',7X,'X',6X,'Y',6X,'Z',/,30X,'(1-SA)',/)

      DO 510 I=1,M                                                 
      READ(60,'(1X,I4,2X,A5,1X,A5,2X,3(1X,F8.3),2X,F4.1,1X,F5.2,1X,
     *F4.2)') RESNUM(I), RESNAM(I), ATNAM(I), X(I), Y(I),
     *Z(I),CHARGE(I), PKINT(I), BURY(I)
      WRITE(50,9020) RESNUM(I), RESNAM(I), ATNAM(I), CHARGE(I), BURY(I),
     *PKINT(I),X(I),Y(I),Z(I)                           
        DO 550 II=1,DIJNUM
          IF((IJNUM(II).EQ.RESNUM(I)).AND.(RESNAM(I)(1:3).EQ.
     *      IJRSNAM(II))) THEN
            WRITE(40,9020)RESNUM(I), RESNAM(I), ATNAM(I), CHARGE(I),
     *           BURY(I), PKINT(I), X(I), Y(I), Z(I)
            IJRSNAM(II)=RESNAM(I)
          ENDIF
550     CONTINUE
        DO 560 II=1,GRPHNUM
          IF((GNUM(II).EQ.RESNUM(I)).AND.(RESNAM(I)(1:3).EQ.GNAM(II))
     *  .AND.((RESNAM(I)(5:5).EQ.GCHN(II)).OR.(GCHN(II).EQ.'Z'))) THEN
            WRITE(30,9021)RESNUM(I), RESNAM(I), ATNAM(I), CHARGE(I), 
     *           BURY(I), PKINT(I), X(I), Y(I), Z(I)
          ENDIF
560     CONTINUE
9021  FORMAT(1X,I3,2X,A5,3X,A5,3X,F4.1,5X,F4.2,5X,F5.2,1X,2X,
     *3(2X,F5.1),5X,"+ + +"/)
9020  FORMAT(1X,I3,2X,A5,3X,A5,3X,F4.1,5X,F4.2,5X,F5.2,1X,2X,
     *3(2X,F5.1),/)
510   CONTINUE                                                 
      WRITE(50,'(///)')
      WRITE(20,'("Low pH is ",F6.2,"   high pH is",F6.2,"   pH interval
     *is ",F6.2," There are ",I4," pH values")') PHL, PHH, PHIN, NN
      WRITE(20,'(I4," ionic strengths explored are: ")')CONCNUM
        DO 511 I = 1, CONCNUM
        WRITE(20,'(F12.8)') SSRION(I)
511     CONTINUE
      WRITE(20,'("There are ",I3," charges in this molecule")') M
      WRITE(20,'(3(3F8.3,F4.2))') (X(I), Y(I), Z(I), BURY(I), I=1,M)

C**********************************************************************
C...
C...   The following block fills in the Rij array.
C...
C**********************************************************************

      DO 520 I=1,M                                                  
        DO 530 J=1,M                                               
        RIJ(I,J)=SQRT((X(I)-X(J))**2+(Y(I)-Y(J))**2+(Z(I)-Z(J))**2) 
530     CONTINUE                                                   
520   CONTINUE                                                    

C***********************************************************************
C...                                                                   
C...  The following block chooses one ionic strength, calls the sub-
C...  routine that performs the actual calculations, and repeats this
C...  process as it goes through the DO loop.
C...                                                               
C***********************************************************************

      BB=B/0.1
      MM=INT(ANINT(BB))*2
     
      IF (POTMDL.EQ.1) THEN
        DO 200 I = 1, N 
        READ(10,'(6E12.5)') (WIJ(I,J),J=1,MM)
200     CONTINUE 
          DO 540 II=1,CONCNUM
            DO 541 I = 1, N
              IF(ABS(SRION(I)-SSRION(II)).LT..00000000001) THEN
                IONST = SRION(I)
                IONC  = SCONC(I) 
                   IS = I 
                CALL TITRATE
                GO TO 540
              ENDIF
541         CONTINUE   
540       CONTINUE
      ELSE 
        DO 199 I = 1, CONCNUM       
        IONST = SRION(I)
         IONC = SCONC(I)
           IS = I
        CALL TITRATE 
199     CONTINUE
      ENDIF
      END


      SUBROUTINE TITRATE                                   
C***********************************************************************
C...                                                                   
C...  Subroutine titrate is the subroutine that performs the calculation
C...  The DO 510 loop covers the main body of the subroutine. It 
C...  performs the iteration. A single pH value is dealt with every
C...  time the program flows past the DO 510 loop. Then delta pK
C...  factors are calculated for each one of the charges, the intrinsic
C...  pKs are corrected, and a new set of charges is calculated. If the
C...  difference between any consecutive set of charges is less than
C...  .01, the iteration is ended and the program flow proceeds on to a
C...  new pH value.
C...                                                                
C...  Some of the parameters of use to this subroutine are: 
C...                                                               
C...  CHARGEI()- The charge of group I at a given pH value
C...  MCHSPH   - The charge of the macromolecule at a specific pH value.
C...  TCHISPH  - Any one element in the MCHSPH array. 
C...                                                           
C...  ITNUMT   - Iteration number total at any given pH.
C...  ITNUM    - Current iteration number. 
C...                                                        
C...  DGISPH   - Delta G for group i at a specific pH. 
C...  DGISPHT  - Temporary DELTA G for group I at a specific pH.
C...  MDGSPH   - Protein DELTA G at a specific pH. 
C...  SDGISPH  - Any element of the MDGSPH array. 
C...                                                           
C...  This subruotine is a version by Bertrand Garcia-Moreno E. of  
C...  the subruotine written by S. J. Shire as modified by Jim Matthew.
C...                                                             
C...  Bertrand Garcia-Moreno E.
C...                                                         
C***********************************************************************

      IMPLICIT NONE
      DIMENSION CHARGEI(400,2), P(400,2), DELTAPK(400), SPH(100),     
     *          DGISPH(400),ITNUMT(100),MCHSPH(100),MDGSPH(100), 
     *          CHARGEZ(400),TCHARGE(400,2),TPH(400,2),PK12(400),    
     *          TDGIJ(400),DGIJ(40,400),DELTAIJ(40,400),GDELTPK(40,100),
     *          GP(40,100),GCHARGI(40,100),GDGISPH(40,100),
     *          GASSCON(40,100),GLOGK(40,100),MIONPOS(100),MIONNEG(100),
     *          MIONTOT(100),MCHRTOT(100)  
      COMMON  /A/ WIJ(20,800),RIJ(400,400),PKINT(400),CHARGE(400),
     *          BURY(400),RESNAM(400),RESNUM(400),M,MM,PHL,PHH,PHIN,
     *          IONST,TEMP,IJNUM(40),IJRSNAM(40),IJATNAM(40),DIJNUM,
     *          DIJLIM,ATNAM(400),GRPHNUM,GNAM(40),GNUM(40),ANSWER1,
     *          ANSWER2,ANSWER3,ANSWER4,IRESNAM(40),IRESNUM(40),IRESI,
     *          SCONC(20),SALTNAM(2),SALTCONC(2,20),ANFCTR(2),CAFCTR(2),
     *          CACHAR(2),ANCHAR(2),POTMDL,SATERM,SAMOD,NEWD,IONC,IS,
     *          GCHN(40),MJHDEB(2)
      REAL LOGK,IONPOS,IONNEG,IONTOT,PHL,PHH,PHIN,PKINT,CHARGE,BURY,
     *          KINTAN(2),KINTCA(2),MIONPOS,MIONNEG,MIONTOT,MCHRTOT,
     *          MCHSPH,MDGSPH,NEWD,TEMP,IONC,IONST,SCONC,CONST2,CACONC,
     *          ANCONC,CACONC1,CACONC2,ALPHA,ANCONC1,ANCONC2,ASSCONS,
     *          AVGBURY,BJ,CHARGEI,CHARGEZ,CHION,CHIONTM,CHTEMP,CONST,
     *          CONST1,CVP,CVP1,CVP2,DELTA,DELTAIJ,DELTAPK,DGIJ,
     *          DGISPH,DGISPHT,DIJLIM,FRACT,GASSCON,GCHARGI,GDELTPK,
     *          GDGISPH,GLOGK,GP,KAN1,KAN2,KAPPA,KCA1,KCA2,KINTAN1,
     *          KINTAN2,KINTCA1,KINTCA2,P,PH,PHI,PK12,RIJ,RIJE,SALTCONC,
     *          SDGISPH,SLOPE,SPH,TCHARGE,TCHISPH,TDGIJ,TDGISPH,TPH,
     *          TRUNC,WIJ,WIJT,Z1,Z2,SUMPK1,SUMPK2,
     *          CONST3,CONST4
      CHARACTER RESNAM*5,ATNAM*5,IJRSNAM*5,IJATNAM*5,GNAM*5,ANSWER1*1,
     *          ANSWER2*1,ANSWER3*1,ANSWER4*1,IRESNAM*5,SALTNAM*7,
     *          SATERM*1,GCHN*1,MJHDEB*1
      INTEGER M, RESNUM, DIJNUM, IJFLAG, FLAGIJ, GRPHNUM, GNUM, IRESNUM,
     *          NN, SAMOD, POTMDL, CACHAR, ANCHAR, IS, ANFCTR,
     *          CAFCTR,ANFCTR1,ANFCTR2,CAFCTR1,CAFCTR2,I,II,III,IJNUM,
     *          IONFLAG,IONFLG1,IRES,IRESI,ITNUM,ITNUMT,J,JJ,K,KK,LL,MM,
     *          PKFLAG

C***********************************************************************
C...                                                             
C...   The number of pH values is calculated, their values are ini-
C...   tialized.
C...                                                              
C***********************************************************************

      PHI=PHH-PHL                                                   
      NN=PHI/PHIN +1                                               
      BJ=0.0                                                      

      DO 500 J=1,NN                                            
      SPH(J)=PHL+(BJ*PHIN)                                    
500   BJ=BJ+1.0                                              

C***********************************************************************
C...                                                               
C...  CONST   converts Wij factors from ergs to kcal/mol.
C...  CONST1  is the 2.303*RT factor. 
C...  CONST2  is the RT factor
C...  TRUNC   describes truncation criterion for proton binding.
C...          Ideally  it should be set to 0.001 but NA and Cl binding
C...          dont converge with this stringent criterion. .01 works.
C...
C...  KINTCA and KINTAN are intrinsic binding constants for cations and
C...          anions (ie binding when Vel=0).
C...                                                            
C...  The values of the concentration of cations and anions is also
C...  calculated at this point. Notice that presently, the program is
C...  only set up to handle the competition between binding of
C...  monovalent and divalent cations. Competition between binding of
C...  cations of any other valency has to be written in explicitly.
C...
C***********************************************************************

      CONST   = ((.314549158)*1.E16)/TEMP
      CONST1  = TEMP*0.0045761                                       
      CONST2  = TEMP*.001987
      CONST3  = 1.43942*1.E13
      CONST4  = 2.30683*1.E-11
      TRUNC   = .01
      TRUNC   = .001
      KINTCA(1) = 1.0
      KINTCA(2) = 1.0
      KINTAN(1)  = 1.0
      KINTAN(2)  = 1.0

      DO 499 J = 1 ,M                                          
             PK12(J) = 0.0                                        
            TPH(J,2) = 0.0
            TPH(J,1) = 0.0
        TCHARGE(J,2) = 0.0
        TCHARGE(J,1) = 0.0
499   CONTINUE

      ANCONC = SALTCONC(1,IS)*ANFCTR(1) + SALTCONC(2,IS)*ANFCTR(2)
      CACONC = SALTCONC(1,IS)*CAFCTR(1) + SALTCONC(2,IS)*CAFCTR(2)

C***********************************************************************
C...                                                               
C...  The DO 510 loop covers the body of the program. It chooses the
C...  pH at which the iteration is performed.
C...  Notice that the values of the variables KINTCA and KINTAN
C...  (Kint for cation and anion binding) have been set equal to 1.
C...  There is room for improvement in this approximation.
C...                                                               
C***********************************************************************

      DO 510 LL=1,NN                                               

        CHARGEZ(LL)=0.0
        DO 520 J=1,M                                              
        TDGIJ(J)=0.0
        DGISPH(J)=0.0
        CHARGEI(J,1)=0.0
        CHARGEI(J,2)=0.0
        DELTAPK(J)=0.0
520     CONTINUE

        IONPOS=0.0                                
        IONNEG=0.0
        IONTOT=0.0
        PH=SPH(LL)                                               

        DO 530 I=1,M                                          
        P(I,1)=PKINT(I)                                        
530     CONTINUE

C***********************************************************************
C...                                                                
C...  The first part of the DO 540 loop that starts below is used to
C...  calculated the charges of every titratable side chain at any
C...  given pH, using the Henderson-Hasselbalch equation. Note that 
C...  single atoms other than site-bound ions that bear a charge larger
C...  than + or - 1 are not allowed to titrate, but rather are assumed
C...  to be fully charged through the pH range investigated.
C...                                                                 
C***********************************************************************

      DO 540 K=1,200                                              
C      write(50,'(/"IT    M     CHARGE       P(I,1)       P(I,2) ",
C     *     "     DPK")')
        IF(K.EQ.1) THEN                                            
          KK=1000                                                 
          IONFLG1=0                                              
        ENDIF                                                   
        SUMPK1=0.0 
        SUMPK2=0.0 
        CHION=0.0                                              
        TCHISPH=0.0                                           
        ITNUM=K                                              
        SDGISPH=0.0                                         

        DO 550 J=1,M                                     

          IF(PKINT(J).NE.0.0.AND.PKINT(J).NE.-1) GO TO 22  
          IF((PKINT(J)+1.0).LT.0001) THEN
            CHARGEI(J,1)=CHARGE(J)
            GO TO 550
          ENDIF 

          IF (CHARGE(J).LT.0.01) THEN
            IF(SALTNAM(2)(1:4).EQ.'NONE'.OR.(ANFCTR(1).EQ.ANFCTR(2).AND.
     *         CAFCTR(1).EQ.CAFCTR(2))) THEN 
              CVP = -1.0*TDGIJ(J)/CONST2                          
              KAN1 = KINTAN(1)*(2.7183**CVP)  
C              IF(CVP.LT.0.0) KAN1 = 0.0                            
              CHARGEI(J,1) = (KAN1*ANCONC)/(1.0+(KAN1*ANCONC))     
              CHARGEI(J,1) = CHARGE(J)*CHARGEI(J,1)                 
            ELSE
              IF (CHARGE(J).EQ.CAFCTR(1)) THEN
                ANFCTR1=ANFCTR(1)
                ANFCTR2=ANFCTR(2)
                CAFCTR1=CAFCTR(1) 
                CAFCTR2=CAFCTR(2) 
                KINTAN1=KINTAN(1)
                KINTAN2=KINTAN(2)
                ANCONC1=ANFCTR1*SALTCONC(1,IS)
                ANCONC2=ANFCTR2*SALTCONC(2,IS)
              ELSE IF(CHARGE(J).EQ.CAFCTR(2)) THEN 
                ANFCTR1=ANFCTR(2)
                ANFCTR2=ANFCTR(1)
                CAFCTR1=CAFCTR(2)
                CAFCTR2=CAFCTR(1)
                KINTAN1=KINTAN(2)
                KINTAN2=KINTAN(1)
                ANCONC1=ANFCTR2*SALTCONC(2,IS)
                ANCONC2=ANFCTR1*SALTCONC(1,IS)
              ELSE
                PRINT *,"The charge of the ions in .sin and .wij dont ",
     *                  "match."
                STOP
              ENDIF
              CVP1 = -1.0 * TDGIJ(J)/CONST2 
              CVP2 = -1.0 * TDGIJ(J)*CAFCTR2/CONST2*CAFCTR1
              KAN1 = KINTAN1*(2.7183**CVP1)
              KAN2 = KINTAN2*(2.7183**CVP2)
              Z1 = (KAN1*ANCONC1)/(1.0 + KAN2*ANCONC2 + KAN1*ANCONC1) 
              Z2 = (KAN2*ANCONC2)/(1.0 + KAN1*ANCONC1 + KAN2*ANCONC2) 
              Z1 = Z1*ANFCTR1
              Z2 = Z2*ANFCTR2
              CHARGEI(J,1) = Z1 + Z2
            ENDIF
          ELSE
            IF(SALTNAM(2)(1:4).EQ.'NONE'.OR.(CAFCTR(1).EQ.CAFCTR(2).AND.
     *        ANFCTR(1).EQ.ANFCTR(2))) THEN 
              CVP = -1.0*TDGIJ(J)/CONST2                          
              KCA1 = KINTCA(1)*(2.7183**CVP)  
C              IF(CVP.LT.0.0) KCA1 = 0.0                            
              CHARGEI(J,1) = (KCA1*CACONC)/(1.0+(KCA1*CACONC))
              CHARGEI(J,1) = CHARGE(J)*CHARGEI(J,1)                 
            ELSE
              IF (CHARGE(J).EQ.ANFCTR(1)) THEN
                CAFCTR1=CAFCTR(1)
                CAFCTR2=CAFCTR(2)
                ANFCTR1=ANFCTR(1)
                ANFCTR2=ANFCTR(2)
                KINTCA1=KINTCA(1)
                KINTCA2=KINTCA(2)
                CACONC1=CAFCTR1*SALTCONC(1,IS)
                CACONC2=CAFCTR2*SALTCONC(2,IS)
              ELSE IF(CHARGE(J).EQ.ANFCTR(2)) THEN 
                CAFCTR1=CAFCTR(2)
                CAFCTR2=CAFCTR(1)
                ANFCTR1=ANFCTR(2)
                ANFCTR2=ANFCTR(1)
                KINTCA1=KINTCA(2)
                KINTCA2=KINTCA(1)
                CACONC1=CAFCTR2*SALTCONC(2,IS)
                CACONC2=CAFCTR1*SALTCONC(1,IS)
              ELSE
                PRINT *,"The charge of the ions in .sin and .wij dont ",
     *                  "match."
                STOP
              ENDIF
              CVP1 = -1.0 * TDGIJ(J)/CONST2 
              CVP2 = -1.0 * TDGIJ(J)*ANFCTR2/CONST2*ANFCTR1
              KCA1 = KINTCA1*(2.7183**CVP1)
              KCA2 = KINTCA2*(2.7183**CVP2)
              Z1 = (KCA1*CACONC1)/(1.0 + KCA2*CACONC2 + KCA1*CACONC1) 
              Z2 = (KCA2*CACONC2)/(1.0 + KCA1*CACONC1 + KCA2*CACONC2) 
              Z1 = Z1*CAFCTR1
              Z2 = Z2*CAFCTR2
              CHARGEI(J,1) = Z1 + Z2
            ENDIF
          ENDIF
          GO TO 550

22        FRACT=(10.**(PH-P(J,1)))                              
          ALPHA=FRACT/(1.+FRACT)                             
          IF (CHARGE(J)) 10,20,30                         

30        CHARGEI(J,1)=(1.-ALPHA)*CHARGE(J)                
          GO TO 550                                     

20        CHARGEI(J,1)=0.                                
          GO TO 550                                   

10        CHARGEI(J,1)=ALPHA* CHARGE(J)                
          GO TO 550                                 

550     CONTINUE                                   

C***********************************************************************
C...                                                                   
C...     The following block calculates effective pK values for each of
C...  the amino acids. The newly defined pK values are then used to cal-
C...  culate a new set of charges, which are then used to calculate a 
C...  new  pK. An iteration ensues.
C...                                                           
C***********************************************************************

          J=1
          FLAGIJ=0
          IJFLAG=0

          DO 560 I=1,M                                                 
            DGISPH(I)=0.                                                
            DGISPHT=0.00                                               
            TDGIJ(I)=00.00                                    
            TDGISPH=00.00
C            IF (CHARGEI(I,1).EQ.0.)GO TO 560                          
            DELTA=0.0                                               
              DO 561 III=1,DIJNUM
                IF((IJNUM(III).EQ.RESNUM(I)).AND.(IJRSNAM(III).EQ.
     *              RESNAM(I)))THEN
                  FLAGIJ=1
                  IJFLAG=IJFLAG + 1
                  GO TO 562
                ENDIF
561           CONTINUE
562         CONTINUE
                 
            DO 570 J=1,M                                         
            DGISPHT=0.0                                         
            TDGISPH=0.0                                        
            DGIJ(III,J)=0.0

            IF(ABS(I-J).LT.00001) GO TO 50
  40          RIJE=RIJ(I,J)                                  
                IF (SATERM.EQ.'Y') THEN
                  IF (SAMOD.EQ.1) THEN
                    AVGBURY=(BURY(I)+BURY(J))/2                 
                  ELSE IF (SAMOD.EQ.2) THEN 
                    AVGBURY=(BURY(I)*BURY(J))
                  ELSE IF (SAMOD.EQ.3) THEN
                    AVGBURY=BURY(J)
                  ENDIF
                ELSE
                  AVGBURY = 1.0
                ENDIF

                IF (POTMDL.EQ.1) THEN
                  II=RIJE*10.                                   
                  IF (II.GT.MM) II=MM                          
                  WIJT=WIJ(IS,II)
                ELSE
                  IF (POTMDL.EQ.2) NEWD=RIJE
                  KAPPA = SQRT(1/(NEWD*TEMP))*50.2912*SQRT(IONST)
                  WIJT = (332.05/(RIJE*NEWD))*(1/2.7183**(KAPPA*RIJE))
                  WIJT = WIJT/CONST3
                ENDIF

             DELTA=WIJT*CHARGEI(J,1)*AVGBURY*(-1.) + DELTA  
             DGISPHT=WIJT*CHARGEI(J,1)*AVGBURY*CHARGEI(I,1)*CONST*CONST1
              IF(PKINT(I).EQ.0.0) THEN                               
              TDGISPH=WIJT*CHARGEI(J,1)*AVGBURY*CHARGE(I)*CONST*CONST1
              DGISPHT=WIJT*CHARGEI(J,1)*AVGBURY*CHARGEI(I,1)
     *                                 *CONST*CONST1
C              write (50,'("TDGISPH =",F9.3)') TDGISPH
              ENDIF                                                    
              GO TO 60                                                
 50           CONTINUE                                               
              TDGISPH=0.0                                           
              DGISPHT=0.0                                          
  60          DGISPH(I)=DGISPH(I)+DGISPHT                         
              TDGIJ(I)=TDGIJ(I)+TDGISPH                          
            
              IF(FLAGIJ.EQ.1) THEN
                DGIJ(IJFLAG,J)=DGISPHT
                DELTAIJ(IJFLAG,J)=CONST*WIJT*CHARGEI(J,1)*AVGBURY*(-1.)
                IF(J.EQ.I) DELTAIJ(IJFLAG,J)=0.0
              ENDIF 

  570       CONTINUE                                            
          DELTAPK(I)=DELTA*CONST                               
          P(I,2)=PKINT(I)+DELTAPK(I)                            
          P(I,1)=(P(I,1)+P(I,2))/2.0
C       write(50,'(I3,2X,I3,2X,F9.3,2X,F9.3,2X,F9.3,2X,F9.3,2X,F9.3)')
C     *      K,I,CHARGEI(I,1),P(I,1),P(I,2),DELTAPK(I),DGISPHT
          IF(K.EQ.1) P(I,1)=P(I,2)
          FLAGIJ=0
 560      CONTINUE                                           

C***********************************************************************
C...                                                                   
C...   This loop calculates the pH- and ionic strength -dependent charge
C...   of the entire macro- molecule due to protons (TCHISPH) and bound
C...   ions (CHION), and establishes whether the iteration is to be
C...   continued or not.
C...                                                                
C***********************************************************************

          DO 580 J=1,M                                               
            IF(PKINT(J).EQ.0.0) THEN                              
              CHION=CHARGEI(J,1) + CHION
              IONFLAG=1                                      
            ELSE                                          
              TCHISPH=CHARGEI(J,1)+TCHISPH                 
            ENDIF                                        
            SUMPK1=SUMPK1+P(I,1)
            SUMPK2=SUMPK2+P(I,2)
          SDGISPH=SDGISPH+DGISPH(J)                    
 580      CONTINUE                                       
          SDGISPH = SDGISPH*0.5
        IF (K.EQ.1) GO TO 539                           
        DO 581 J=1,M
            IF(ABS(ABS(P(J,1))-ABS(P(J,2))).GT.TRUNC) GO TO 539
581     CONTINUE
        DO 582 J=1,M
          IF(ABS(ABS(CHARGEI(J,1))-ABS(CHARGEI(J,2))).GT.TRUNC) 
     *         GO TO 539
582     CONTINUE
        GO TO 70   
539     CHIONTM=CHION                                       
        PKFLAG=0
        DO 541 J=1,M
          CHARGEI(J,2)=CHARGEI(J,1)
541     CONTINUE 

540     CHTEMP=TCHISPH                                     

C***********************************************************************
C...                                                                   
C...   The algorithm flows into the 70 CONTINUE line when the iteration
C...   has been completed.
C...
C...   The pK1/2 of every single group is also calculated at this stage.
C...                                                                   
C***********************************************************************

 70   CONTINUE                                                         
      DO 610 J=1,M                                                     
        IF((ABS(CHARGE(J)).NE.1).OR.(PKINT(J).EQ.0.0)) THEN           
          PK12(J)=-20.0                                              
        ELSE IF(CHARGE(J).EQ.-1.) THEN                              
          IF(TPH(J,2).NE.0.0)GO TO 80                              
          IF(ABS(CHARGEI(J,1)).LT..5) THEN                          
            TCHARGE(J,1)=CHARGEI(J,1)                              
            TPH(J,1)=PH                                         
            GO TO 80                                           
          ELSE IF(ABS(CHARGEI(J,1)).GT..5) THEN                 
            TCHARGE(J,2)=CHARGEI(J,1)                          
            TPH(J,2)=PH                                     
            SLOPE=(TCHARGE(J,2)-TCHARGE(J,1))/(TPH(J,2)-TPH(J,1)) 
              PK12(J)=TPH(J,1)+(-.5-TCHARGE(J,1))/SLOPE          
            IF(TPH(J,1).EQ.0.0) PK12(J)=0.0                     
          ENDIF                                                
80    CONTINUE                                                
        ELSE IF(CHARGE(J).EQ.1) THEN                         
          IF(TPH(J,2).NE.0.0)GO TO 90                       
          IF(ABS(CHARGEI(J,1)).GT..5) THEN                   
            TCHARGE(J,1)=CHARGEI(J,1)                       
            TPH(J,1)=PH                                  
            GO TO 90                                    
          ELSE IF(ABS(CHARGEI(J,1)).LT..5) THEN          
            TCHARGE(J,2)=CHARGEI(J,1)                   
            TPH(J,2)=PH                              
            SLOPE=(TCHARGE(J,2)-TCHARGE(J,1))/(TPH(J,2)-TPH(J,1)) 
            PK12(J)=TPH(J,1)+(.5-TCHARGE(J,1))/SLOPE             
            IF(TPH(J,1).EQ.0.0) PK12(J)=0.0
          ENDIF                                                 
 90     CONTINUE                                               
        ENDIF                                               
610    CONTINUE                                              

C**********************************************************************
C...
C...  Results are printed
C...
C**********************************************************************
      IF (ANSWER1.EQ.'Y') THEN
        WRITE(50,'(/"IONIC STRN. =",F12.8,"  IONIC CONC. =",F12.8,
     *        "     PH = ",F6.2, "     # OF ITERATIONS =",I3)')IONST,
     *        IONC, PH, ITNUM
        WRITE(50,'("GRP  ",3X,"RES",3X,"PKINT",5X,"DPK",5X,"PK",5X,  
     *       "ZAA",5X,"DELTA G",4X,"G",6X,"Ka",4X,"pKion")')
      ENDIF
        IF(DIJNUM.NE.0) WRITE(40,'(////////5X,"IONIC STRN. = ",F12.8,
     *  "  IONIC CONC. =",F12.8,"     PH = ",F6.2,/)') IONST, IONC, PH
        DO 590 I=1,M                                     
          IF(PKINT(I).EQ.0.0) THEN                        
            ASSCONS=2.7183**(-1*TDGIJ(I)/(.001987*TEMP))             
            LOGK=ALOG(ASSCONS)                            
              IF(ANSWER1.EQ.'Y') THEN
                IF(ANSWER4.EQ.'Y') THEN
                 WRITE(50,'(A5,I5,4X,F5.2,3(3X,F5.2),3X,F6.2,1X,F6.2,1X,
     *                F8.2,1X,F6.2)') RESNAM(I),RESNUM(I),PKINT(I),
     *                DELTAPK(I),P(I,1),CHARGEI(I,1), DGISPH(I),
     *                TDGIJ(I),ASSCONS,LOGK
                ELSE IF (ANSWER4.EQ.'N') THEN
                  DO  95 IRES=1,IRESI
                    IF((RESNAM(I).EQ.IRESNAM(IRES)).AND.(RESNUM(I).EQ.
     *                IRESNUM(IRES))) THEN
                      WRITE(50,'(A5,I5,4X,F5.2,3(3X,F5.2),3X,F6.2,1X,
     *                     F6.2,1X,F8.2,1X,F6.2)') RESNAM(I),RESNUM(I),
     *                     PKINT(I),DELTAPK(I),P(I,1),CHARGEI(I,1), 
     *                     DGISPH(I),TDGIJ(I),ASSCONS,LOGK
                    ENDIF
95                CONTINUE
                ENDIF           
              ENDIF
            IF(CHARGE(I).GT.0.0) THEN                               
              IONPOS=IONPOS+CHARGEI(I,1)                             
            ELSE IF(CHARGE(I).LT.0.0) THEN                        
              IONNEG=IONNEG+CHARGEI(I,1)                           
            ENDIF                                               
          ELSE                                               
            IF(ANSWER1.EQ.'Y') THEN
              IF(ANSWER4.EQ.'Y') THEN
                WRITE(50,'(A5,I5,4X,F5.2,3(3X,F5.2),3X,F6.2,2X,F6.2)')
     *          RESNAM(I),RESNUM(I), PKINT(I), DELTAPK(I), P(I,1),
     *          CHARGEI(I,1), DGISPH(I),ABS(CHARGEI(I,1))*(1-BURY(I))
              ELSE IF (ANSWER4.EQ.'N') THEN
                DO  96 IRES=1,IRESI
                  IF((RESNAM(I).EQ.IRESNAM(IRES)).AND.(RESNUM(I).EQ.
     *              IRESNUM(IRES))) THEN
                    WRITE(50,'(A5,I5,4X,F5.2,3(3X,F5.2),3X,F6.2,F6.2)')
     *              RESNAM(I),RESNUM(I),PKINT(I), DELTAPK(I),P(I,1),
     *              CHARGEI(I,1),DGISPH(I),PK12(I),
     *              ABS(CHARGEI(I,1))*(1-BURY(I))
                  ENDIF
96              CONTINUE
              ENDIF           
            ENDIF 
          ENDIF
         DO 592 III=1,GRPHNUM
          IF((GNUM(III).EQ.RESNUM(I)).AND.(RESNAM(I)(1:3).EQ.GNAM(III))
     *  .AND.((RESNAM(I)(5:5).EQ.GCHN(III)).OR.(GCHN(III).EQ.'Z'))) THEN
           GDELTPK(III,LL)=DELTAPK(I)
           GP(III,LL)=P(I,1)
           GCHARGI(III,LL)=CHARGEI(I,1)
           GDGISPH(III,LL)=DGISPH(I)
             IF(PKINT(I).EQ.0.0) THEN
               GASSCON(III,LL)=ASSCONS
               GLOGK(III,LL)=LOGK
             ENDIF
           ENDIF
592      CONTINUE
         DO 591 III=1,DIJNUM
           IF((RESNUM(I).EQ.IJNUM(III)).AND.(RESNAM(I).EQ.IJRSNAM(III)))
     *        THEN
              WRITE(40,'(//,1X,"GROUP",1X,"ATOM ","RES# ",1X,"GROUP",1X,
     *             "ATOM ","RES# ",3X," Zi  ",2X," Zj  "," Rij ",4X,
     *             " DGij",3X,"DpKij")')
              DO 599 J=1,M
                IF(ABS(DIJLIM).LT.ABS(DGIJ(III,J)))THEN 
                  WRITE(40,'(1X,A5,A5,I5,2X,A5,A5,I5,2x,3(F5.2,2X),
     *                  2(F7.3,1X))') RESNAM(I), ATNAM(I), RESNUM(I),
     *                  RESNAM(J),ATNAM(J), RESNUM(J), CHARGEI(I,1),
     *                  CHARGEI(J,1), RIJ(I,J), DGIJ(III,J),
     *                  DELTAIJ(III,J)
                ENDIF
599           CONTINUE
            ENDIF
591       CONTINUE
   
        CHARGEZ(LL)=ABS(CHARGEI(I,1))*(1-BURY(I))+CHARGEZ(LL)
         
590     CONTINUE                                              
        IONTOT=IONPOS+IONNEG                                  
      WRITE(20,'("PH=",F5.2,8X,"IONIC STRN. =",F12.8,"  IONIC CONC. =",
     *      F12.8)') PH, IONST, IONC  
      WRITE(20,'(13(F6.2))') (CHARGEI(JJ,1),JJ=1,M)               
      MCHSPH(LL) = TCHISPH                              
      MDGSPH(LL) = SDGISPH                             
      ITNUMT(LL) = ITNUM                              
      MIONPOS(LL) = IONPOS
      MIONNEG(LL) = IONNEG
      MIONTOT(LL) = IONTOT
      MCHRTOT(LL) = IONTOT + TCHISPH
 510  CONTINUE                                     

C***********************************************************************
C...                                                                
C...  THIS IS THE END OF THE LOOP THAT DOES THE ITERATION.         
C...                                                              
C***********************************************************************

      WRITE(50,'(//)')
      IF(ANSWER2.EQ.'Y'.OR.ANSWER3.EQ.'Y')
     *              WRITE(50,'(10X,"AA",14X,"RES",14X,"PK12")')
C      WRITE(30,'(////"IONIC STRN.",F12.8,"   IONIC CONC. =",F12.8,/)')
C     *      IONST, IONC
      IF(GRPHNUM.NE.0) THEN 
C        WRITE(30,'(/,1X,"NAME",2X,"NUM",3X," PH",3X, "  CHARGE",1X,
C     *    " DELTA G",3X,"pKi",3X,"DELTA pK",4X,"Ka  ",4X,"pKion-pK12")')
      ENDIF

      DO 620 I=1,M                                                
        DO 621 III=1,GRPHNUM
          IF((GNUM(III).EQ.RESNUM(I)).AND.(RESNAM(I)(1:3).EQ.GNAM(III))
     *  .AND.((RESNAM(I)(5:5).EQ.GCHN(III)).OR.(GCHN(III).EQ.'Z'))) THEN
            DO 622 J=1,NN
              IF(PKINT(I).EQ.0.0) THEN
C                WRITE(30,'(A5,I5,1X,F5.2,4X,F5.2,4X,F5.2,3X,F5.2,3X,
C     *               F5.2,1X,F9.2,4X,F5.2)') RESNAM(I),RESNUM(I),SPH(J),
                WRITE(30,'(F8.4,I5,1X,F5.2,4X,F5.2,4X,F5.2,3X,F5.2,3X,
     *               F5.2,1X,F9.2,4X,F5.2)') IONST,RESNUM(I),SPH(J),
     *               GCHARGI(III,J),GDGISPH(III,J),GP(III,J),
     *               GDELTPK(III,J),GASSCON(III,J),GLOGK(III,J)
              ELSE
C                WRITE(30,'(A5,I5,1X,F5.2,4X,F5.2,4X,F5.2,3X,F5.2,3X,
C     *                F5.2,14X,F5.2)')RESNAM(I),RESNUM(I),SPH(J),
                WRITE(30,'(F8.4,I5,1X,F5.2,4X,F5.2,4X,F5.2,3X,F5.2,3X,
     *                F5.2,1X,F9.2,4X,F5.2)')IONST,RESNUM(I),SPH(J),
     *                GCHARGI(III,J),GDGISPH(III,J),GP(III,J),
     *                GDELTPK(III,J),PK12(I),PK12(I)
              ENDIF
622         CONTINUE
          ENDIF
621     CONTINUE
        IF(ANSWER2.EQ.'Y') THEN
          IF(ANSWER3.EQ.'Y') THEN 
            IF((RESNAM(I)(1:3).EQ.'HIS').OR.(RESNUM(I).EQ.1)) THEN 
              WRITE(50,'(10X,A5,10X,I3,10X,F8.4)')RESNAM(I),RESNUM(I),
     *             PK12(I)
            ENDIF
          ELSE
            WRITE(50,'(10X,A5,10X,I3,10X,F8.4)') RESNAM(I),RESNUM(I),
     *           PK12(I)
          ENDIF
        ENDIF
620   CONTINUE                                                  
C      WRITE(30,'("      ITERA  CATIONS  ANIONS  IONS      H+  "
C     *      "  TOTAL   DG-TOTAL  DG-H+  DG-IONS",/," pH   TIONS   BOU"
C     *      "ND   BOUND   BOUND   BOUND   CHARGE  KCAL/MOL ")')
      WRITE(50,'(/T10,"PROTEIN CHARGE AND SUMMED FREE ENERGY AS A ", 
     *     "FUNCTION OF PH",/,8X,"IONIC STRN. =",F12.8,
     *     "  IONIC CONC. =",F12.8,/)') IONST, IONC
      WRITE(50,'("      ITERA  CATIONS  ANIONS  IONS      H+  "
     *      "  TOTAL   DG-TOTAL  DG-H+  DG-IONS",/," pH   TIONS   BOU"
     *      "ND   BOUND   BOUND   BOUND   CHARGE  KCAL/MOL ")')

        DO 600 I=1,NN                                           
C         WRITE(30,'(F5.2,2X,I3,3X,F6.2,2X,F6.2,2X,F6.2,2X,F6.2,2X,F6.2,
C     *        3X,F7.2,2X,F7.2,1X,F7.2)')SPH(I),ITNUMT(I),MIONPOS(I),
C     *        MIONNEG(I),MIONTOT(I),MCHSPH(I),MCHRTOT(I),MDGSPH(I)
         WRITE(50,'(F5.2,2X,I3,3X,F6.2,2X,F6.2,2X,F6.2,2X,F6.2,2X,F6.2,
     *        3X,F7.2,2X,F7.2,1X,F7.2)')SPH(I),ITNUMT(I),MIONPOS(I),
     *        MIONNEG(I),MIONTOT(I),MCHSPH(I),MCHRTOT(I),MDGSPH(I),
     *        CHARGEZ(I)
        WRITE(70,'(F5.2,2X,F8.4,3X,F6.2,2X,F6.2,2X,F6.2,2X,F6.2,2X,F6.2,
     *        3X,F7.2,2X,F7.2,1X,F7.2)')SPH(I),IONST,MIONPOS(I),
     *        MIONNEG(I),MIONTOT(I),MCHSPH(I),MCHRTOT(I),MDGSPH(I),
     *        CHARGEZ(I)
 600    CONTINUE                                              
        WRITE(50,'(///)')
      RETURN                                                 
      END                                                   
