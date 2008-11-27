C      PROGRAM STATACC (INFILE1,INFILE2,OUTFILE1,OUTFILE2,INFILE3)
C...   Program call commented out in i686 Linux port (MJH, 5/3/05)
C**********************************************************************
C...                                                                  
C...   This program calculates the static accessibilities of atoms in
C...   proteins and nucleic acids following the methodology of Lee and
C...   Richards (J. Mol. Biol. 1971, vol 55, pg 379) as implemented by
C...   J. B. Matthew.
C...
C...   Variable names:
C...   RESNAM    residue name
C...   RESNUM    residue number
C...   ATNAM     atom name
C...   ATNUM     atom number
C...   PRBRAD    radius of probe
C...   DISLIM    largest possible distance between two atoms
C...   DIS       distance between 2 atoms
C...   RADIUS    effective radius of an atom
C...   ATL       list of atom names
C...   RADL      list of van der Waals radii for atoms in ATL
C...   PACS      percent static accessibility
C...   ACCESIB   accessible surface area
C...   NUMSLIC   number of slices across the central atom
C...   YRADIUS   radii of circles defined by atoms crossed by planes
C...   DEG       degrees
C...   PKL       list of pKint that correspond to various charged at.
C...   SASL      list of SA values used for normalization
C...
C...   The following files are used by this program:
C...   TAPE 10 = .chr
C...   TAPE 20 = .cor
C...   TAPE 30 = .saout
C...   TAPE 50 = .elc
C...   TAPE 100 = input from the UNIX shell
C...
C...   Bertrand Garcia-Moreno E.
C...
C...   May 1988
C...
C...   Fortran 77 on a HP 9000-825SRX
C...
C**********************************************************************

C        Copyright 1988, Bertrand Garcia-Moreno E.
C        This program is distributed under General Public License v. 3.  See the
C        file COPYING for a copy of the license.  

      DIMENSION RESNUM(200),RESNAM(200),ATNAM(200),X(200),Y(200),   
     *          Z(200),RADIUS(200),DIS(200),ATL(15),RADL(15)
      INTEGER RESNUM, NUM
      CHARACTER RESNAM*5,ATNAM*5,HEADER*80,ATL*4,INFILE1*20,INFILE2*20,
     *          INFILE3*20,OUTFILE1*20,OUTFILE2*20            
      REAL DIS,X,Y,Z,PRBRAD,DISLIM,RADIUS,RADL,PACS,ACCESSIB

      DATA ATL(1),ATL(2),ATL(3),ATL(4),ATL(5),ATL(6),ATL(7),ATL(8),
     *ATL(9),ATL(10),ATL(11),ATL(12),ATL(13), ATL(14),ATL(15),
     *RADL(1),RADL(2),RADL(3),RADL(4),RADL(5),RADL(6),RADL(7),RADL(8),
     *RADL(9),RADL(10),RADL(11),RADL(12),RADL(13),RADL(14),RADL(15)
     */' O  ',' N  ',' CA ',' CAL',' S  ','FE ','CL ','CU ','K  ',
     *'NH4','ZN ','MG ','MN ','NA ',' O',
     *1.52, 1.55, 1.70, 0.99, 2.19, 0.74, 1.81, 0.72, 1.33, 1.43,
     *0.74, 0.66, 0.80,0.99,1.50/	

C......................................................................
C...
C...  Tape 10 contains the coordinates of the atoms whose accessibility
C...  we are calculating.
C...  Tape 20 contains the coordinates of the structure to which the
C...  atoms listed in Tape 10 belong.
C...  Tape 30 is the file on which the output listing the contacts
C...  defining the accessibilities is written.
C...  Tape 40 is the file created by the UNIX shell interactive inter-
C...  face used to drive this program. It is used to pass parameter
C...  values to the program.
C...  Tape50 is the file created by this program which becomes the 
C...  input to satkel.f
C...
C......................................................................
  
      PRBRAD=1.4
      DIS(1)=0.0
C     SURFSUM=0.0   
   
C...  Following 5 lines of code added by Mike Harms, 12/09/2004
      CALL GETARG(1,INFILE1)
      CALL GETARG(2,INFILE2)
      CALL GETARG(3,OUTFILE1)
      CALL GETARG(4,OUTFILE2)
      CALL GETARG(5,INFILE3)

      OPEN(UNIT=10,FILE=INFILE1,STATUS='OLD',IOSTAT=IOS)
      OPEN(UNIT=20,FILE=INFILE2,STATUS='OLD',IOSTAT=IOS)
      OPEN(UNIT=30,FILE=OUTFILE1,STATUS='NEW',IOSTAT=IOS)
      OPEN(UNIT=50,FILE=OUTFILE2,STATUS='NEW',IOSTAT=IOS)
      OPEN(UNIT=100,FILE=INFILE3,STATUS='OLD',IOSTAT=IOS)
C     READ(100,'(A1)') ANSWER
C     IF((ANSWER.EQ.'n').OR.(ANSWER.EQ.'N')) READ(100,*) PRBRAD

C...  Above line causing compiler crash using g77, i686 linux.  Commented out.
C...  As far as I can tell, this only looks for alternate radius sizes.  
C...  MJH - 5/03/2005

101   CONTINUE

    
      DISLIM= 2*(2.2 + PRBRAD)

C...................................................................... 
C...                                                            
C...  Do 500 loop reads and writes headers
C...                                                          
C......................................................................

        DO 500 I=1,5                                         
        READ(10,'(A80)') HEADER                                
        WRITE(30,'(A80)') HEADER
        WRITE(50,'(A80)') HEADER                              
500     CONTINUE                                         

C.......................................................................
C...                                                               
C...  The program loops around the 100   READ statement for every atom
C...  whose accessibility we are going to calculate. The coordinates of
C...  the "central" atom (the atoms whose SA we are to calculate) are
C...  stored as the first element of the X, Y, Z, etc arrays.
C...
C.......................................................................

100   READ(10,9000,END=9999) RESNUM(1),RESNAM(1),ATNAM(1),X(1),Y(1),Z(1)
      IF((RESNAM(1).EQ.'     ').OR.(ATNAM(1).EQ.'     ')) GO TO 100
      WRITE(30,9010)RESNUM(1),RESNAM(1),ATNAM(1)
      WRITE(30,9020)                                         
      CALL DISTANS(RESNUM,RESNAM,ATNAM,X,Y,Z,NUM,DIS,DISLIM)       

C.......................................................................
C...                                                              
C...  Subroutine Distans has found all the atoms around the central
C...  atom that will influence its SA. The atoms are now assigned an
C...  effective radius equal to the sum of the atoms van der Waals
C...  radius and the radius of the probe. The default value of the
C...  probe radius is 1.4, the radius of a water molecule.
C...
C.......................................................................

      DO 510 I=1,NUM                                              
      RADIUS(I)= 1.8 + PRBRAD
      
      DO 520 II=1,15
      IF(ATNAM(I)(1:4).EQ.ATL(II)) RADIUS(I)=RADL(II) + PRBRAD
C      IF((ATNAM(I)(1:2).EQ.ATL(II)).AND.(RESNAM(I).EQ.'PTP ')) 
C     *    RADIUS(I)=RADL(II) + PRBRAD
C     IF(RESNAM(I).EQ.'PTP  ') RADIUS(I)=1.52 + PRBRAD
520   CONTINUE
      IF(RESNAM(I)(1:3).EQ.'PTP') RADIUS(I)=1.52 + PRBRAD
      IF((ATNAM(I)(1:4).EQ.' N  ').AND.(RESNUM(I).EQ.1))RADIUS(I)=1.8+
     *                                                  PRBRAD

C......................................................................
C...                                                           
C...  Polar (noncarbon) atoms that are within 3.7 angstroms of the
C...  central atom are identified.
C...                                                           
C......................................................................

        IF(DIS(I).LE.(RADIUS(1)+RADIUS(I))) THEN
          IF((DIS(I).LT.3.7).AND.(RESNUM(1).NE.RESNUM(I))) THEN  
            IF((ATNAM(I)(2:3).EQ.'CL').OR.
     *         (ATNAM(I)(2:3).EQ.'CU').OR.
     *         (ATNAM(I)(2:4).EQ.'CAL')) THEN                 
              WRITE(30,9030) RESNUM(I), RESNAM(I), ATNAM(I),       
     *                   X(I), Y(I), Z(I), RADIUS(I), DIS(I)   
            ELSE IF(ATNAM(I)(2:2).EQ.'C') THEN                
              WRITE(30,9040) RESNUM(I), RESNAM(I), ATNAM(I),    
     *                   X(I), Y(I), Z(I), RADIUS(I), DIS(I)
            ELSE                                           
              WRITE(30,9030) RESNUM(I), RESNAM(I), ATNAM(I),           
     *                    X(I), Y(I), Z(I), RADIUS(I), DIS(I)      
            ENDIF                                                 
          ELSE                                                   
            WRITE(30,9040) RESNUM(I),RESNAM(I),ATNAM(I),X(I),Y(I),Z(I),
     *                RADIUS(I),DIS(I)                               
          ENDIF                                                     
        ENDIF
510   CONTINUE                                                     

C.......................................................................
C...
C...   The subroutines that calculate the accessibility and that nor-
C...   malize the accessibility data are called.
C...   The data is written and the files are closed. So long!
C...
C.......................................................................

      CALL ACSCALC(NUM,X,Y,Z,RADIUS,PACS,ACCESIB)          
      CALL ACSNORM (NUM,RESNUM,RESNAM,ATNAM,X,Y,Z,PACS,BURY)  
C     SURFSUM=SURFSUM+ACCESIB                                      


      WRITE (30,9055)                                                  
      WRITE (30,9060) ACCESIB                                        
      WRITE (30,9070) PACS                                          
C     WRITE (30,9071) SURFSUM                                     
      WRITE (30,9072) BURY
      GO TO 100                                               

9999  CONTINUE                                              
9000  FORMAT(5X,I4,3X,A5,3X,A5,2X,3(2X,F8.3),6X,I5)
9010  FORMAT('************************',2X,I4,2X,A5,2X,A5,2X,
     *'***********************')                             
9020  FORMAT(25X,'X',8X,'Y',10X,'Z',7X,'R(EFF)',4X,'R(IJ)',/)
9030  FORMAT(I4,2X,A5,2X,A5,5(2X,F8.2),' HB OR SB')   
9040  FORMAT(I4,2X,A5,2X,A5,5(2X,F8.2))             
9055  FORMAT('***************************************************',
     *         '***************',//)                              
9060  FORMAT(3X,'THE TOTAL ACCESSIBLE SURFACE AREA IS=',F9.2)    
9070  FORMAT(3X,'THE PERCENT ACCES.   SURFACE AREA IS=',F9.2)    
9071  FORMAT(3X,'THE ACCUMULATIVE SURFACE AREA IS=',F12.2,/////) 
9072  FORMAT(3X,'THE    (1-SA)    PARAMETER  IS      =',F9.2,/////)
9090  FORMAT('YOU MUST MUST SEARCH FOR DIS LARGER THAN 6.4')    
      STOP                                                     
      END                                                     



      SUBROUTINE DISTANS(RESNUM,RESNAM,ATNAM,X,Y,Z,NUM,DIS,DISLIM)  

C***********************************************************************
C...
C...   This subroutine pulls out all the atoms that are within a
C...   distance of the central atoms. The distance is defined by
C...   the sum of the probe radius and the radius of the largest
C...   atom that enters the calculation.
C...                                                              
C***********************************************************************

      DIMENSION RESNUM(200),RESNAM(200),ATNAM(200),X(200),Y(200),   
     *          Z(200),DIS(200)                        
      CHARACTER RESNAM*5,ATNAM*5,HEADER*80              
      REAL X, Y, Z, DISLIM, DIS, DX, DY, DZ
      INTEGER RESNUM,K,NUM
      K=2                                          

C......................................................................
C...                                                           
C...  Tape20 is accessed once per central atom.
C...                                                         
C......................................................................

        DO 500 I=1,5                                        
        READ(20,'(A80)') HEADER                               
500     CONTINUE      

C......................................................................
C...                   
C...   Distances are checked and atoms to be used in the calculation
C...   of the SA are selected. Contacts of less than 2.0 angstroms
C...   between atoms are identified.
C...
C.......................................................................

100   READ(20,9000,END=9998)RESNUM(K),RESNAM(K),ATNAM(K),X(K),Y(K),Z(K)
      IF(K.EQ.199) WRITE(30,9010) RESNUM(1),RESNAM(1),ATNAM(1)        
      IF((RESNUM(K).EQ.RESNUM(1)).AND.(ATNAM(K).EQ.ATNAM(1)).AND.
     *(RESNAM(K).EQ.RESNAM(1)))GO TO 100
      IF((RESNAM(K).EQ.'     ').AND.(ATNAM(K).EQ.'     ')) GO TO 100
      DX=X(K)-X(1)                               
      DY=Y(K)-Y(1)                              
      DZ=Z(K)-Z(1)                             
      DIS(K)=SQRT(DX**2+DY**2+DZ**2)          
      IF(DIS(K).GT.DISLIM) GO TO 100
      IF((DIS(K).LT.2.0).AND.(RESNUM(1).NE.RESNUM(K))) WRITE(30,9020)
     *RESNUM(1), RESNUM(K)
      K=K+1                                                   
      GO TO 100                                              

9998  CONTINUE      

C.......................................................................
C...                                                               
C...    NUM is a variable that describes the number of atoms found to
C...    influence the SA of the central atom
C...                                                              
C.......................................................................

      NUM=K-1                                                       

9000  FORMAT(5X,I4,3X,A5,3X,A5,2X,3(2X,F8.3))
9010  FORMAT('BEWARE',2X,I4,2X,A5,2X,A5,2X,'HAS MORE THAN 199 CON', 
     *        'TACTS LESS THAN 6.4 ANGSTROMS')                    
9020  FORMAT('CLOSE CONTACT BEWTEEN RES',1X,I4,1X,'AND RES',1X,I4) 
      REWIND(20)                                                 
      RETURN                                                    
      END   



      SUBROUTINE ACSCALC(NUM,X,Y,Z,RADIUS,PACS,ACCESIB)  
C***********************************************************************
C...                                                                   
C...     This subroutine is the one that actually performs the calcula-
C...   tion of static accessibilities according to the methodology of
C...   Lee and Richards, following the algorithm of Matthew.
C...                                                               
C***********************************************************************
 
      DIMENSION X(200),Y(200),Z(200),RADIUS(200),YRADIUS(200)   
      INTEGER FLAG,NUMSLIC,NUM
      REAL SLICE,X,Y,Z,PACS,ACCESIB,RADIUS,YRADIUS,YC,Y1,DEG,D,XT,ZT,
     *ACCES, ACC
      ACCES=0.0                                              
      FLAG=0     

C.......................................................................
C...                                                          
C...   The number of slices across the central atoms assuming slice
C...   thickness of 0.25 angs is calculated. Slicing is done along the
C...   Y axis. Y1 refers to the value of Y(1) as one moves to slices
C...   away from the central slice. YC refers to the distance from the
C...   central slice to the other slices
C...
C.......................................................................

200   CONTINUE                                                    
      IF(FLAG.EQ.1) THEN                                         
        NUMSLIC=0                                               
        SLICE=(RADIUS(1)/.25)                                      
        NUMSLIC=IFIX(SLICE)                                       
        Y1=Y(1) - .25                                        
        YC=.25                                              
      ELSE                                                 
       SLICE=(RADIUS(1)/.25)+1.0                             
        NUMSLIC=IFIX(SLICE)                                  
        Y1=Y(1)                                         
        YC=0.0                                         
      ENDIF                                           

C.......................................................................
C...                                                               
C...   The DO 500 is the heart of the program. The length of arc
C...   available to contact with water is calculated for each slice .
C...   The first time around the DO 500 loop the lengths of arc on one
C...   hemisphere of the atom are calculated. The second time around
C...   the ones on the other hemisphere are calculated.
C... 
C.......................................................................

210     DO 500 N=1,NUMSLIC                                        
        IF(YC.GE.RADIUS(1)) GO TO 500     

C.......................................................................
C...                                                            
C...    The radius of the circles defined by the slicing of the atoms 
C...    in each slice are calculated in the DO 510 loop
C...                                                         
C.......................................................................

          DO 510 J=1,NUM                                    
          IF(ABS(Y(J)-Y1).GT.RADIUS(J)) GO TO 100                    
          YRADIUS(J)=SQRT(RADIUS(J)**2-ABS(Y(J)-Y1)**2)             
          GO TO 510                                     
100       YRADIUS(J)=0.0                               
510       CONTINUE                                    

C.......................................................................
C...
C...   Coordinates (XT,ZT) at each of 360 points on the perimeter of
C...   central atom at slice N  are calculated. If the distance between
C...   this point and the center of the surrounding atoms is larger
C...   than the radius of the surrounding atoms, the length of arc
C...   equivalent to a degree is available for contact with water.
C...
C.......................................................................

          DEG=0.0                                    
          D=0.0                                     

          DO 520 J=1,360                                        
          D=D+1.0                                              
          XT=YRADIUS(1)*COS(D)+X(1)                           
          ZT=YRADIUS(1)*SIN(D)+Z(1)                          
            DO 530 L=2,NUM                                  
            IF(SQRT((XT-X(L))**2+(ZT-Z(L))**2).LE.YRADIUS(L)) GO TO 520
530         CONTINUE                                    
          DEG=DEG+1.0                                  
520       CONTINUE                                    

C.......................................................................
C...
C...   The accessible surface area is calculated according to the for-
C...   mula given by Richards or by Matthews (pg 412 Biochem. Biophys.
C...   Res. Comm. (1978) vol 81, #2)
C...
C.......................................................................

          ACC=(RADIUS(1)*.25*(3.1416*YRADIUS(1)*DEG/180.)/
     *        SQRT(RADIUS(1)**2-(YC**2)))
          IF(YRADIUS(1).EQ.0.0) ACC=0.0
          ACCES=ACCES+ACC                           

C..............................................................  
C...
C...   The calculations proceeds on to the next slice.
C...                                                      
C..............................................................      

          IF(FLAG.EQ.0) THEN
            Y1=Y1+.25                                           
            YC=YC+.25                                          
                IF(N.EQ.NUMSLIC) THEN                        
                  FLAG=1                                    
                  GO TO 200                                
                ENDIF                                    
          ELSE                                            

C.......................................................................
C...                                                            
C...   This part of the program integrates arc lengths on the other
C...   hemisphere of the atom
C...                                                        
C......................................................................


            Y1=Y1-.25              
            YC=YC+.25             
          ENDIF                  
500     CONTINUE                

C.......................................................................
C...                                                              
C...   The total accessible surface area of an atom is expressed in
C...   terms of percent accessible surface area.
C...
C.......................................................................

      PACS=(100.0*ACCES)/(4.0*3.1416*(RADIUS(1)**2))     
      ACCESIB=ACCES                                     
      RETURN                                              
      END               




      SUBROUTINE ACSNORM(NUM,RESNUM,RESNAM,ATNAM,X,Y,Z,PACS,BURY)
C.......................................................................
C...                                                            
C...   This subroutine normalizes the accessibilities calculated in
C...   ACSCALC, it assigns a pKint to a charged residue, and it writes
C...   TAPE50, which after editing becomes the input file used by
C...   satkel.f in the calculation of electrostatic interactions.
C...   In order to add new atoms to this subroutine, modify the DATA
C...   statement by increasing the ATL, PKL and SASL arrays and by 
C...   raising the limits in the DO 500 loop and the DIMENSION
C...   statements accordingly.
C...
C.......................................................................

      DIMENSION RESNUM(200),RESNAM(200),ATNAM(200),X(200),Y(200),     
     *          Z(200),ATL(17),PKL(17),SASL(17)
      INTEGER RESNUM,FLAG
      REAL X, Y, Z, PKL, SASL, CHARGE, PKINT
      CHARACTER RESNAM*5,ATNAM*5,ATL*4
      DATA ATL(1),ATL(2),ATL(3),ATL(4),ATL(5),ATL(6),ATL(7),ATL(8),
     *ATL(9),ATL(10),ATL(11),ATL(12),ATL(13),ATL(14),ATL(15),ATL(16),
     *ATL(17),
     *PKL(1),PKL(2),PKL(3),PKL(4),PKL(5),PKL(6),PKL(7),PKL(8),PKL(9),
     *PKL(10),PKL(11),PKL(12),PKL(13),PKL(14),PKL(15),PKL(16),PKL(17),
     *SASL(1),SASL(2),SASL(3),SASL(4),SASL(5),SASL(6),SASL(7),SASL(8),
     *SASL(9),SASL(10),SASL(11),SASL(12),SASL(13),SASL(14),SASL(15),
     *SASL(16),SASL(17)
     */' OD1',' OD2',' NZ ',' ND1',' NE2',' OE1',' OE2',' OH ',' NH2',
     *' NH1',' OXT',' O  ',' O1',' O2',' SG ','CL ','NA ',
     *4.00,4.00,10.40,6.00,6.60,4.50,4.50,10.0,12.0,12.0,3.60,3.60,
     *4.00,4.00,9.00,0.00,0.00,
     *50.7,50.7,55.8,27.6,31.2,50.9,50.9,51.6,51.3,51.3,00.0,00.0,
     *46.6,46.6,44.2,100.00,100.00/
      FLAG=0
        IF((RESNUM(1).EQ.1).AND.(ATNAM(1).EQ.' N  '))  THEN   
          PKINT= 7.00                                            
          CHARGE=1.00                                           
          WRITE(50,9000) RESNUM(1),RESNAM(1),ATNAM(1),X(1),Y(1), 
     *              Z(1),CHARGE,PKINT,'0000'                       
          RETURN                                                   
        ENDIF

          DO 500 I=1,17
          IF(ATNAM(1)(1:4).EQ.ATL(I))THEN
            PKINT=PKL(I)
            SASTAND=SASL(I)
            IF(ATNAM(1)(2:2).EQ.'N') CHARGE=1.00
            IF(ATNAM(1)(2:2).EQ.'O') CHARGE=-1.0
            IF(ATNAM(1)(1:2).EQ.'CL') CHARGE=-1.0
            IF(ATNAM(1)(1:2).EQ.'NA') CHARGE= 1.0
            FLAG=1
          ENDIF
500       CONTINUE     
          
          IF(FLAG.NE.1) THEN
          WRITE(50,9010) RESNUM(1),RESNAM(1),ATNAM(1),X(1),   
     *                    Y(1), Z(1) , '0000', '00000','0000'
          RETURN                                            
          ENDIF                                              
C...                                                      
        IF(PACS.LE.0.05) THEN                            
          BURY=0.99                                     
        ELSE                                           
          IF(SASTAND.EQ.0.0) THEN
            BURY=0.0
          ELSE
            BURY=1.0-(PACS/SASTAND)                     
            IF(BURY.GT.0.99) BURY=0.99                 
            IF(BURY.LE.0.01) BURY=0.01                
          ENDIF                                      
        ENDIF
      WRITE (50,9020) RESNUM(1), RESNAM(1), ATNAM(1), X(1), Y(1),    
     *               Z(1), CHARGE, PKINT, BURY                      
9000  FORMAT(1X,I4,2X,A5,1X,A5,2X,3(1X,F8.3),2X,F4.1,1X,F5.2,      
     *       1X,A4)                                               
9010  FORMAT(1X,I4,2X,A5,1X,A5,2X,3(1X,F8.3),2X,A4,1X,A5,        
     *       1X,A4)                                             
9020  FORMAT(1X,I4,2X,A5,1X,A5,2X,3(1X,F8.3),2X,F4.1,1X,F5.2,  
     *        1X,F4.2)                                        
      RETURN                                                 
      END                                                   
