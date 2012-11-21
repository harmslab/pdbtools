       PROGRAM MKWIJ INFILE1
C*********************************************************************
C...
CMJH    NOTES ON LINUX (x86) PORT 
CMJH    This code as originaally written on an HP 9000-825.  At some
CMJH    point prior to May, 2005 it was recompiled for SGI using 
CMJH    the MIPSpro 7.2.1 F77 compiler.  This version has been 
CMJH    slightly altered to allow it to compile in a Linux x86 
CMJH    environment using gcc 3.4.3.  A comparison of the output from
CMJH    the SGI and linux binaries is in the test/ folder.  No  
CMJH    original code has been deleted, but rather commented out with
CMJH    the 'CMJH' comment line.  Specific changes are noted in line
CMJH    with the same comment line.  
CMJH
CMJH    Michael J. Harms, May 19, 2005
C...
C...    MKWIJ uses the model for dielectric properties of proteins by
C...    Tanford and Kirkwood to calculate the electrostatic free energy
C...    of interaction between a pair of positive unit charges. 
C...    For a description of the model see J. Am. Chem. Soc. (1957) 79,
C...    5333-5339  and references therein.
C...    
C...    Parameters required by the algorithm:
C...
C...    DINT  - value for the dielectric constant inside the protein
C...    DEXT  - value for the dielectric constant outside the protein
C...    TEMP  - temperature in degrees Kelvin
C...    N     - number of ionic strengths to be explored
C...    SRION - value of a particular ionic strength (Molar)
C...    B     - the radius of the protein's sphere equivalent in
C...            angstroms. The largest value of B that the program can
C...            handle is 40.0 angstroms.           
C...    A     - the ionic exclusion radius; A=B+2.0 angstroms
C...    G     - the depth of burial of charges relative to the
C...            interface betweena the regions of DEXT and DINT
C...            For purposes of calculations with the SA-TK
C...            algorithm, G=0.0.
C...            
C...   Peculiarities of the algorithm:
C...   
C...   The code that follows was written in ANSI-standard FORTRAN 77
C...   for an HP 9000-825 computer. The program need not operate
C...   with variables of double precision in machines with word size
C...   larger than 32 bits. Double precision in 32 bit machines is
C...   necessary for reasons stated in documentation below. The 
C...   intrinsic fortran functions used in the program might not be
C...   available in the same way in all Fortran compilers: 
C...   IDINT   -double prec. to integer conversion
C...   DNINT   -double prec. nearest whole number
C...   DFLOAT  -integer to double prec. conversion
C...   DABS    -absolute value of a double prec. variable
C...   DSQRT   -square root of a double prec. variable
C...   DLOG    -natural log of a double prec. variable
C...
C...   Default Parameters:
C...
C...   The default values for the parameters of the algorithm are
C...   those required by calculations on hemoglobin.
C...
C...
C...  February 1988
C... 
C...  Code by Bertrand Garcia-Moreno E. apres program ELECT written by
C...  Steve J. Shire and George H. Hanania in the laboratory of
C...  Frank R. N. Gurd.
C...
C...  Bertrand Garcia-Moreno E.
C...  Department of Biology
C...  The Johns Hopkins University
C...  Baltimore, MD 21218 (301) 339-7239
C...
C***********************************************************************

C        Copyright 1988, Bertrand Garcia-Moreno E.
C        This program is distributed under General Public License v. 3.  See the
C        file COPYING for a copy of the license.  

      IMPLICIT DOUBLE PRECISION (A-G, O-Z),INTEGER (I-N)
      CHARACTER HPROTEIN*20, FILENAM*20,DEFAULT*3, INFILE1*20, 
     *          SALTNAM*7
      INTEGER ANFCTR, CAFCTR
      DIMENSION WIJ(20,800),CTHETA(800),CIJ(800),RIJ(800),SRION(20),
     *          SCONC(20),SALTCONC(2,20),ANFCTR(2),CAFCTR(2), 
     *          SALTNAM(2)

      DATA ((WIJ(I,J),J=1,800),I=1,20),
     *     (CTHETA(I),I=1,800),
     *     (CIJ(I),I=1,800),
     *     (RIJ(I),I=1,800),
     *     (SRION(I),I=1,20)/16000*0.0,800*0.0,800*0.0,800*0.0,20*0.0/

C***********************************************************************
C...
C...   The default values of various parameters are those used in
C...   calculations of hemoglobin.
C...
C***********************************************************************

      HPROTEIN='HEMOGLOBIN'
      SRION(1)=0.00
      SRION(2)=0.01
      SRION(3)=0.10
      SRION(4)=0.18
      DINT=4.0                                                  
      DEXT=78.5                                                 
      TEMP=298.16                                               
      B=27.0
      A=29.0
      G=0.0
      N=4

C***********************************************************************
C...
C...   If the default values are not going to be used, new values
C...   are read
C...
C***********************************************************************
      
CMJH  The CALL argument had to be added because the PROGRAM statement 
CMJH  was commented out.  
      CALL GETARG(1,INFILE1)

CMJH  This forces the first salt concentration to be 0.0.  This is
CMJH  required because of the problem initializing SRION on line ~148.
      SRION(1) = 0.0
      SCONC(1) = 0.0
      SALTCONC(1,1) = 0.0
      SALTCONC(2,1) = 0.0
         
      OPEN (UNIT=10,FILE=INFILE1,STATUS='OLD')
        READ(10,'(A20)') FILENAM
      OPEN (UNIT=20,FILE=FILENAM,STATUS='NEW',IOSTAT=IOS)
CMJH  IF(IOS.NE.0) STOP 'Error opening Tape20'
CMJH  Crashes for some reason.  Just a file existance check.  
      READ(10,'(A3)') DEFAULT
      IF((DEFAULT(1:1).EQ.'y').OR.(DEFAULT(1:1).EQ.'Y')) GO TO 100 
      READ (10,*) HPROTEIN
      READ (10,*) DINT, DEXT, TEMP
      READ (10,*) A, B, G
      READ (10,*) N
        N = N + 1
CMJH    SRION(I) = 0.0
CMJH    For some reason this variable initialization casused
CMJH    segmentation fault.  Maybe array overrun? 
        DO 110 I = 2,N
        READ (10,'(F12.8,1X,F12.8,2(1X,A7,1X,F10.7,1X,I1,1X,I1))')
     *        SRION(I), SCONC(I), SALTNAM(1), SALTCONC(1,I),  ANFCTR(1),
     *        CAFCTR(1), SALTNAM(2), SALTCONC(2,I), ANFCTR(2), CAFCTR(2)
110     CONTINUE
      IF(B.GT.40.0) STOP 'Protein radius too large for arrays'
      IF (N.GT.20) STOP 'Too many ionic strengths for arrays'
100   CONTINUE
      WRITE (20,*) 'Wij factors from Tanford-Kirkwood algorithm'
      WRITE (20,*) 'calculated for the study of ',HPROTEIN
      WRITE (20,*) 'The parameters for the calculations are:'
      WRITE (20,'("Dint = ",F5.1)') DINT
      WRITE (20,'("Dext = ",F5.1)') DEXT
      WRITE (20,'("Temp = ",F6.2)') TEMP
      WRITE (20,'("A= ",F5.2," angstroms")') A
      WRITE (20,'("B= ",F5.2," angstroms")') B
      WRITE (20,'("G= ",F5.2," angstroms")') G
      WRITE (20,'("The Wij factors were calculated at ",I3," ionic",
     *             " strengths (M):")') N
        DO 120 I=1,N
        WRITE (20,'(F12.8,1X,F12.8,2(1X,A7,1X,F10.7,1X,I1,":",I1))')
     *        SRION(I), SCONC(I), SALTNAM(1), SALTCONC(1,I),  ANFCTR(1),
     *        CAFCTR(1), SALTNAM(2), SALTCONC(2,I), ANFCTR(2), CAFCTR(2)
120     CONTINUE

C***********************************************************************
C...
C...  Angstroms converted into cms. The values of Wij will be calcula-
C...  ted in ergs, a unit of energy of particular value in the world
C...  of biophysico-chemical research!!!
C...
C***********************************************************************

      A = A * 1E-8
      B = B * 1E-8
      G = G * 1E-8
 
C**********************************************************************
C...
C...    Compute parameters used in the calculations of Bij and Cij.
C...
C***********************************************************************

      DIFF=B-G                                                         
      ROE=DIFF*DIFF/(B*B)                                              
      FACTA=B/DINT                                                     
      DELTA=DINT/DEXT                                                  
      SIGMA=DIFF*DIFF/(A*A)                                            

C***********************************************************************
C...
C...      Compute Aij and Bij terms in steps of 0.1 angstroms. Notice
C...   that Aij and Bij terms are ionic strength-independent. MM is the
C...   number of 0.1 angstrom steps in the length equivalent to the
C...   largest distance of interaction for any pair of values of B and
C...   G. 
C...
C***********************************************************************

      RIJ(1)=.1E-8                                                     
      RM=(DIFF*1.E+08)/.1                                              
      MM=IDINT(DNINT(RM))*2      
        DO 140 II=1,MM     
        CTHETA(II)=(2.*DIFF*DIFF-RIJ(II)*RIJ(II))/(2.*DIFF*DIFF)       
        AIJ=FACTA*(1./RIJ(II))                                         
        BIJ=1.  

C*********************************************************************
C...
C...   Iteration for the calculation of Bij terms. Conversion of the
C...   iteration is set at the level of the 5th significant figure.
C...
C*********************************************************************

          DO 150 I=1,500                                          
          IF(I.EQ.500) PRINT *,' Bij iteration for Rij= ',II/10.,
     *    ' did not converge'                    
          RI=DFLOAT(I)                                   
          YY=CTHETA(II)                                  
          BIJS=BIJ+(ROE**RI*FCNT(YY,I))/((RI+1.)*(RI+1.))  
          IF (DABS(BIJS-BIJ)-1.E-5) 160,160,151    
  151     BIJ=BIJS
  150     CONTINUE
  160   SCONT=DSQRT(1.-2.*ROE*CTHETA(II)+ROE*ROE)    
        BIJ=(BIJ*DELTA*DELTA/DINT)+((DELTA-3.*DELTA*DELTA)/(DINT*ROE))
     1  *DLOG((SCONT+ROE-CTHETA(II))/(1.-CTHETA(II)))+(1.-2.*DELTA+2.*
     2  DELTA*DELTA)/(DINT*SCONT)                      

C**********************************************************************
C...
C...        Wij terms at I=0.0 M are computed and stored. The value
C...        of the constant 2.30686E-19 (ergs*cm),  was obtained using:
C...        q=electronic charge=1.60209E-19 C
C...        e=permittivity=8.85419E-12 C2/Nm2
C...
C...        Wij = q**2/4*pi*e
C...
C...
C*********************************************************************

        WIJ(1,II)=((AIJ-BIJ)/B)*2.30686E-19                          
        RIJ(II+1)=RIJ(II)+.1E-8                                        
  140   CONTINUE                                                       
 
C**********************************************************************
C...
C...       Calculation of Cij terms for Wij factors at N nonzero 
C...       ionic strengths
C...       For your information, the value of kappa was computed with
C...       the following values for the physical constants:
C...
C...       q=1.60209E-19 C
C...       N=6.02252E23  molecules/mol
C...       V=1000 cm3
C...       k=1.38054E-23 J/K
C...       e=8.85419E-12 C2/Nm2
C...
C...       kappa= sqrt((8*pi(q**2)*N)/(4*pi*e*k*V))*sqrt(I/Dext*T)
C...
C**********************************************************************

          DO 170 K=2,N
          RION=SRION(K)   
          GKAPPA=(50.29/DSQRT(DEXT*TEMP))*DSQRT(RION)*1.E08  
          Q=GKAPPA*A 
          RIJ(1)=.1E-8    
            DO 180 II=1,MM
            SCIJ=Q/(1.+Q)   
            YY=CTHETA(II) 

C***********************************************************************
C...
C...   Iteration for computation of Cij terms. Criterion for conversion
C...   of the iteration is set at the level of the 5th significant
C...   figure. The upper limit of the DO 190 loop is set at 63 because
C...   when I=63, the number returned by function FUN is larger than
C...   the HP 9000-825 can handle.
C...  
C**********************************************************************

                DO 190 I=1,63  
                RI=DFLOAT(I)   
                IF(I.EQ.63) PRINT *,' Cij iteration for Rij= ',II/10.,
     *          ' did not converge'                    
                SCIJS=SCIJ+((2.*RI+1.)/(2.*RI-1.))*((DEXT/((RI+1.)*DEXT+
     1          (RI*DINT)))**2)*(Q*Q*(SIGMA**RI)*FCNT(YY,I))/((FUN(Q,I+1
     2          )/FUN(Q,I-1))+((RI*(DEXT-DINT))/((RI+1.)*DEXT+RI*DINT))*
     3          ((B/A)**(2.*RI+1.))*(Q*Q/(4.*RI*RI-1.)))               
                IF (DABS(SCIJS-SCIJ)-1.E-5) 191,191,192
192             SCIJ=SCIJS                                             
190             CONTINUE
191        CONTINUE
           CIJ(II)=SCIJS/DEXT
180        CONTINUE               

C***********************************************************************
C...
C...     Wij factors at non-zero ionic strength are calculated.
C...                                                       
C***********************************************************************

        DO  200 II=1,MM                                   
        WIJ(K,II)=WIJ(1,II)-(CIJ(II)/A)*2.30686E-19   
200     CONTINUE                                       
170   CONTINUE                                      

C**********************************************************************
C...
C...    Write the Wij factors unto tape 20 and exit the program
C...
C*********************************************************************

        DO 210 I=1,N                             
        WRITE(20,'(6E12.5)') (WIJ(I,J),J=1,MM)        
210     CONTINUE                              
      CLOSE (10)
      CLOSE (20)
      END

*********************   FUNCTION FCNT   ******************************
 
      DOUBLE PRECISION FUNCTION FCNT(Q,L)

C**********************************************************************
C...
C...     This function uses a recursive formula to calculate the 
C...     Legendre Polynomial of the nth order, Pn(x), where
C...     x=cos(theta). The special cases where n=1 (Pn=1) and when
C...     n=2 (Pn=cos(theta)) are incorporated explicitly.
C...                         
C*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-G, O-Z), INTEGER (I-N)
      P1=1.              
      P2=Q              
        IF ( L .EQ. 1 ) THEN
          FCNT=P2
          RETURN
        ELSE IF ( L .GT. 1) THEN
          DO 100 K = 2, L
            FL=DFLOAT(K)
            P3=((2.0*FL-1.0)*Q*P2-(FL-1)*P1)/(FL)
            P1=P2
            P2=P3
  100     CONTINUE
          FCNT=P2
          RETURN
        ELSE
          FCNT=1.0
          RETURN
        ENDIF
      END                            

********************   FUNCTION FUN   **********************************

      DOUBLE PRECISION FUNCTION FUN(Q,I)             

C*********************************************************************
C...
C...     This function is used to calculate the value of the polyno-
C...     mials employed in the computation of the Cij terms
C...                                                    
C***********************************************************************

      IMPLICIT DOUBLE PRECISION (A-G, O-Z), INTEGER (I-N)
      FUN=1.0                                       
      IF(I-1) 10,20,20                                
 10   RETURN                                      
 20     DO 100 J=1,I                                 
        RJ=DFLOAT(J)
        FUN1=((2.**RJ)*(Q**RJ))*((FACT(I)*FACT(2*I-J))) 
        FUN2=(FACT(J)*FACT(2*I)*FACT(I-J))             
 100    FUN=FUN+(FUN1/FUN2)                           
      RETURN                                       
      END                                         

***********************************************************************

      DOUBLE PRECISION FUNCTION FACT (NIM)                        

C**********************************************************************
C...
C...     This function computes the factorial of a number. If the
C...     number (NIM) is larger that 34, the value of the factorial
C...     becomes larger than the largest number a 32 bit machine can
C...     handle. The use of the double precision variable does away
C...     with this problem, which would arise if the conversion in
C...     DO 190 loop inside the main program, were larger than 13,
C...     which it is at many values of rij.
C...
C**********************************************************************

      IMPLICIT DOUBLE PRECISION (A-G, O-Z), INTEGER (I-N)
      IF(NIM) 10,20,30
 10   FACT=0.0                              
      RETURN                               
 20   FACT=1.0                            
      RETURN                             
 30   FACT=1.0                          
      DO 100 I=1,NIM                     
      FA=DFLOAT(I)                     
 100  FACT=FACT*FA                   
      RETURN                        
      END                          
