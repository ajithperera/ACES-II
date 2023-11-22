      SUBROUTINE CONSTANTS
C     
C     This routine sets the constants and fractions needed for
C     the B-LYP functional.
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      LOGICAL YESNO
C
      CHARACTER*5 PARAM2
      CHARACTER*6 PARAM1
      CHARACTER*80 FNAME
C
      COMMON /PAR/ PI,PIX4
      COMMON /CNST/ ALPHA,BETA,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,
     $ X14,X15
      COMMON /FRAC/ F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12,F13,F14,F15
C     
C     Check for file containing user values for functionals
      CALL GFNAME('CONSTANT',FNAME,ILENGTH)
      INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=YESNO)
      IF(YESNO)THEN
      OPEN(UNIT=11,FILE=FNAME(1:ILENGTH),STATUS='OLD',
     &        FORM='FORMATTED',ACCESS='SEQUENTIAL')
C
C     The value of alpha for X-alpha (1.d+00 yields the standard X-alpha
C     value, 2/3 yields the standard LDA value)
C
         READ(11,1010)PARAM1,XX
C
 1010 FORMAT(A6,F7.5)
C
         IF(XX.EQ.0.D+00)THEN
            ALPHA=2.D+00/3.D+00
         ELSE
            ALPHA=XX
         ENDIF
C
C     The value which scales the size of the correlation hole in
C     the LYP functional (BETA=1.D+00 yields the value determined by
C     Colle and Salvetti)
C
         READ(11,1011)PARAM2,XX
         IF(XX.EQ.0.D+00)THEN
            BETA=1.D+00
         ELSE
            BETA=XX
         ENDIF
C
 1011 FORMAT(A5,F8.5)
C
      ENDIF
      IF(ALPHA.EQ.0.D+00)ALPHA=2.D+00/3.D+00
      IF(BETA.EQ.0.D+00)BETA=1.D+00
C
C     The Becke XC  functional parameter (0.0042d+00 is the value 
C     determined by Becke)
      X3=0.0042D+00
C
C     Lee-Yang-Parr functional parameters
C     A
      X4=0.04918D+00
C     B
      X5=(1.D+00/(BETA*BETA))*0.132D+00
C     C
      X6=(1.D+00/BETA)*0.2533D+00
C     D
      X7=(1.D+00/BETA)*0.349D+00
C     A*B
      X8=X4*X5
C
C     Constant for X-alpha
      X9=ALPHA*(-2.25D+00*(3.D+00/(4.D+00*PI))**(1.D+00/3.D+00))
C
C     Constant for term 2 in the LYP functional
      X10=(2.D+00**(11.D+00/3.D+00))*0.3D+00*
     $   (3.D+00*PI*PI)**(2.D+00/3.D+00)
C
C     6*B for denominator in Becke term
      X11=6.D+00*X3
C
C     4*A for term 1 in LYP
      X12=4.D+00*X4
C
C     Constant for X-alpha derivative
      X13=ALPHA*(-2.D+00*(3.D+00/(4.D+00*PI))**(1.D+00/3.D+00))
      X14=6.D+00*X3*X3
      X15=2.D+00*X3
C
C     Values of some fractions to be used multiple times
      F1=4.D+00/3.D+00
      F2=1.D+00/3.D+00
      F3=8.D+00/3.D+00
      F4=1.D+00/9.D+00
      F5=11.D+00/3.D+00
      F6=0.D+00
      F7=0.D+00
      F8=0.D+00
      F9=0.D+00
      F10=0.D+00
      F11=0.D+00
      F12=0.D+00
      F13=0.D+00
      F14=0.D+00
      F15=0.D+00
C     
      RETURN
      END
