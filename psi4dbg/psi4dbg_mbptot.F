
      SUBROUTINE PSI4DBG_MBPTOUT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FLAGS/ IFLAGS(100)
cAP - 500 is the maximum number of CC iterations (not basis functions)
      COMMON /ENERGY/ ENERGY(500,2),IXTRLE(500)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      CHARACTER*13 LABEL2(22)
      CHARACTER*12 LABEL3(6,2)
      CHARACTER*9 LABEL1(6)
      DATA LABEL1 /'D-MBPT(2)','D-MBPT(3)','D-MBPT(4)',
     &             'Q-MBPT(4)','S-MBPT(4)','T-MBPT(4)'/
C
      DATA LABEL2 /'MBPT(2)      ','MBPT(3)      ',
     &             'SDQ-MBPT(4)  ','MBPT(4)      ',
     &             'LCCD         ','LCCSD        ',
     &             'LCCSDT       ','CCD          ',
     &             'CCSD         ','CCSD+T(CCSD) ',
     &             'CCSD+T*(CCSD)','CCSDT-1      ',
     &             'CCSDT-1b     ','CCSDT-2      ',
     &             'CCSDT-3      ','CCSDT-4      ',
     &             'CCSDT        ','UCCSD(4)     ',
     &             'UCCSDT(4)    ','QCISD(T)     ',
     &             'CCSD(T)      ','CCD+ST(CCD)  '/
C
      DATA LABEL3 /'   S-MBPT(2)','   S-MBPT(3)',
     &             '   L-MBPT(4)','  NL-MBPT(4)',
     &             ' FT1-MBPT(4)','   T-MBPT(4)',
     &             '   D-MBPT(2)','   D-MBPT(3)',
     &             'WT12-MBPT(4)','WT12-MBPT(4)',
     &             ' WT2-MBPT(4)','            '/
C
      CALL GETREC(20,'JOBARC','SCFENEG ',IINTFP,Z)
      WRITE(*,200)
      WRITE(*,10)
      WRITE(*,200)
C
CJDW 3/98. Change for non-HF.
C
      IF(IFLAGS(11).LT.2.AND.IFLAGS(22).EQ.0.AND.IFLAGS(38).EQ.0) THEN
C
       DO 15 I=1,6
        IF(ENERGY(I,1).NE.0.0)THEN
         Z=Z+ENERGY(I,1)
         WRITE(*,20)LABEL1(I),ENERGY(I,1),Z
        ENDIF
15     CONTINUE
C
CJDW 3/98. Change for non-HF.
C
      ELSEIF((IFLAGS(11).EQ.2.OR.IFLAGS(38).NE.0).AND.
     &                           IFLAGS(22).EQ.0) THEN
C
       DO 30 I=1,6
        DO 31 J=1,2
         IF(ENERGY(I,J).NE.0.0) THEN
          Z=Z+ENERGY(I,J)
          WRITE(*,22)LABEL3(I,J),ENERGY(I,J),Z
         ENDIF
31      CONTINUE
30     CONTINUE
C
      ELSE IF(IFLAGS(22).NE.0) THEN
C
       DO 40 I=1,6
        IF(ENERGY(I,2).NE.0.0) THEN
         Z=Z+ENERGY(I,2)
         WRITE(*,22)LABEL3(I,2),ENERGY(I,2),Z
        ENDIF
40     CONTINUE
C
      ENDIF
      WRITE(*,200)
      WRITE(*,100)LABEL2(1),Z
C
CJDW 3/27/97. Finally, write the energy to the JOBARC file.
C
      CALL PUTREC(1,'JOBARC','TOTENERG',IINTFP,Z)
cmn again for vibronic purposes:
      CALL PUTREC(1,'JOBARC','TOTENER2',IINTFP,Z)
      CALL PUTREC(1,'JOBARC','PARENERG',IINTFP,Z)
C
200   FORMAT(T3,55('-'))
100   FORMAT(T9,'Total ',A,' energy: ',F16.10,' a.u.')
10    FORMAT(T8,'Correction',T26,'Increment',T45,'Cumulative')
20    FORMAT(T9,A9,T22,F16.10,T41,F16.10)
22    FORMAT(T6,A12,T22,F16.10,T41,F16.10)
      RETURN
      END
