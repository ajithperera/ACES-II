      SUBROUTINE PRINTSUM(ROOT,FSTR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL MBPT2,CC
      DIMENSION TML(3),TMR(3),DIPSTR(3),OSCSTR(3),POL(3,3)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/POLAR/POLTOT(3,3)
      COMMON/REFTYPE/MBPT2,CC
      DATA ONE /1.0D0/
      DATA FACT0,FACT1,FACT2 /0.6666666666666666D0,27.2113957D0,
     &   8065.54093D0/
C
      ITHR=3*IINTFP
      CALL GETREC(20,'JOBARC','TMRIGHT ',ITHR,TMR)
      CALL GETREC(20,'JOBARC','TMLEFT  ',ITHR,TML)
      CALL GETREC(20,'JOBARC','TOTENERG',IINTFP,E)
C
      CALL ZERO(POL,9)
      DO 10 IXYZ=1,3
       DIPSTR(IXYZ)=TMR(IXYZ)*TML(IXYZ)
       OSCSTR(IXYZ)=FACT0*ROOT*DIPSTR(IXYZ)    
10    CONTINUE
      FSTR=SNRM2(3,OSCSTR,1)
      IF(ABS(ROOT).GT.1.D-5)THEN
       DO 15 I=1,3
        DO 16 J=1,3
         POL(I,J)=POL(I,J)+2.0D0*(TML(I)*TMR(J))/ROOT
16      CONTINUE
15     CONTINUE
      ENDIF 
      CALL SAXPY(9,ONE,POL,1,POLTOT,1)
      WAVLEN=1.0D7/(ROOT*FACT1*FACT2)
      WAVNUM=ROOT*FACT1*FACT2
      WRITE(6,1007)
      WRITE(6,1006)
      WRITE(6,1007)
      WRITE(6,1002)(TMR(I),I=1,3)
      WRITE(6,1003)(TML(I),I=1,3)
      WRITE(6,1004)(DIPSTR(I),I=1,3)
      WRITE(6,1005)(OSCSTR(I),I=1,3)
      WRITE(6,1007)

      WRITE(6,1000)ROOT*FACT1,WAVLEN,WAVNUM
      IF(CC)THEN
       WRITE(6,1008)ROOT+E
      ELSE
       WRITE(6,1018)ROOT+E
      ENDIF
      WRITE(6,1001)FSTR
      WRITE(6,1009)POL(1,1),POL(2,2),POL(3,3),POL(1,2),POL(1,3),
     &             POL(2,3)
      WRITE(6,1010)POLTOT(1,1),POLTOT(2,2),POLTOT(3,3),
     &             POLTOT(1,2),POLTOT(1,3),POLTOT(2,3)
C
1009  FORMAT(T3,'Contribution to electric polarizability: ',/,
     &       T5,'XX: ',F15.9,' YY: ',F15.9,' ZZ: ',F15.9,/,
     &       T5,'XY: ',F15.9,' XZ: ',F15.9,' YZ: ',F15.9)
1010  FORMAT(T3,'Cumulative electric polarizability: ',/,
     &       T5,'XX: ',F15.9,' YY: ',F15.9,' ZZ: ',F15.9,/,
     &       T5,'XY: ',F15.9,' XZ: ',F15.9,' YZ: ',F15.9)
1000  FORMAT(T3,'Transition energy ',F10.4,' eV  (',F9.4,' nm; ',F11.3,
     &          ' cm-1)')
1008  FORMAT(T3,'Total EOM-CCSD electronic energy ',F20.12,' a.u.')
1018  FORMAT(T3,'Total EOM-MBPT(2) electronic energy ',F20.12,' a.u.')
1001  FORMAT(T3,'Norm of oscillator strength : ',F12.8) 
1002  FORMAT(T3,'Right Transition Moment',T28,F12.8,T42,F12.8,T56,F12.8)
1003  FORMAT(T3,'Left Transition Moment ',T28,F12.8,T42,F12.8,T56,F12.8)
1004  FORMAT(T3,'Dipole Strength        ',T28,F12.8,T42,F12.8,T56,F12.8)
1005  FORMAT(T3,'Oscillator Strength    ',T28,F12.8,T42,F12.8,T56,F12.8)
1006  FORMAT(T34,'X',T48,'Y',T62,'Z')
1007  FORMAT(70('-'))
C
      RETURN
      END
