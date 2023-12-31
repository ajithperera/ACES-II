C
C THIS ROUTINE WRITES OUT A SUMMARY OF THE MOLECULAR ORBITALS
C USED IN THE CORRELATED PART OF THE CALCULATION.
C
      SUBROUTINE ORBSUM(SCR,ISCR,NORBS,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POP,VRT,DIRPRD
      CHARACTER*5 SPCASE(2)
      CHARACTER*1 SP(2)
      DIMENSION SCR(NORBS),ISCR(NORBS)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
C
      DATA SPCASE/'alpha','beta '/
      DATA SP    /'A','B'/      
C
      CALL ACES_COM_SYM
C
      NCOL=NORBS/2 + MOD(NORBS,2)
C
      WRITE(6,600)
C
C LOOP OVER SPIN CASES
C
      DO 10 ISPIN=1,1+IUHF
       IOFF=1
       DO 11 IRREP=1,NIRREP
        N=POP(IRREP,ISPIN)
        DO 12 I=1,N
         ISCR(IOFF)=IRREP
         IOFF=IOFF+1
12      CONTINUE
11     CONTINUE
       DO 21 IRREP=1,NIRREP
        N=VRT(IRREP,ISPIN)
        DO 22 I=1,N
         ISCR(IOFF)=IRREP
         IOFF=IOFF+1
22      CONTINUE
21     CONTINUE
C       
       WRITE(6,1000)
       IF(IUHF.NE.0)WRITE(6,500)SPCASE(ISPIN)(1:LINBLNK(SPCASE(ISPIN)))
       WRITE(6,*)
       WRITE(6,2000)
       WRITE(6,1000)
       CALL GETREC(20,'JOBARC','SCFEVAL'//SP(ISPIN),NORBS*IINTFP,SCR)
C
       IFIRST=1
       DO 20 IFIRST=1,NCOL
        IOTHER=IFIRST+NCOL
        IF(IOTHER.LE.NORBS)THEN
         WRITE(6,3000)IFIRST,SCR(IFIRST),ISCR(IFIRST),
     &                IOTHER,SCR(IOTHER),ISCR(IOTHER)
        ELSE
         WRITE(6,3001)IFIRST,SCR(IFIRST),ISCR(IFIRST)
        ENDIF
20     CONTINUE
       WRITE(6,1000)
10    CONTINUE
C
3000  FORMAT(T3,I4,T10,F15.7,T32,I2,T40,I4,T47,F15.7,T69,I2)
3001  FORMAT(T3,I4,T10,F15.7,T32,I2)
1000  FORMAT(72('-'))
500   FORMAT(T26,'* Spin case ',A,' *')
600   FORMAT(T3,'Summary of active molecular orbitals: ')
2000  FORMAT(T3,'Index',T14,'Eigenvalue',T28,'Symmetry',
     &       T40,'Index',T51,'Eigenvalue',T65,'Symmetry')
C
      RETURN
      END
