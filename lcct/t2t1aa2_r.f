C
C     THIS SUBROUTINE CALCULATES THE FOLLOWING CONTRIBUTIOBN TO
C     THE T1 INCREMENT (LAMBDA = .FALSE.)
C
C      Z(I,A) = - SUM M<N,E T(MN,AE) <MN//IE>
C
C OR   Z(i,a) = - SUM m<n,e T(mn,ea) <mn//ej)
C
C     OR THE L1 INCREMENT (LAMBDA = .FALSE.)
C
C      Z(I,A) = - SUM M<N,E L(MN,AE) <MN//IE>
C
C OR   Z(i,a) = - SUM m<n,e L(mn,ea) <mn//ej)
C
      SUBROUTINE T2T1AA2_R(W,ICORE,MAXCOR,POP,VRT,ISPIN,LISTT2,LSTOFF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT,DISSYT,DISSYW
      DIMENSION ICORE(MAXCOR),POP(8),VRT(8)
      DIMENSION W(1)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /FLAGS/ IFLAGS(100)
C
      DATA AZERO,ONE,ONEM /0.0D0,1.0D0,-1.0D0/
C
C  The iflags(39) flag, reserved for "standard" or "semicanonical"
C  orbitals, is used here to determine which list is needed for the
C  ROHF-MBPT code, since this routine is called by two others at
C  different times.  The value of "99" is set before the second call
C  to E3S in ROHFPT.
C
      LISTT = (LISTT2 - 1) + ISPIN
CSSS      LISTT = LISTT2  + ISPIN

      LISTW=ISPIN+6+LSTOFF
C
C    LOOP OVER IRREPS
C
      DO 1000 IRREP=1,NIRREP
       NVRTSQ=0
       DO 1001 IRREPJ=1,NIRREP
        NVRTSQ=NVRTSQ+VRT(IRREPJ)*VRT(DIRPRD(IRREP,IRREPJ))
1001   CONTINUE
       DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
       DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
       NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
       NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
       I001=1
       I002=I001+IINTFP*NUMSYT*NVRTSQ
       I003=I002+IINTFP*MAX(NUMSYW*DISSYW,NUMSYT*DISSYT)
       I004=I003+IINTFP*MAX(NUMSYT,NUMSYW,DISSYT,DISSYW)
       I005=I004+IINTFP*MAX(NUMSYT,NUMSYW,DISSYT,DISSYW)
       I006=I005+IINTFP*MAX(NUMSYT,NUMSYW,DISSYT,DISSYW)
       IF(I006.LT.MAXCOR) THEN
        CALL GETLST(ICORE(I002),1,NUMSYT,1,IRREP,LISTT)
        CALL TRANSP(ICORE(I002),ICORE(I001),NUMSYT,DISSYT)
        CALL SYMEXP(IRREP,VRT,NUMSYT,ICORE(I001))
        CALL GETLST(ICORE(I002),1,NUMSYW,2,IRREP,LISTW)
        CALL SYMTR1(IRREP,POP,VRT,DISSYW,ICORE(I002),ICORE(I003),
     &              ICORE(I004),ICORE(I005))
        IOFFT=0
        IOFFW=0
        IOFFS=1
        DO 100 IRREPJ=1,NIRREP
         NVRTJ=VRT(IRREPJ)
         NOCCJ=POP(IRREPJ)
         IRREPI=DIRPRD(IRREPJ,IRREP)
         NVRTI=VRT(IRREPI)
C
         CALL XGEMM('T','N',NVRTJ,NOCCJ,NVRTI*DISSYW,ONE,
     &              ICORE(I001+IOFFT*IINTFP),NUMSYT*NVRTI,
     &              ICORE(I002+IINTFP*IOFFW),DISSYW*NVRTI,ONE,
     &              W(IOFFS),NVRTJ)
        IOFFT=IOFFT+NVRTI*NVRTJ*NUMSYT
        IOFFW=IOFFW+NVRTI*NOCCJ*DISSYW
        IOFFS=IOFFS+NVRTJ*NOCCJ
100    CONTINUE
      ELSE
       STOP 'W5AA'
      ENDIF
1000  CONTINUE
      RETURN
      END
