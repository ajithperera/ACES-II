C     
      SUBROUTINE S1HHS1(ICORE, MAXCOR, ISIDE, ISPIN)
C
C THE S1 CONTRIBUTION TO S1 INVOLVING U(HH) IS CALCULATED
C
C  Z[I,P] = -SUM Hbar(m,I) S(m,P)  [ISIDE = 1]
C            J
C
C  Z[J,P] = -SUM S(P,m) * Hbar(j,m) [ISIDE = 2]
C            I
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE, ONEM
C
      DIMENSION ICORE(MAXCOR)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SINFO/ NS(8), SIRREP
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
C
      DATA ONE /1.0D0/
      DATA ONEM /-1.0D0/
      DATA LISTU0 /91/
C
CSSS      Write(6,"(a)") "-----------------S1HHS1-------------------"
      LISTS1IN = LS1IN
      LISTS1EX = LS1OUT
      I000 = 1
      I010 = I000 + NFMI(ISPIN) * IINTFP
      I020 = I010 + POP(SIRREP,ISPIN) * NS(SIRREP) * IINTFP
      I030 = I020 + POP(SIRREP,ISPIN) * NS(SIRREP) * IINTFP
      CALL GETLST(ICORE(I000), 1, 1, 1, ISPIN,LISTU0)
      CALL GETLST(ICORE(I010), 1, 1, 1, ISPIN, LISTS1IN)
#ifdef _DEBUG_LVL0
      ndim= POP(SIRREP,ISPIN)*NS(SIRREP)*IINTFP
      call checksum("s1(h)-i :",icore(i010),Ndim)
      call checksum("fmi     :",icore(i000),Ndim)
#endif 
      CALL GETLST(ICORE(I020), 1, 1, 1, ISPIN, LISTS1EX)
C
C  CALCULATE OFFSET IN U IRREPM=SIRREP=IRREPI)
C
      ICOUNT = I000
      DO 10 AIRREP = 1, SIRREP-1
         ICOUNT = ICOUNT + POP(AIRREP,ISPIN) *
     $                     POP(AIRREP,ISPIN) *IINTFP
 10   CONTINUE
C 
C NOW DO MULTIPLICATION
C
      NROW = POP(SIRREP,ISPIN)
      NCOL = NS(SIRREP)
      NSUM = POP(SIRREP,ISPIN)
      IF (ISIDE.EQ.1) THEN
         CALL XGEMM('T','N',NROW,NCOL,NSUM,ONEM,ICORE(ICOUNT),
     $      NROW,ICORE(I010),NSUM,ONE,ICORE(I020),NROW)
      ELSE
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONEM,ICORE(ICOUNT),
     $      NROW,ICORE(I010),NSUM,ONE,ICORE(I020),NROW)
      ENDIF
C
      CALL PUTLST(ICORE(I020), 1,1,1, ISPIN,LISTS1EX)
#ifdef _DEBUG_LVL0
      ndim= POP(SIRREP,ISPIN)*NS(SIRREP)*IINTFP
      call checksum("s1(h)-o :",icore(i020),Ndim)
#endif 
C
      RETURN
      END
