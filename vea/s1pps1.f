C     
      SUBROUTINE S1PPS1(ICORE, MAXCOR, ISIDE, ISPIN)
C
C THE S1 CONTRIBUTION TO S1 INVOLVING U(PP) IS CALCULATED
C
C  Z[A,P] = SUM U(A,C) S(C,P)  [ISIDE = 1]
C            C
C
C  Z[C,P] = SUM S(A,P) * U(A,C) [ISIDE = 2]
C            A
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE
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
      DATA LISTU0 /92/
C
      LISTS1IN = LS1IN
      LISTS1EX = LS1OUT
      I000 = 1
      I010 = I000 + NFEA(ISPIN) * IINTFP
      I020 = I010 + VRT(SIRREP,ISPIN) * NS(SIRREP) * IINTFP
      I030 = I020 + VRT(SIRREP,ISPIN) * NS(SIRREP) * IINTFP
      CALL GETLST(ICORE(I000), 1, 1, 1, ISPIN,LISTU0)
      CALL GETLST(ICORE(I010), 1, 1, 1, ISPIN, LISTS1IN)
      CALL GETLST(ICORE(I020), 1, 1, 1, ISPIN, LISTS1EX)
C
C  CALCULATE OFFSET IN U
C
      ICOUNT = I000
      DO 10 AIRREP = 1, SIRREP-1
         ICOUNT = ICOUNT + VRT(AIRREP,ISPIN) *
     $      VRT(AIRREP,ISPIN) *IINTFP
 10   CONTINUE
C 
C NOW DO MULTIPLICATION
C
      NROW = VRT(SIRREP,ISPIN)
      NCOL = NS(SIRREP)
      NSUM = VRT(SIRREP,ISPIN)
      IF (ISIDE.EQ.1) THEN
         CALL XGEMM('N', 'N', NROW, NCOL, NSUM,ONE, ICORE(ICOUNT),
     $      NROW, ICORE(I010), NSUM, ONE, ICORE(I020), NROW)
      ELSE
         CALL XGEMM('T', 'N', NROW, NCOL, NSUM,ONE, ICORE(ICOUNT),
     $      NROW, ICORE(I010), NSUM, ONE, ICORE(I020), NROW)
      ENDIF
C
      CALL PUTLST(ICORE(I020), 1,1,1, ISPIN,LISTS1EX)
C
      RETURN
      END