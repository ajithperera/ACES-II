C
      SUBROUTINE ALTSYMPCK2(WIN, WOUT, NSIZIN, NSIZOT, ISCR, IRREPX, 
     &                      SPTYPE, RERDTP)
C
C This routine accepts a symmetry packed four-index list returns the
C same list but with an alternative scheme of symmetry packing.
C
C The list (A<B,I<J) is presumed to be packed (AB,IJ). This routine returns
C the list (AI,BJ) packed as (AI,BJ) or (AJ,BI). (These possibilities
C differ by a sign).
C
C IMPORTANT:
C
C  The space allocated for both WIN and WOUT in the
C  calling routine must be equal to the size of WOUT.
C
C INPUT: 
C
C  WIN    - The symmetry packed (AB,IJ) list.
C  NSIZIN - The total size of the sym. packed input vector.
C  NSIZOT - The toatl size of the sym. packed output vector.
C  SPTYPE - The spin type for the input list.
C           'AAAA' for (AB,IJ)
C           'BBBB' for (ab,ij)
C  RERDTP - 'AIBJ' for (AI,BJ) positive, 'AJBI' for (AJ,BI) positive.
C
C OUTPUT: 
C
C  WOUT   - The symmetry packed (AI,bj) or (bj,AI) list.
C       
C SCRATCH:
C
C  ISCR   - Scratch area to hold the symmetry vectors and inverse
C           symmetry vectors which are needed. 
C           (SIZE: NVRT*(NVRT-1)/2 + NOCC*(NOCC-1)/2 + NOCC*NVRT)
C         
      IMPLICIT INTEGER (A-Z)
      CHARACTER*4 REORTP,SPTYPE,RERDTP
      DOUBLE PRECISION WIN(NSIZIN),WOUT(NSIZOT),ISCR(*)
C
      COMMON/MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/ NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/FILES/ LUOUT,MOINTS
C
      IF (SPTYPE .EQ. 'AAAA') THEN
         ISPIN = 1
      ELSEIF(SPTYPE.EQ.'BBBB')THEN
         ISPIN = 2
      ELSE
         WRITE(LUOUT, 1000) SPTYPE
 1000    FORMAT(T3, 'ALTSYMPCK2, Unknown spin type ',A4,' passed in.')
         CALL ERREX
      ENDIF
C
      IF (RERDTP .EQ. 'AIBJ') THEN
         REORTP = '1324'
      ELSEIF (RERDTP .EQ. 'AJBI') THEN
         REORTP = '1423'
      ELSE
         WRITE(LUOUT, 1001) RERDTP
 1001    FORMAT(T3, 'ALTSYMPCK2, Unknown reorder type ',A4,' passed 
     &          in.')
         CALL ERREX
      ENDIF
C
C EXPAND LIST FROM A<B,I<J TO AB,IJ
C
      IOFFT = 1
      IOFFF = 1
C
      DO 10 IRREPR = 1, NIRREP
C
         IRREPL = DIRPRD(IRREPX, IRREPR)
C
         DISSZT = IRPDPD(IRREPL, ISPIN)
         NUMDST = IRPDPD(IRREPR, 2+ISPIN)
         DISSZF = IRPDPD(IRREPL, 18+ISPIN)
         NUMDSF = IRPDPD(IRREPR, 20+ISPIN)
C
         CALL SCOPY(NUMDST*DISSZT, WIN(IOFFT), 1, WOUT(IOFFF), 1)
         CALL SYMEXP(IRREPR, POP(1, ISPIN), DISSZT, WOUT(IOFFF)) 
         CALL SYMEXP2(IRREPL, VRT(1,ISPIN), DISSZF, DISSZT,
     &                NUMDSF, WOUT(IOFFF), WOUT(IOFFF))
C
         IOFFF = IOFFF + NUMDSF*DISSZF
         IOFFT = IOFFT + NUMDST*DISSZT
 10   CONTINUE
C
      ITOTSZF = IOFFF - 1
      CALL SCOPY(ITOTSZF, WOUT, 1, WIN, 1)
C     
      CALL SSTGEN(WIN, WOUT, ITOTSZF, VRT(1,ISPIN), VRT(1,ISPIN),
     &            POP(1,ISPIN), POP(1,ISPIN), ISCR, IRREPX, REORTP)
C     
      RETURN
      END
