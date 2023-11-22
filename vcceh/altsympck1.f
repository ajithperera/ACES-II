C      
      SUBROUTINE ALTSYMPCK1(WIN, WOUT, NSIZIN, NSIZOT, ISCR,
     &                      IRREPX, SPTYPE)
C
C This routine accepts a symmetry packed four-index list and returns
C the same list but with an alternative scheme of symmetry packing.
C
C The list (Ab,Ij) is presumed to be packed (Ab,Ij). This routine returns
C the list packed in two possible ways.
C
C INPUT: 
C
C  WIN    - The symmetry packed (Ab,Ij) list.
C  NSIZIN - The total size of the sym. packed input vector.
C  NSIZOT - The total size of the sym. packed output vector.
C  SPTYPE - The spin type for the input list.
C           'AABB' for (AI,bj) returned.
C           'BBAA' for (bj,AI) returned.
C
C OUTPUT: 
C
C  WOUT   - The symmetry packed (AI,bj) or (bj,AI) list.
C       
C SCRATCH:
C
C  ISCR   - Scratch area to hold the symmetry vectors and inverse
C           symmetry vectors which are needed.
C
C           Size: (NVRTA*NVRTB + NOCCA*NOCCB + NVRTA*NOCCB +
C                  NVRTA*NOCCA + NVRTB*NOCCB)
C         
      IMPLICIT INTEGER (A-Z)
      CHARACTER*4 REORTP,SPTYPE
      DOUBLE PRECISION WIN(NSIZIN),WOUT(NSIZOT),ISCR(*)
      DIMENSION IOFFTAR(8)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
C
      REORTP = '1324'
C
      CALL SSTGEN(WIN, WOUT, NSIZIN, VRT(1,1), VRT(1,2), POP(1,1),
     &            POP(1,2), ISCR, IRREPX, REORTP)
C
      IF (SPTYPE .EQ. 'BBAA') THEN
C
         IOFFTAR(1) = 1
C
         DO 5 IRREPR = 1, NIRREP - 1
C
            IRREPL = DIRPRD(IRREPX, IRREPR)
C     
            DISSIZ = IRPDPD(IRREPL, 10)
            NUMDIS = IRPDPD(IRREPR, 9)
C
            IOFFTAR(IRREPR + 1) = IOFFTAR(IRREPR) + DISSIZ*NUMDIS
C
 5       CONTINUE
C
         IOFF1 = 1
C
         DO 10 IRREPR = 1, NIRREP
C
            IRREPL = DIRPRD(IRREPX, IRREPR)
C
            NUMDIS = IRPDPD(IRREPR, 10)
            DISSIZ = IRPDPD(IRREPL, 9)
C 
            CALL TRANSP(WOUT(IOFF1), WIN(IOFFTAR(IRREPL)), NUMDIS,
     &                  DISSIZ)
C
            IOFF1=IOFF1+NUMDIS*DISSIZ
C
 10      CONTINUE
C
         TOTSIZ = IOFF1 - 1
C
         CALL SCOPY(TOTSIZ, WIN, 1, WOUT, 1)
C
      ENDIF
C
      RETURN
      END
