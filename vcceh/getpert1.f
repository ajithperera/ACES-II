C
      SUBROUTINE GETPERT1(ICORE, MAXCOR, MXCOR, IUHF, IRREPX, 
     &                    TLIST, IOFFSET, IOFFT1)
C
C  This routine loads the T1AA And T1BB Vectors into the very top of
C  Core, and returns a pointer array giving the offsets where each
C  irrep of the two spin cases begins.  For RHF, the addresses of the
C  AA and BB T1 vectors are identical and only one is held.
C
C  ICORE   - The core vector (T1 returned at top)
C  MAXCOR  - The total core size 
C  MXCOR   - The size of core below the start of the T1 vectors.
C  IUHF    - The UHF/RHF flag
C  IOFFT1  - A two dimensional array giving the address of
C            the beginning of each irrep in the T1 vector.
C            For example, IOFFT1(3, 2) gives the address of
C            the first element of the third irrep for spin
C            case 2 (Bb).
C
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR),IOFFT1(8,2),ITSTART(2)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,NF1BB,NF2BB
C
C Common blocks used in the quadratic term
C 
      COMMON /QADINTI/ INTIMI, INTIAE, INTIME
      COMMON /QADINTG/ INGMNAA, INGMNBB, INGMNAB, INGABAA,
     &                 INGABBB, INGABAB, INGMCAA, INGMCBB,
     &                 INGMCABAB, INGMCBABA, INGMCBAAB,
     &                 INGMCABBA
      COMMON /APERTT1/ IAPRT1AA             
      COMMON /BPERTT1/ IBPRT1AA
      COMMON /APERTT2/ IAPRT2AA1, IAPRT2BB1, IAPRT2AB1, IAPRT2AB2, 
     &                 IAPRT2AB3, IAPRT2AB4, IAPRT2AA2, IAPRT2BB2,
     &                 IAPRT2AB5 
      MXCOR = MAXCOR
C
C  Compute offsets for beginning of T1aa and T1bb (1000 loop) and offsets for
C  beginning of irreps (2000 loop).
C
      DO 1000 ISPIN = 2, 2 - IUHF, -1
C
         NLIST = 22 + ISPIN
         TLIST2 = ISPIN
C
         IF (IUHF .EQ. 0) THEN
            NLIST = 23
            TLIST2 = 1 
         ENDIF
C
         T1SIZ = IRPDPD(IRREPX, ISYTYP(1, NLIST))
         ITSTART(ISPIN) = MXCOR - T1SIZ*IINTFP + 1
C      
         IF(IUHF .EQ. 0) ITSTART(1) = ITSTART(2)
         MXCOR = MXCOR - T1SIZ*IINTFP
         IOFF = ITSTART(ISPIN)
C     
         DO 2000 IRREPR = 1, NIRREP
            IRREPL = DIRPRD(IRREPR, IRREPX)
            IOFFT1(IRREPL, ISPIN) = IOFF
            IF (IUHF .EQ. 0) IOFFT1(IRREPL, 1) = IOFF
            IOFF = IOFF + VRT(IRREPL, ISPIN)*POP(IRREPR, ISPIN)*IINTFP
 2000    CONTINUE
C     
         CALL GETLST(ICORE(ITSTART(ISPIN)),1, 1, 1, TLIST2, TLIST)
C     
 1000 CONTINUE
C      
      RETURN
      END 
