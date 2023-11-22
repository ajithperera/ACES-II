C
      SUBROUTINE MKDBLRNG2(ICORE, MAXCOR, LISTG, LISTT, CALMOD,
     &                     MAXSIZ, LSTSCR, ADD41, IRREPX)
C
      IMPLICIT INTEGER (A-Z)
      LOGICAL ADD41
      DOUBLE PRECISION ONE,TWO,HALF
      CHARACTER*3 CALMOD
      DIMENSION ICORE(MAXCOR)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
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
C
      DATA ONE  /1.0D+00/
      DATA TWO  /2.0D+00/
      DATA HALF /0.5D+00/
C
C Loop over irreps.
C     
      DO 100 IRREPEM = 1, NIRREP
         IRREPGEM = IRREPEM
         IRREPGBJ = DIRPRD(IRREPGEM, IRREPX)
         IRREPTEM = IRREPGEM
         IRREPTAI = DIRPRD(IRREPTEM, IRREPX)
         IRREPQBJ = IRREPGBJ
         IRREPQAI = IRREPQBJ
C
         DISSYG = IRPDPD(IRREPGEM, ISYTYP(1, LISTG))
         NUMSYG = IRPDPD(IRREPGBJ, ISYTYP(2, LISTG))
         DISSYT = IRPDPD(IRREPTAI, ISYTYP(1, LISTT))
         NUMSYT = IRPDPD(IRREPTEM, ISYTYP(2, LISTT))
         DISSYQ = DISSYT
         NUMSYQ = NUMSYG
         MAXSIZ = MAX(NUMSYQ*DISSYQ, MAXSIZ)
C
         I000 = 1
         I010 = I000 + IINTFP*NUMSYQ*DISSYQ
         I020 = I010 + IINTFP*MAX(NUMSYT*DISSYT, NUMSYQ*DISSYQ)
         I030 = I020 + IINTFP*DISSYG*NUMSYG
C     
C Can we do it in core?
C     
         IF (I030 .LT. MAXCOR) THEN
C     
C Form the symmetry packed Q vector
C     
            CALL RNGCONTR(ICORE(I020), ICORE(I010), ICORE(I000),
     &                    DISSYG, NUMSYG, DISSYT, NUMSYT, DISSYQ,
     &                    NUMSYQ, LISTG, LISTT, IRREPGBJ, IRREPTEM,
     &                    CALMOD, .FALSE.)
         ELSE
            CALL INSMEM('MKDBLRNG2', I030, MAXCOR)
         ENDIF
C
         CALL PUTLST(ICORE(I000), 1, NUMSYQ, 1, IRREPQBJ, LSTSCR)
C
 100  CONTINUE
C     
C Now we have (aI,Bj) ordered piece. Change the ordering to (aj,BI)
C add that to the original (aI,Bj). Note that in SSTRNG call we 
C use spin-case 'AAAA' instead of 'ABAB', but for a closed shell
C this does not make any difference. (This routine calls only 
C for closed shell case).
C
      ISIZE = ISYMSZ(ISYTYP(1, LSTSCR), ISYTYP(2, LSTSCR))
C
      I000 = 1
      I010 = I000 + ISIZE*IINTFP
      I020 = I010 + ISIZE*IINTFP
C
      CALL GETALL(ICORE(I000), ISIZE, 1, LSTSCR)
      CALL SSTRNG(ICORE(I000), ICORE(I010), ISIZE, ISIZE, ICORE(I020),
     &            'AAAA')
      CALL SAXPY(ISIZE, HALF, ICORE(I010), 1, ICORE(I000), 1)
C     
C  Add in piece from T1*T1*Hbar calculated in T1T1INDBL1
C     
      IF(ADD41) THEN
         CALL GETALL(ICORE(I010), ISIZE, 1, 41)
         CALL SAXPY(ISIZE, ONE, ICORE(I010), 1, ICORE(I000), 1)
      ENDIF
C     
C Symmetrize the quantity
C     
      IOFF = I000
C
      DO 200 IRREP = 1, NIRREP
         NUMSYQ = IRPDPD(IRREP, ISYTYP(2, LISTG))
         CALL MPMT(ICORE(IOFF), NUMSYQ)
         IOFF = IOFF + NUMSYQ*NUMSYQ*IINTFP
 200  CONTINUE
C     
      CALL PUTALL(ICORE(I000), 1, 1, LSTSCR)
C     
      RETURN
      END
