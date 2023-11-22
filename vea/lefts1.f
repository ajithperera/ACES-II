      SUBROUTINE LEFTS1(ICORE,MAXCOR,IUHF, ISPIN)
C
C   THIS SUBROUTINE CALCULATES THE CI-COEFFICIENTS OF THE
C   LEFT HAND EIGENVECTOR DERIVING FROM S2 * EXP [- T1].
C
C CONTRACTION :
C
C     Z(b,p) = S1(b,p) +  SUM [T1(k,c) (2 S(cb, kp) - S(bc, kp)) [ RHF CASE ]
C              k,c
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      LOGICAL SKIP, print
C
      DIMENSION ICORE(MAXCOR), NUMSZS(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SINFO/ NS(8), SIRREP
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
C
      DATA ONE  /1.0/
      DATA ZILCH/0.0/
      DATA ONEM/-1.0/
C
      print = .false.
      CALL GETLST(ICORE, 1, 1, 1, ISPIN, LS1IN)
      CALL PUTLST(ICORE, 1, 1, 1, ISPIN, LS1OUT)
       if (print) then
          write(6,*) ' lefts1, s1 coefficients '
          call output(icore,1,vrt(sirrep,ispin),1,1,
     $       vrt(sirrep,ispin),1,1)
       endif
C
      DO 5 XIRREP = 1, NIRREP
         MIRREP = DIRPRD(XIRREP, SIRREP)
         NUMSZS(XIRREP) = POP(MIRREP, ISPIN) * NS(SIRREP)
 5    CONTINUE      
C
      DO 25 SINSPIN = 1, 1+IUHF
         MSPIN = SINSPIN
         LISTS2IN = LS2IN(ISPIN, SINSPIN + 1 - IUHF)
       CALL GETLEN(LENS, VRT(1,SINSPIN), VRT(1,ISPIN),
     $      POP(1,SINSPIN), NS)
       I000=1
       I010=I000+LENS*IINTFP
       IF ((IUHF.NE.0) .AND. (SINSPIN.EQ.ISPIN)) THEN
          CALL GETEXP2(ICORE(I000), LENS, NUMSZS, IRPDPD(1, 18+ISPIN),
     $       LISTS2IN, VRT(1, ISPIN), IRPDPD(1, ISYTYP(1, LISTS2IN)))
       ELSE
          CALL GETALLS2(ICORE(I000), LENS,POP(1,SINSPIN),NS(1),1,
     $       LISTS2IN)
       ENDIF
C
C SPINADAPT THE S COEFFICIENTS BY IRREP
C
       IF (IUHF.EQ.0) THEN
          skip = .false.
          if (skip) then
             write(6,*)' no spinadaptation in s2uhps1 '
          else
          ICOUNT = I000
          DO 10 XIRREP = 1,NIRREP
             MIRREP = DIRPRD(XIRREP,SIRREP)
             DISSYS=IRPDPD(XIRREP,ISYTYP(1,LISTS2IN))
             NUMDSS=POP(MIRREP,1) * NS(SIRREP)
             I020 = I010 + NUMDSS * IINTFP
             I030 = I020 + NUMDSS * IINTFP
             CALL SPINAD3(XIRREP, VRT(1,1), DISSYS, NUMDSS,
     $          ICORE(ICOUNT), ICORE(I010), ICORE(I020))
             ICOUNT = ICOUNT + NUMDSS * DISSYS * IINTFP
 10       CONTINUE
       endif
       ENDIF
C
C REORDER THE S COEFFICIENT ACCORDING TO S(AB, MP) -> S(AM, BP)
C
       I020 = I010 + LENS*IINTFP
        IF(I020.GT.MAXCOR)THEN
C
C OUT-OF-CORE ALGORITHM
C
         WRITE(6,*)' out-of-core AB not coded, s2uhps1 '
         WRITE(6,1000) I020, MAXCOR
 1000    FORMAT('  REQUIRED MEMORY :', I12, '  AVAILABLE : ', I12)
         call errex
C
        ENDIF
C
       CALL SCOPY(LENS,ICORE(I000),1,ICORE(I010),1)
       CALL SSTGEN(ICORE(I010), ICORE(I000), LENS, VRT(1,MSPIN),
     $    VRT(1,ISPIN), POP(1,MSPIN), NS, ICORE(I020), 1, '1324')
C
C ICORE(I000) CONTAINS ALL [SPIN ADAPTED] S COEFFICIENTS [2 S(AB, MP) -
C       S(BA, MP)], ORDERED AS AM;BP
C
       I020 = I010 + NT(MSPIN) * IINTFP
       CALL GETLST(ICORE(I010), 1,1,1,MSPIN, 90)
       BSPIN = ISPIN
C
C ONLY THE FIRST IRREP OF S VECTOR IS NEEDED
C
       DISSYS= IRPDPD(1, 15 + SINSPIN)
       NUMDSS=VRT(SIRREP,ISPIN) * NS(SIRREP)
C
       CALL GETLST(ICORE(I020), 1, 1, 1, ISPIN, LS1OUT)
       if (print) then
          write(6,*) ' lefts1, s1 coefficients '
          call output(icore(I020),1,vrt(sirrep,ispin),1,1,
     $       vrt(sirrep,ispin),1,1)
       endif
C
C PERFORM MATRIX MULTIPLICATION
C
C         Z(1,BP) = SUM U(CK,1) * S(CK, BP) 
C                      CK
C  
          NROW=1
          NCOL=NUMDSS
          NSUM=DISSYS
            CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(I010),NSUM,
     &               ICORE(I000),NSUM,ONE,ICORE(I020),NROW)
C
            CALL PUTLST(ICORE(I020),1, 1, 1, ISPIN, LS1OUT)
       if (print) then
          write(6,*) ' lefts1, s1 coefficients '
          call output(icore(I020),1,vrt(sirrep,ispin),1,1,
     $       vrt(sirrep,ispin),1,1)
       endif
C
 25   CONTINUE
C
      RETURN
      END
