      SUBROUTINE S2HPPHS2(ICORE,MAXCOR,IUHF, ISIDE, ISPIN)
C
C   THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF S2 TO
C   S2 INVOLVING THE HPPH BLOCK OF EFFECTIVE INTEGRALS W (RINGS).
C
C CONTRACTIONS :
C
C     Z(Ab,Ip) = SUM [<kA|cI> (2 S(cb, kp) - S(bc, kp)) [ RHF CASE ]
C                k,c
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      logical prints, print, skip
C
      DIMENSION ICORE(MAXCOR)
      DIMENSION LISTW0(2,2), NUMSZS(8)
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
C     56: KaCi, ordered as C,k ; a,i
C
C  PARAMETERS DETERMINING THE LISTS THAT ARE USED
C   SINSPIN = SPIN OF C,K,   SOUTSPIN: SPIN OF A AND I,  ISPIN: SPIN OF B
C   WE HAVE W(SINSPIN, SOUTSPIN)      
C
      prints = .false.
      LISTW0(1,1) = 54
      LISTW0(1,2) = 56
      LISTW0(2,1) = 57
      LISTW0(2,2) = 55
C
C LISTS 54 AND 55 CONTAIN   - W(MBEJ)  !!!
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
             write(6,*)' no spinadaptation in s2hpphs2'
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
       CALL SCOPY(LENS,ICORE(I000),1,ICORE(I010),1)
       CALL SSTGEN(ICORE(I010), ICORE(I000), LENS, VRT(1,MSPIN),
     $    VRT(1,ISPIN), POP(1,MSPIN), NS, ICORE(I020), 1, '1324')
C
C ICORE(I000) CONTAINS ALL [SPIN ADDAPTED] S COEFFICIENTS [2 S(AB, MP) -
C       S(BA, MP)], ORDERED AS AM;BP
C
C
C  NOW CARRY OUT THE SUMMATION WHICH USES THE [ SPIN ADAPTED ]
C TRANSPOSED S-ELEMENTS.
C
       DO 50 SOUTSPIN = 1, 1 + IUHF
          CALL GETLEN(LENS, VRT(1, SOUTSPIN), VRT(1,ISPIN),
     $       POP(1, SOUTSPIN), NS(1))
       I020 = I010 + LENS *IINTFP
       CALL ZERO(ICORE(I010),LENS)
       IF (ISIDE.EQ.1) THEN
          LISTW = LISTW0(SINSPIN, SOUTSPIN + 1 - IUHF)
       ELSE
          LISTW = LISTW0(SOUTSPIN, SINSPIN + 1 - IUHF)
       ENDIF
CIUHF=1 LISTW(1,1),LISTW(1,2), LIST(2,1),LISTW(2,2)
CIUHF=0 LISTW(1,2)
       BSPIN = ISPIN
C
C LOOP OVER IRREPS OF RHS OF S VECTOR
C
       ICOUNT1 = I000
       ICOUNT2 = I010
       DO 100 XIRREP=1,NIRREP
        BIRREP=DIRPRD(XIRREP,SIRREP)
        IRREPW=XIRREP
        DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW)) 
        NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW)) 
        DISSYS= IRPDPD(IRREPW, 15 + SINSPIN)
        NUMDSS=VRT(BIRREP,BSPIN) * NS(SIRREP)
        I030=I020+IINTFP*DISSYW*NUMDSW
        IF(I030.GT.MAXCOR)THEN
C
C OUT-OF-CORE ALGORITHM
C
         WRITE(6,*)' out-of-core AB not coded, s2hpphs2'
         WRITE(6,1000) I030, MAXCOR
 1000    FORMAT('  REQUIRED MEMORY :', I12, '  AVAILABLE : ', I12)
         call errex
C
        ENDIF
C
C
C DO IN-CORE ALGORITHM
C
C
C READ W INTO W(CK, AI)
C     
         CALL GETLST(ICORE(I020),1,NUMDSW,1,IRREPW,LISTW)
         IF ((LISTW.EQ.54) .OR. (LISTW.EQ.55)) THEN
            CALL SSCAL(NUMDSW*DISSYW, ONEM, ICORE(I020),1)
         ENDIF
C
C PERFORM MATRIX MULTIPLICATION
C
C         Z(AI, BP) = SUM W(CK, AI) * S(CK, BP)         [ISIDE = 1]
C                      CK
C         Z(CK, BP) = SUM W(CK, AI) * S(AI, BP)         [ISIDE = 2]
C                      AI
C
         IF (ISIDE.EQ.1) THEN
            NROW = NUMDSW
         ELSE
            NROW=DISSYW
         ENDIF
          NCOL=NUMDSS
          NSUM=DISSYS
          print = prints .and. (nrow*ncol*nsum .gt. 0)
          if (print) then
             write(6,*)' sinspin', sinspin
             write(6,*)' soutspin', soutspin
             write(6,*)'  input w-coefficients, irrepw =',irrepw
             call output(icore(i020), 1, dissyw, 1, numdsw,
     $          dissyw, numdsw, 1)
C             write(6,*)'  input s-coefficients'
C             call output(icore(icount1), 1, dissys, 1, numdss,
C     $          dissys, numdss,1)
          endif
          IF (ISIDE.EQ.1) THEN
            CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(I020),NSUM,
     &               ICORE(ICOUNT1),NSUM,ONE,ICORE(ICOUNT2),NROW)
          ELSE
            CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(I020),NROW,
     &               ICORE(ICOUNT1),NSUM,ONE,ICORE(ICOUNT2),NROW)
          ENDIF
C          if (print) then
C             write(6,*)' output s-coefficients'
C             call output(icore(icount2), 1, nrow,
C     $          1, numdss, nrow, numdss, 1)
C          endif
          ICOUNT1 = ICOUNT1 + DISSYS * NUMDSS * IINTFP
          ICOUNT2 = ICOUNT2 + NROW * NUMDSS * IINTFP
100    CONTINUE
C
C     REORDER S COEFFICIENTS BACK TO ORIGINAL FORM S[CK, BP] -> S[CB, KP]
C
       I030 = I020 + LENS *IINTFP
       CALL SSTGEN(ICORE(I010), ICORE(I020), LENS, VRT(1,SOUTSPIN),
     $  POP(1, SOUTSPIN), VRT(1, BSPIN), NS(1), ICORE(I030), 1, '1324')
C
C ADD THE S INCREMENTS TO LISTS2EX
C
       LISTS2EX = LS2OUT(ISPIN, SOUTSPIN + 1 - IUHF)
       IF ((IUHF.NE.0) .AND. (SOUTSPIN.EQ.ISPIN)) THEN
          CALL ASSYMALL(ICORE(I020), LENS, NUMSZS, IRPDPD(1, 18+ISPIN),
     $        VRT(1, ISPIN), ICORE(I030), MAXCOR-I030+1)
          CALL GETEXP2(ICORE(I010), LENS, NUMSZS, IRPDPD(1, 18+ISPIN),
     $       LISTS2EX, VRT(1, ISPIN), IRPDPD(1, ISYTYP(1, LISTS2EX)))
       ELSE
          CALL GETALLS2(ICORE(I010), LENS,POP(1,SOUTSPIN),NS(1),
     $       1,LISTS2EX)
       ENDIF       
C
       CALL SAXPY (LENS,ONE,ICORE(I020),1,ICORE(I010),1)
C
       IF ((IUHF.NE.0).AND. (SOUTSPIN.EQ.ISPIN)) THEN
          CALL PUTSQZ(ICORE(I010), LENS, NUMSZS, IRPDPD(1, 18+ISPIN),
     $       LISTS2EX, VRT(1, ISPIN), IRPDPD(1, ISYTYP(1, LISTS2EX)),
     $       ICORE(I020), MAXCOR-I020+1)
       ELSE
          CALL PUTALLS2(ICORE(I010),LENS,POP(1,SOUTSPIN), NS(1),1,
     $       LISTS2EX)
       ENDIF
C
 50    CONTINUE
 25   CONTINUE
      print = .false.
      if (print) then
         write(6,*)' hpph contributions are calculated'
      endif
C
      RETURN
      END
