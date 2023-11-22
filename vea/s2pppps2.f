      SUBROUTINE S2PPPPS2(ICORE,MAXCOR,IUHF,ISIDE, ISPIN)
C
C  THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF S2 TO
C  S2 INVOLVING THE PPPP BLOCK OF EFFECTIVE INTEGRALS W.
C
C CONTRACTION :
C
C     Z(Ab,Mp) = SUM ,<Ab|Cd> S(Cd,Mp) [SPIN ADAPTED RHF]
C                Cd
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      logical print, prints
C
      DIMENSION ICORE(MAXCOR)
      DIMENSION LISTW0(2,2)
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
C     233: W(AbCd), ordered according to Ab;Cd
C
      prints = .FALSE.
      LISTW0(1,1) = 231
      LISTW0(1,2) = 233
      LISTW0(2,1) = 233
      LISTW0(2,2) = 232
C
C  IF ISPIN =  AND IMIXSPIN =2 : SPECIAL CASE
C
C
C LOOP OVER IRREPS OF RHS OF S VECTOR
C
      DO 50 IMIXSPIN = 1, 1+IUHF
         LISTW = LISTW0(ISPIN, IMIXSPIN + 1 - IUHF)
         MSPIN = IMIXSPIN
         LISTS2EX = LS2OUT(ISPIN, IMIXSPIN + 1 - IUHF)
         LISTS2IN = LS2IN(ISPIN, IMIXSPIN + 1 - IUHF)
       DO 100 MIRREP=1,NIRREP
        XIRREP=DIRPRD(MIRREP,SIRREP)
        IRREPW=XIRREP
        DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW)) 
        NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW)) 
        DISSYS=IRPDPD(XIRREP,ISYTYP(1,LISTS2IN))
        DISSYSEX = IRPDPD(XIRREP, ISYTYP(1, LISTS2EX))
        NUMDSS=POP(MIRREP, MSPIN) * NS(SIRREP)
        print = (prints .and. ((numdss*dissys) .GT. 0))
        I000=1
        I010=I000
        I020=I010+IINTFP*DISSYS*NUMDSS
        I030=I020+IINTFP*DISSYSEX*NUMDSS
        ITOP = I030
        IF ((ISPIN.EQ.1) .AND. (IMIXSPIN.EQ.2)) THEN
           MAXS = MAX(NUMDSS, DISSYS)
           I040 = I030 + MAXS * IINTFP
           I050 = I040 + MAXS * IINTFP
           I060 = I050 + MAXS * IINTFP
           ITOP = I060
        ENDIF
        I070=ITOP+IINTFP*DISSYW*NUMDSW
C
        IF(I070.GT.MAXCOR)THEN
C
C OUT-OF-CORE ALGORITHM
C
         CORLFT = MAXCOR - ITOP + 1
C         write(6,*)' out of core algorithm in pppp ', corlft
C
C  WORKARRAY AVAILABLE FOR W-INTEGRALS
C
         NINCOR = CORLFT/(DISSYW*IINTFP)
         IF (NINCOR.EQ.0) THEN
            WRITE(6,*)' OUT-OF-CORE ALGORITHM S2PPPPS2 FAILS'
            CALL ERREX
         ENDIF
C
         CALL GETLST(ICORE(I010),1,NUMDSS,1,XIRREP,LISTS2IN)
         if (print) then
            write(6,*)' s-coefficients'
            call output(icore(i010), 1, dissys, 1, numdss,
     $         dissys, numdss,1)
         endif         
C
C IF ISPIN = 1 AND IMIXSPIN = 2 THEN SPECIAL MULTIPLICATION:
C
C  Z(aB, iP) = SUM[D,c] (W(Ba, Dc) S(cD, ip)  does not match
C  therefore interchange cD -> S(Dc, ip).
C  Result Z(Ba,iP) -> interchange Ba -> S(aB, iP)
C
         IF ((ISPIN.EQ.1) .AND. (IMIXSPIN.EQ.2)) THEN
            CALL SYMTR3(IRREPW, VRT(1,IMIXSPIN),
     $         VRT(1, ISPIN), DISSYS, NUMDSS, ICORE(I010),
     $         ICORE(I030), ICORE(I040), ICORE(I050))
         if (print) then
            write(6,*)' s-coefficients after symtr3'
            call output(icore(i010), 1, dissys, 1, numdss,
     $         dissys, numdss,1)
         endif         
         ENDIF
         IF (ISIDE .EQ. 1) THEN
            CALL SCOPY(NUMDSS*DISSYS,ICORE(I010),1, ICORE(I020),1)
            CALL TRANSP(ICORE(I020), ICORE(I010), NUMDSS, DISSYS)
            CALL ZERO(ICORE(I020), NUMDSS*DISSYS)
            ISTART = I010
            IFIRST = 1
            NLEFT = NUMDSW
 500        NREAD = MIN(NINCOR,NLEFT)
            CALL GETLST(ICORE(ITOP), IFIRST, NREAD, 1, IRREPW, LISTW)
               NROW = DISSYW
               NCOL = NUMDSS
               NSUM = NREAD
               CALL XGEMM('N','T',NROW,NCOL,NSUM,ONE,ICORE(ITOP),NROW,
     &            ICORE(ISTART),NCOL, ONE, ICORE(I020),NROW)
               ISTART = ISTART + NREAD * NUMDSS * IINTFP
               IFIRST = IFIRST + NREAD
               NLEFT = NLEFT - NREAD
               IF (NLEFT .NE. 0)  GOTO 500
            IF ((IFIRST.NE. NUMDSW+1) .OR. (ISTART.NE.I020)) THEN
               WRITE(6,*)' SOMETHING WRONG IN S2PPPPS2'
               WRITE(6,*) IFIRST, NUMDSW+1, ISTART, I020
               CALL ERREX
            ENDIF
C
         ELSEIF (ISIDE .EQ. 2) THEN
            CALL ZERO(ICORE(I020), NUMDSS*DISSYS)
            ISTART = I020
            IFIRST = 1
            NLEFT = NUMDSW
 600        NREAD = MIN(NINCOR,NLEFT)
            CALL GETLST(ICORE(ITOP), IFIRST, NREAD, 1, IRREPW, LISTW)
               NROW = NUMDSS
               NCOL = NREAD
               NSUM = DISSYW
               CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(I010),NSUM,
     &            ICORE(ITOP),NSUM, ONE, ICORE(ISTART),NROW)
               ISTART = ISTART + NREAD * NUMDSS * IINTFP
               IFIRST = IFIRST + NREAD
               NLEFT = NLEFT - NREAD
               IF (NLEFT .NE. 0 ) GOTO 600
            IF ((IFIRST.NE. NUMDSW+1) .OR. (ISTART.NE.I030)) THEN
               WRITE(6,*)' SOMETHING WRONG IN S2PPPPS2'
               WRITE(6,*) IFIRST, NUMDSW, ISTART, I030
               CALL ERREX
            ENDIF
            CALL SCOPY(NUMDSS*DISSYS,ICORE(I020),1, ICORE(I010),1)
            CALL TRANSP(ICORE(I010), ICORE(I020), DISSYS, NUMDSS)
         ENDIF
C
C  Iside = 1 or Iside =2, out of core algorithm
C
        ELSE
C
C
C DO IN-CORE ALGORITHM
C
C
C READ W, S2IN, S2EX
C     
         CALL GETLST(ICORE(ITOP),1,NUMDSW,1,IRREPW,LISTW)
         if (print) then
            write(6,*)' W-integrals pppp-block'
            write(6,*)'irrepw = ', irrepw
            write(6,*)' imixspin =' , imixspin
            call output(icore(iTOP), 1, dissyw, 1, numdsw,
     $         dissyw, numdsw,1)
         endif
         CALL GETLST(ICORE(I010),1,NUMDSS,1,XIRREP,LISTS2IN)
         if (print) then
            write(6,*)' s-coefficients'
            call output(icore(i010), 1, dissys, 1, numdss,
     $         dissys, numdss,1)
         endif         
C
C IF ISPIN = 1 AND IMIXSPIN = 2 THEN SPECIAL MULTIPLICATION:
C
C  Z(aB, iP) = SUM[D,c] (W(Ba, Dc) S(cD, ip)  does not match
C  therefore interchange cD -> S(Dc, ip).
C  Result Z(Ba,iP) -> interchange Ba -> S(aB, iP)
C
         IF ((ISPIN.EQ.1) .AND. (IMIXSPIN.EQ.2)) THEN
            I040 = I030 + NUMDSS * IINTFP
            I050 = I040 + NUMDSS * IINTFP
            CALL SYMTR3(IRREPW, VRT(1,IMIXSPIN),
     $         VRT(1, ISPIN), DISSYS, NUMDSS, ICORE(I010),
     $         ICORE(I030), ICORE(I040), ICORE(I050))
         if (print) then
            write(6,*)' s-coefficients after symtr3'
            call output(icore(i010), 1, dissys, 1, numdss,
     $         dissys, numdss,1)
         endif         
         ENDIF
C PERFORM MATRIX MULTIPLICATION
C
          NROW=DISSYW
          NCOL=NUMDSS
          NSUM=NUMDSW
          IF (ISIDE.EQ.1) THEN
             CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(ITOP),NROW,
     &               ICORE(I010),NSUM, ZILCH, ICORE(I020),NROW)
          ELSE
C
C ISIDE = 2, LEFT HAND SIDE
C
             CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(ITOP),NROW,
     &               ICORE(I010),NSUM,ZILCH,ICORE(I020),NROW)          
          ENDIF
         if (print) then
            write(6,*)'output s-coefficients'
            call output(icore(i020), 1, dissys, 1, numdss,
     $         dissys, numdss,1)
         endif
         ENDIF
C
C END IN CORE-ALGORITHM
C
          IF ((ISPIN.EQ.1) .AND. (IMIXSPIN.EQ.2)) THEN
             CALL SYMTR3(IRREPW, VRT(1,ISPIN),
     $         VRT(1, IMIXSPIN), DISSYS, NUMDSS, ICORE(I020),
     $         ICORE(I030), ICORE(I040), ICORE(I050))
         if (print) then
            write(6,*)'output s-coefficients after symtr3'
            call output(icore(i020), 1, dissys, 1, numdss,
     $         dissys, numdss,1)
         endif         
          ENDIF
C
          CALL GETLST(ICORE(I010),1,NUMDSS,1,XIRREP,LISTS2EX)
          CALL SAXPY(NUMDSS*DISSYS, ONE, ICORE(I020), 1, ICORE(I010),1)
          CALL PUTLST(ICORE(I010),1, NUMDSS,1,XIRREP,LISTS2EX)
C
 100   CONTINUE
 50    CONTINUE
C
       RETURN
       END
