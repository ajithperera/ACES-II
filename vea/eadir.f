
C
C *******************************************************
C  THE FOLLOWING PROCEDURES FORM THE CORE OF THE EA-EOM PROGRAM. THEY EVALUATE THE
C  PRODUCT OF THE EA-EOM MATRIX TIMES AN ARBITRARY VECTOR RESIDING ON LS1IN AND LS2IN.
C *******************************************************
C
      SUBROUTINE EADIR(ICORE,MAXCOR,IUHF, ISIDE, KSPIN)
C
C   THIS SUBROUTINE CALCULATES THE MATRIX VECTOR PRODUCT
C   OF THE EA-EOM MATRIX AND A VECTOR
C
C  INPUT : COMMON SLISTS
C         LS1IN, 
C         LS2IN : THESE LISTS CONTAIN
C                   THE P AND 2PH COMPONENTS OF THE INPUT VECTOR
C         LS2IN IS LABELED BY (ISPIN, IMIXSPIN) -> S(Ab, Ip)
C         A AND I HAVE IMIXSPIN, B HAS ISPIN
C             ISIDE:   = 1:  EAMAT * V 
C                      = 2:  V * EAMAT
C
C          THE COMMONBLOCK SINFO CONTAINS FURTHER INPUT INFORMATION
C
C           NS:  THE NUMBER OF INPUTVECTORS OF REPRESENTATION SIRREP
C         LENS:  THE LENGTH OF THE 2PH PART OF THE S-VECTOR
C
C  OUTPUT: LISTS1EX, LISTS2EX: LISTS THAT CONTAIN THE RESPECTIVE COMPONENTS 
C                OF THE EA-EOM MATRIX TIMES THE INPUT VECTOR
CEND
      IMPLICIT INTEGER (A-Z)
      logical skip, skip2, skipdir, EXCICORE, LEFTHAND, SINGONLY,
     $   DROPCORE, ADDST
      DOUBLE PRECISION ONE, ONEM
      DOUBLE PRECISION TIN, TOUT,TIMDUM
      DOUBLE PRECISION TSSTGEN, TEADAVID, TEADIR
C
      DIMENSION ICORE(MAXCOR)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SINFO/ NS(8), SIRREP
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
      COMMON/EACALC/LEFTHAND, EXCICORE, SINGONLY, DROPCORE
      COMMON/COREINFO/IREPCORE, SPINCORE, IORBCORE, IORBOCC
      COMMON /TIMSUB/ TSSTGEN, TEADAVID, TEADIR
      COMMON/STINFO/ITOTALS, STMODE, LISTST, STCALC, NSIZEST
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW      
C
      DATA ONE, ONEM /1.0D0, -1.0D0/
C
      CALL TIMER(1)
      ADDST = .FALSE.
      IF (EXCICORE) THEN
         ADDST = ((STMODE.EQ.1) .AND. (KSPIN .EQ. SPINCORE))
      ENDIF
      IF (ADDST) THEN
         ISTART = 1
         IEND = 2
      ELSE
         ISTART = KSPIN
         IEND = KSPIN
      ENDIF
C
      DO 10 ISPIN = 1,1+IUHF
         CALL SETZERO(IUHF, ISPIN, ICORE, MAXCOR, LS1OUT,
     $      LS2OUT)
 10   CONTINUE
      IF ((STMODE .NE. 0) .AND. (ISIDE .EQ. 2)) THEN
         CALL SETZERO(IUHF, 3-KSPIN, ICORE, MAXCOR, LS1IN,
     $      LS2IN)
      ENDIF
C
      DO 20 ISPIN = ISTART, IEND
C
         IF (ADDST .AND. (ISPIN .NE. SPINCORE)) THEN
C
C  EXTRA CONTRIBUTIONS IN THE SPIN-FLIP / SPIN-FLIP TERMS
C            
         ENDIF
C
      skipdir = .false.
      if (skipdir) then
         write(6,*) ' direct terms in eadir are skipped'
      else
      skip = .false.
      if (skip) then
      write(6,*) ' s1 in s1 contribution is skipped'
      else
      CALL S1PPS1(ICORE, MAXCOR, ISIDE, ISPIN)
      endif
C
      skip2 = (.false.)
      if (skip2) then
         write(6,*) 's1 in s2 contribution is skipped'
      else
      CALL S1PPPHS2(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)
      endif
C
      skip2 = (.false.)
      if (skip2) then
         write(6,*) 'one-particle s1 in s2 contribution is skipped'
      else
      CALL S1UHPS2(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)
      endif      
C
      skip2 = (.false.)
      if (skip2) then
         write(6,*) ' phpp s2 -> s1 contribution is skipped'
         else
      CALL S2PHPPS1(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)
      endif
C
      skip2 = (.false.)
      if (skip2) then
         write(6,*) ' one-particle s2 -> s1 contribution is skipped'
      else
         CALL S2UHPS1(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)
      endif
C
      skip = .false.
      if (skip) then
         write(6,*) ' all s2 in s2 contributions are skipped in eadir'
      else
      skip2 = .false.
      if (skip2) then
         write(6,*)' phph-contributions are skipped'
      else
      CALL S2PHPHS2(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)
      endif
C
C THE FOLLOWING SUBROUTINE ALSO EVALUATES PART OF THE S1 CONTRIBUTION
C
      skip2 = .false.
      if (skip2) then
         write(6,*)'hpph-contributions are skipped'
      else
      CALL S2HPPHS2(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)
      endif
C
      skip2 = .false.
      if (skip2) then
         write(6,*)' pppp-contributions are skipped'
      else
         MAXCORS = MAXCOR
      CALL S2PPPPS2(ICORE, MAXCORS, IUHF, ISIDE, ISPIN)
      endif
C
      skip2 = .false.
      if (skip2) then
         write(6,*)' 3-particle-contributions are skipped'
      else
      CALL S2CON3S2(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)
      endif
C
      endif
      endif
C
      skip = .false.
      if (skip) then
         write(6,*) 'singlet - triplet terms are skipped'
      else
         IF (EXCICORE) THEN
            IF ((STMODE .NE. 0) .AND.
     $         (KSPIN .EQ. SPINCORE)) THEN
               CALL SINGTRIP(ICORE, MAXCOR, IUHF, ISIDE, ISPIN)
            ENDIF
         ENDIF
      endif
C
 20   CONTINUE
C
      CALL TIMER(1)      
      TEADIR = TEADIR + TIMENEW
C
      RETURN
      END
