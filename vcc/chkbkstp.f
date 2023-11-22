      SUBROUTINE CHKBKSTP(SCR, MAXCOR, IUHF, ICONVG)
C
C This routine shceks to see if we should stop the CC iterations and
C update the orbitals in a Brueckner calculation.  It will also dump out
C the T2 coefficients to be used as an initial guess for the next CC
C calculation.
C
C SG 7/21/98
C
      IMPLICIT NONE
C
      INTEGER MAXCOR, IUHF, ICONVG
      DOUBLE PRECISION SCR(MAXCOR)
C
      INTEGER IFLAGS
      COMMON /FLAGS/ IFLAGS(100)
      INTEGER POP, VRT, NT, NFMI, NFEA
      INTEGER IINTLN, IFLTLN, IINTFP, IALONE, IBITWD
      COMMON /MACHSP/ IINTLN, IFLTLN, IINTFP, IALONE, IBITWD
      COMMON /SYM/ POP(8,2), VRT(8,2), NT(2), NFMI(2), NFEA(2)
      INTEGER IRPDPD, ISYTYP, ID
      COMMON /SYMPOP/ IRPDPD(8,22), ISYTYP(2,500), ID(18)
C
      INTEGER ISYMSZ,i
      DOUBLE PRECISION FNDLRGAB,d2
C
      DOUBLE PRECISION TEN
      PARAMETER (TEN = 10.0D0)
C
      INTEGER I000, I010, I020, I030, I040, I050, LEN, NAMLEN
      DOUBLE PRECISION TOL, T2CNV
      LOGICAL YESNO, PCCD
      CHARACTER *80 FULNAM
C
      IF (ICONVG .EQ. 0) RETURN
C
      TOL = TEN**(-IFLAGS(76))
C
      I000 = 1
      I010 = I000 + NT(1)
      IF (IUHF.NE.0) I020 = I010 + NT(2)
      LEN = NT(1) + IUHF * NT(2)
      CALL GETLST(SCR(I000), 1, 1, 1, 1, 90)
      IF (IUHF.NE.0) CALL GETLST(SCR(I010), 1, 1, 1, 2, 90)
C
C Test to see if the largest T1 coefficient is larger than the tolerance.
C If it is, then set up to rotate orbitals and restart CC iterations.
C
      CALL GETREC(20, 'JOBARC', 'T2CNVCRT', IINTFP, T2CNV)
      d2 = FNDLRGAB(SCR(I000), LEN) 
      IF (FNDLRGAB(SCR(I000), LEN) .GT. MAX(TOL, TEN*TEN*T2CNV)) THEN
        ICONVG = 0
        IF (IUHF .NE. 0) THEN
          I030 = I020 + ISYMSZ(ISYTYP(1,45), ISYTYP(2,45))
        ELSE
          I030 = I010
        ENDIF
        I040 = I030 + ISYMSZ(ISYTYP(1,44), ISYTYP(2,44))
        I050 = I040 + ISYMSZ(ISYTYP(1,46), ISYTYP(2,46))
        IF (I050 .GT. MAXCOR) THEN
          WRITE(6,1000)
 1000     FORMAT(T3, '@CHKBKSTP-I, not enough memory to dump T',
     &       ' as initial guess')
        ELSE
          CALL ZERO(SCR(I000), LEN)
          IF (IUHF .NE. 0) THEN
            CALL GETALL(SCR(I020), ISYMSZ(ISYTYP(1,45), ISYTYP(2,45)),
     &         1, 45)
          ENDIF
          CALL GETALL(SCR(I030), ISYMSZ(ISYTYP(1,44), ISYTYP(2,44)),
     &       1, 44)
          CALL GETALL(SCR(I040), ISYMSZ(ISYTYP(1,46), ISYTYP(2,46)),
     &       1, 46)
          CALL DMPTGSS(SCR(I000), I050-I000, 'TGUESS  ')
        ENDIF
      ELSE
C
C Make sure that there is not an old TGUESS file lying around
C
        CALL GFNAME('TGUESS  ', FULNAM, NAMLEN)
        INQUIRE(FILE=FULNAM(1:NAMLEN), EXIST=YESNO)
        IF (YESNO) THEN
          OPEN(UNIT=94, FILE=FULNAM(1:NAMLEN), STATUS='OLD',
     &       FORM='UNFORMATTED', ACCESS='SEQUENTIAL')
          CLOSE(UNIT=94, STATUS='DELETE')
        ENDIF
      ENDIF
C
      RETURN
      END

