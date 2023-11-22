      SUBROUTINE INITDAV(NEA, SIRREP, ISPIN, IUHF, SCR, MAXCOR)
C
C  PARAMETERS FOR DAVIDSON PROCEDURE ARE INITIALIZED. ULTIMATELY THIS
C  MUST BE REWRITTEN AND SPECIFIED BY INPUT
C
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION SCR, EIGVAL, THRESH, R, P, ONE, ROOT
      LOGICAL EMINFOL,EVECFOL, LEFTHAND, EXCICORE
C
      PARAMETER (MAXORD=100)
      PARAMETER (MAXROOT=100)
C
      DIMENSION SCR(MAXCOR)
C
      COMMON/LISTDAV/LISTC, LISTHC, LISTH0
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
      COMMON/EXTINF/NDIMR,IOLDEST
      COMMON/EXTINF2/ROOT
      COMMON/EXTINF3/IROOT,LOCROOT,ITROOT
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/RMAT/ R(10000), P(10000)
      COMMON/FLAGS/IFLAGS(100)
      COMMON/ROOTS/EIGVAL(100,8,3), OSCSTR(100,8,3)
      COMMON/CNVRGE/EMINFOL,EVECFOL
      COMMON/EACALC/LEFTHAND, EXCICORE, SINGONLY, DROPCORE
      COMMON/STINFO/ITOTALS, STMODE, LISTST, STCALC, NSIZEST
C
      DATA ONE /1.0D0/
C
      THRESH = -100.0
      STCALC = 2
      NSIZEC = NEA
      NDIMR = 1
      IOLDEST = 1
C
      IROOT = 0
C
      CALL ZERO(R, 10000)
      CALL ZERO(P, 10000)
C
C  CALCULATE THE NEEDED EXCITATION PATTERNS
C
C    COLUMN 5 OF H0  CONTAINS THE INTERESTING COMPONENTS
C    COLUMN 6 OF H0  CONTAINS THE ACTUAL INCLUDED COMPONENTS
C    COLUMN 7 OF H0  CONTAINS THE CLOSED SHELL COMPONENTS OF THE
C                     OTHER SPIN-TYPE (ONLY IF STMODE .NE. 0)
C
      IF (EXCICORE) THEN
          CALL CALCEXCP(ISPIN, IUHF, SCR, MAXCOR, STMODE.NE.0,
     $      .TRUE., 5, LISTH0)
         CALL CALCEXCP(ISPIN, IUHF, SCR, MAXCOR, STMODE.NE.0,
     $      DROPCORE, 6, LISTH0)
         IF (STMODE .NE. 0) THEN
            CALL CALCEXCP(3-ISPIN, IUHF, SCR, MAXCOR, STMODE.NE.0,
     $         .TRUE., 1, LISTST)
         ENDIF
      ELSE
         CALL CALCEXCP(ISPIN, IUHF, SCR, MAXCOR, .FALSE.,
     $      .FALSE., 5, LISTH0)
         CALL CALCEXCP(ISPIN, IUHF, SCR, MAXCOR, .FALSE.,
     $      .FALSE., 6, LISTH0)
      ENDIF
C
      CALL CALCH0(ISPIN, IUHF, SCR, MAXCOR)
C
C LOAD GUESS FOR FIRST VECTOR
C
C  THE INITIAL GUESS IS BASED ON THE DIAGONAL MATRIX H0
C
        I000 = 1
        I010 = I000 + NSIZEC
        I020 = I010 + NSIZEC
        CALL GETLST(SCR(I000), 4,1,1,1, LISTH0)
        CALL FNDMINE(NSIZEC, SCR(I000), EIGVAL(1,SIRREP,ISPIN),
     $     SCR(I010),0, ROOT, ILOC, 1.D-4, THRESH)
        SCR(ILOC) = 1.0D30
        CALL PUTLST(SCR(I000), 4, 1,1,1, LISTH0)
C
        CALL ZERO(SCR,NSIZEC)
        SCR(ILOC) = ONE
        CALL NORMVEC(SCR(I000), NSIZEC, SCR(I010), MAXCOR-I010+1,
     $     IUHF, 2, SIRREP)
        WRITE(6,'(1X,A,E12.5,I6)')' GUESS FOR FIRST EIGENVALUE',
     $     ROOT, ILOC
C
        CALL PUTS(SCR(I000),NSIZEC,ISPIN, IUHF, LS1IN, LS2IN)
C
C  THIS GUESS MAY ALSO BE USED FOR LEFTHAND, AND IS PUT ON 
C  COLUMN 3 OF LISTH0
C
        CALL PUTLST(SCR(I000),3,1,1,1,LISTH0)
C
        RETURN
        END
