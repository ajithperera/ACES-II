      SUBROUTINE T1RING(ICORE,MAXCOR,IUHF,LAMBDA)
C
C DRIVER FOR W(MBEJ) <- T1 CONTRIBUTIONS.  SELECTS BETWEEN INCORE
C AND OUT OF CORE ALGORITHMS AND CALLS THE APPROPRIATE ROUTINES
C
C To use incore algorithms, there should be at least enough memory to keep 
C I(mb,ef), T(j,f) and W(mb,ej) at the same time. 
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      LOGICAL LAMBDA,AOBASIS,AOLOG
C
      DIMENSION ICORE(MAXCOR)
C
      COMMON /MACHSP/ IINTLN,IFLTLn,iintfp,ialone,ibitwd
      common /sym/ pop(8,2),vrt(8,2),nt(2),nfea(2),nfmi(2)
      common /sympop/ irpdpd(8,22),isytyp(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /AOLOG/ AOLOG
C
      DATA ONE  /1.0/   
      DATA ONEM /-1.0/
      DATA ZILCH /0.0/
C
      AOBASIS=IFLAGS(93).EQ.2.AND.(.NOT.LAMBDA)
C
C IF CALCULATION IN AOBASIS FIRST READ THE F(e;a)
C
      AOLOG=.FALSE.
      IF (AOBASIS) THEN
         AOLOG=.TRUE.
         I000=1
         DO 1111 ISPIN=1,1+IUHF
            CALL UPDMOI(1,NFMI(ISPIN),6+ISPIN,92,0,0)
            CALL GETLST(ICORE(I000),1,1,1,ISPIN,92)
            CALL PUTLST(ICORE(I000),1,1,1,6+ISPIN,92)
 1111    CONTINUE
      ENDIF
C
C SPIN CASES AAAA AND BBBB
C
      INEED=-1
      DO 10 ISPIN=1,1+IUHF
         LSTTAR=53+ISPIN
         LSTINT1=6+ISPIN
         LSTINT2=26+ISPIN
         FULTYP1=18+ISPIN
         FULTYP2=20+ISPIN
C
         T1SIZE=NT(1)+IUHF*NT(2)
         ABFULL=IRPDPD(1,FULTYP1)
         MNFULL=IRPDPD(1,FULTYP2)
C
         DO 20 IRREP=1,NIRREP
            TARDSZ=IRPDPD(IRREP,ISYTYP(1,LSTTAR))
            TARDIS=IRPDPD(IRREP,ISYTYP(2,LSTTAR))
C
            TARSIZ=TARDSZ*TARDIS
            FULDSZ1=IRPDPD(IRREP,FULTYP1)
            FULDSZ2=IRPDPD(IRREP,FULTYP2)
            INTDSZ1=IRPDPD(IRREP,ISYTYP(1,LSTINT1))
            INTDIS1=IRPDPD(IRREP,ISYTYP(2,LSTINT1))
            INTDSZ2=IRPDPD(IRREP,ISYTYP(1,LSTINT2))
            INTDIS2=IRPDPD(IRREP,ISYTYP(2,LSTINT2))
            I1=TARDIS*TARSIZ+T1SIZE
            I2=3*TARDSZ
C
            I3A=FULDSZ1*INTDIS1+INTDSZ1
            I3B=FULDSZ2*INTDIS2+INTDSZ2
C
            I4A=2*ABFULL
            I4B=2*MNFULL
C
            INEEDA=I1+MAX(I2,I3A,I4A)
            INEEDB=I1+MAX(I2,I3B,I4B)
            INEED =MAX(INEED,INEEDA,INEEDB)
 20      CONTINUE
 10   CONTINUE
C
      IF(INEED.LT.MAXCOR/IINTFP)THEN
C
         CALL T1RA01(ICORE,MAXCOR,IUHF,LAMBDA)
C
      ELSE
C
         CALL T1RA02(ICORE,MAXCOR,IUHF,LAMBDA)
C
      ENDIF
C
C SPIN CASES ABBA AND BAAB
C
      INEED=-1
      DO 110 ISPIN=1,1+IUHF
         LSTTAR=57-ISPIN
         LSTINT1=31-ISPIN
         LSTINT2=8+ISPIN+(1-IUHF)
C
         T1SIZE=NT(1)+IUHF*NT(2)
         FULTYP1=18+ISPIN
         FULTYP2=20+ISPIN
         ABFULL=IRPDPD(1,FULTYP1)
C
         DO 120 IRREP=1,NIRREP
            TARDSZ=IRPDPD(IRREP,ISYTYP(1,LSTTAR))
            TARDIS=IRPDPD(IRREP,ISYTYP(2,LSTTAR))
C
            TARSIZ=TARDSZ*TARDIS
            INTDSZ1=IRPDPD(IRREP,ISYTYP(1,LSTINT1))
            INTDIS1=IRPDPD(IRREP,ISYTYP(2,LSTINT1))
            INTDSZ2=IRPDPD(IRREP,ISYTYP(1,LSTINT2))
            INTDIS2=IRPDPD(IRREP,ISYTYP(2,LSTINT2))
            I1=TARDIS*TARSIZ+T1SIZE
            I2=3*TARDSZ
C
            I3A=INTDSZ1*INTDIS1+3*MAX(INTDSZ1,INTDIS1)
            I3B=INTDSZ2*INTDIS2+3*MAX(INTDSZ2,INTDIS2)
C
            I4=2*ABFULL
C
            I5=TARDIS*TARDSZ+3*TARDSZ

            INEEDA=I1+MAX(I2,I3A,I4,I5)
            INEEDB=I1+MAX(I2,I3B,I4,I5)
            INEED =MAX(INEED,INEEDA,INEEDB)
 120     CONTINUE
 110  CONTINUE
C
      IF(INEED.LT.MAXCOR/IINTFP)THEN
C
         CALL T1RB01(ICORE,MAXCOR,IUHF,LAMBDA)
C
      ELSE
C
         CALL T1RB02(ICORE,MAXCOR,IUHF,LAMBDA)
C
      ENDIF
C
C STORE THE VALUES OF SUM(m,f) t(m;f) <ma||fe> IF AOBASIS IS TRUE
C DO THIS BY SUBTRACTING THE ACTUAL F(e;a)'s 
C BY THE F(e;a)'s READ IN BEFORE
C
      IF (AOBASIS) THEN
         I000=1
         DO 2222 ISPIN=1,1+IUHF
            I010=I000+NFMI(ISPIN)*IINTFP
            CALL GETLST(ICORE(I000),1,1,1,ISPIN,92)
            CALL GETLST(ICORE(I010),1,1,1,6+ISPIN,92)
            CALL SAXPY(NFMI(ISPIN),ONEM,ICORE(I010),1,ICORE(I000),1)
            CALL PUTLST(ICORE(I000),1,1,1,6+ISPIN,92)
 2222    CONTINUE
      ENDIF
C
C SPIN CASES ABAB AND BABA (FOR UHF ONLY)
C
c     IF(IUHF.EQ.0 )RETURN
C
      INEED=-1
      DO 210 ISPIN=1,2
         LSTTAR=55+ISPIN
         LSTTMP1=37+ISPIN
         LSTTMP2=40-ISPIN
         LSTINT1=28+ISPIN
         LSTINT2=8+ISPIN
C
         T1SIZE=NT(1)+IUHF*NT(2)
         DO 220 IRREP=1,NIRREP
            TARDSZ=IRPDPD(IRREP,ISYTYP(1,LSTTAR))
            TARDIS=IRPDPD(IRREP,ISYTYP(2,LSTTAR))
C
            TARSIZ=TARDSZ*TARDIS
            TMPDSZ1=IRPDPD(IRREP,ISYTYP(1,LSTTMP1))
            TMPDIS1=IRPDPD(IRREP,ISYTYP(2,LSTTMP1))
            TMPDSZ2=IRPDPD(IRREP,ISYTYP(1,LSTTMP2))
            TMPDIS2=IRPDPD(IRREP,ISYTYP(2,LSTTMP2))
            INTDSZ1=IRPDPD(IRREP,ISYTYP(1,LSTINT1))
            INTDIS1=IRPDPD(IRREP,ISYTYP(2,LSTINT1))
            INTDSZ2=IRPDPD(IRREP,ISYTYP(1,LSTINT2))
            INTDIS2=IRPDPD(IRREP,ISYTYP(2,LSTINT2))
            I1=TARDIS*TARSIZ+T1SIZE
            I2=3*TARDSZ
C
            I3A=INTDSZ1*INTDIS1+3*MAX(INTDSZ1,INTDIS1)
            I3B=INTDSZ2*INTDIS2+3*MAX(INTDSZ2,INTDIS2)
C
            I4=TARDIS*TARDSZ+3*TARDSZ
C
            I5A=3*MAX(TMPDSZ1,TMPDIS1)
            I5B=3*MAX(TMPDSZ2,TMPDIS2)

            INEEDA=I1+MAX(I2,I3A,I4,I5A)
            INEEDB=I1+MAX(I2,I3B,I4,I5B)
            INEED =MAX(INEED,INEEDA,INEEDB)
 220     CONTINUE
 210  CONTINUE
C
      IF(INEED.LT.MAXCOR/IINTFP)THEN
C
         CALL T1RC01(ICORE,MAXCOR,IUHF,LAMBDA)
C
      ELSE
C
         CALL T1RC02(ICORE,MAXCOR,IUHF,LAMBDA)
C
      ENDIF
C
      RETURN
      END
