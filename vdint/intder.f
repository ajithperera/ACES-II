      SUBROUTINE INTDER(IPRINT,NOINT,SEGMEN,WORK1,LWORK1)
C
C     TUH
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      LOGICAL NOINT,SEGMEN
      DIMENSION WORK1(LWORK1)
      COMMON /INTADR/ IWKAO, IWKSO, IWKHHS, IWK1HH, IWK1HC, IWKLST
      LOGICAL         AEQB,   CEQD,   DIAGAB, DIAGCD, DIACAC,
     *                ONECEN, PQSYM,  DTEST,
     *                TPRIAB, TPRICD, TCONAB, TCONCD
      COMMON /INTINF/ THRESH,
     *                NHKTA,  NHKTB,  NHKTC,  NHKTD,
     *                MAXAB,  MAXCD,  JMAX0,
     *                KHKTA,  KHKTB,  KHKTC,  KHKTD,
     *                KHKTAB, KHKTCD, KHABCD,
     *                MHKTA,  MHKTB,  MHKTC,  MHKTD,
     *                MULA,   MULB,   MULC,   MULD,
     *                NORBA,  NORBB,  NORBC,  NORBD, NORBAB, NORBCD,
     *                NUCA,   NUCB,   NUCC,   NUCD,  NUCAB,  NUCCD,
     *                NSETA,  NSETB,  NSETC,  NSETD,
     *                ISTEPA, ISTEPB, ISTEPC, ISTEPD,
     *                NSTRA,  NSTRB,  NSTRC,  NSTRD,
     *                AEQB,   CEQD,
     *                DIAGAB, IAB0X,  IAB0Y,  IAB0Z,
     *                DIAGCD, ICD0X,  ICD0Y,  ICD0Z,
     *                DIACAC, ONECEN, PQSYM,  IPQ0X, IPQ0Y, IPQ0Z,
     *                TPRIAB, TPRICD, TCONAB, TCONCD,
     *                MAXDER, DTEST
      LOGICAL CROSS1, CROSS2
      COMMON /CROSSD/ CROSS1, CROSS2
      LOGICAL         DC101, DC1H1, DC1E1, DC2H1, DC2E1,
     *                DC102, DC1H2, DC1E2, DC2H2, DC2E2,
     *                DPATH1, DPATH2
      COMMON /SUBDIR/ DC101, DC1H1, DC1E1, DC2H1, DC2E1,
     *                DC102, DC1H2, DC1E2, DC2H2, DC2E2,
     *                DPATH1, DPATH2, NTOTAL
      LOGICAL TKTIME
C
C     Initialize NOINT
C
      NOINT = .FALSE.
C
C     ******************************************************
C     ******************************************************
C     ********** Section One: Hermitian Integrals **********
C     ******************************************************
C     ******************************************************
C
      JMAX   = JMAX0 + MAXDER
      JMAXP  = JMAX + 1
      ISTEPT = NUCAB*NUCCD
      ISTEPU = ISTEPT*JMAXP
      ISTEPV = ISTEPU*JMAXP
      NRTUV  = ISTEPV*JMAXP
CSSS      Write(6,*) "in intder", JMAX0, MAXDER, NUCAB, NUCCD

C
C     *******************************
C     ***** RJ(0,0,0) Integrals *****
C     *******************************
C
C     Work space allocation
C
      IWKHHS = IWKAO + MAX(2*ISTEPV,NTOTAL)
      IWK1HH = IWKAO + NTOTAL
      IHHLST = IWK1HH + ISTEPT*(JMAX + 1)*(JMAX + 2)*(JMAX + 3)/6
      IWKMAX = MAX(IHHLST,IWKHHS + 2*NRTUV)
      IF (IPRINT .GE. 10) THEN
         WRITE (LUPRI,'(/1X,A,I10)') ' IWKSO  ', IWKSO
         WRITE (LUPRI,'(/1X,A,I10)') ' IWKAO  ', IWKAO
         WRITE (LUPRI,'(/1X,A,I10)') ' IWKHHS ', IWKHHS
         WRITE (LUPRI,'(/1X,A,I10)') ' IWK1HH ', IWK1HH
         WRITE (LUPRI,'(/1X,A,I10)') ' IHHLST ', IHHLST
         WRITE (LUPRI,'(/1X,A,I10)') ' IWKMAX ', IWKMAX
      END IF
      IF (IWKMAX .GT. LWORK1) THEN
       WRITE (LUPRI, 1000) IWKMAX, LWORK1
       CALL ERREX
      END IF
C
c      IF (TKTIME) TIMSTR = SECOND()
      CALL R000(JMAX,NOINT,ISTEPT,ISTEPU,ISTEPV,NRTUV,NUCAB,NUCCD,
     *          THRESH,ONECEN,IPRINT)
c      IF (TKTIME) THEN
c         TIMEND = SECOND()
c         TIME = TIMEND - TIMSTR
c         TR000X(JMAX) = TR000X(JMAX) + TIME
c         TR000 = TR000 + TIME
c         TIMSTR = TIMEND
c      END IF
      IF (NOINT) RETURN
C
C  RJ(T,U,V) Integrals 
C
      CALL HERI(JMAX,ISTEPT,ISTEPU,ISTEPV,NRTUV,IPQ0X,IPQ0Y,IPQ0Z,
     *          IPRINT,WORK1,LWORK1)
C
C  Section Two: Cartesian Integrals 
C
C   Clear space for final Cartesian integrals 
C
      CALL ZERO(WORK1(IWKAO + 1),NTOTAL)
C
C  Evaluation of integrals 
C
C     Differentiation Path 1
CSSS      Write(6,*) DPATH1, DPATH2, DC2H2, DC2E2, NTOTAL 
C
      IF (DPATH1) THEN
         IPATH = 1
         IWKLST = IHHLST
         CALL C1DRIV(NHKTA,NHKTB,KHKTA,KHKTB,KHKTAB,MHKTA,MHKTB,
     *               NORBA,NORBB,NUCA,NUCB,NSETA,NSETB,
     *               MAXAB,MAXCD,NUCAB,NUCCD,NORBAB,JMAX,
     *               ISTEPA,ISTEPB,ISTEPT,ISTEPU,ISTEPV,
     *               TPRIAB,TCONAB,DIAGAB,PQSYM,IPQ0X,IPQ0Y,IPQ0Z,
     *               IAB0X,IAB0Y,IAB0Z,IPATH,MAXDER,DC101,DC1H1,DC1E1,
     *               CROSS1,SEGMEN,IPRINT,WORK1,LWORK1)
c         IF (TKTIME) THEN
c            TIMEND = SECOND()
c            TCERI1 = TCERI1 + TIMEND - TIMSTR
c            TIMSTR = TIMEND
c         END IF
         CALL C2DRIV(NHKTC,NHKTD,KHKTC,KHKTD,KHKTAB,KHKTCD,KHABCD,
     *               MHKTC,MHKTD,NORBC,NORBD,NUCC,NUCD,NSETC,NSETD,
     *               NUCCD,NORBAB,NORBCD,ISTEPC,ISTEPD,
     *               TPRICD,TCONCD,DIAGCD,
     *               PQSYM,IPQ0X,IPQ0Y,IPQ0Z,ICD0X,ICD0Y,ICD0Z,
     *               IPATH,MAXDER,DC2H1,DC2E1,CROSS1,DTEST,SEGMEN,
     *               IPRINT,WORK1,LWORK1,1)
c         IF (TKTIME) THEN
c            TIMEND = SECOND()
c            TCERI2 = TCERI2 + TIMEND - TIMSTR
c            TIMSTR = TIMEND
c         END IF
      END IF
C
C     Differentiation Path 2
C
      IF (DPATH2) THEN
         IPATH = 2
         IWKLST = IHHLST
         IF (DPATH1) CALL HERSWP(JMAX,NUCAB,NUCCD,IPRINT,WORK1,LWORK1)
         CALL C1DRIV(NHKTC,NHKTD,KHKTC,KHKTD,KHKTCD,MHKTC,MHKTD,
     *               NORBC,NORBD,NUCC,NUCD,NSETC,NSETD,
     *               MAXCD,MAXAB,NUCCD,NUCAB,NORBCD,JMAX,
     *               ISTEPC,ISTEPD,ISTEPT,ISTEPU,ISTEPV,
     *               TPRICD,TCONCD,DIAGCD,PQSYM,IPQ0X,IPQ0Y,IPQ0Z,
     *               ICD0X,ICD0Y,ICD0Z,IPATH,MAXDER,DC102,DC1H2,DC1E2,
     *               CROSS2,SEGMEN,IPRINT,WORK1,LWORK1)
c         IF (TKTIME) THEN
c            TIMEND = SECOND()
c            TCERI1 = TCERI1 + TIMEND - TIMSTR
c            TIMSTR = TIMEND
c         END IF
CSSS         Write(6,*) "Entering C2DRIV ----2", IPATH
         CALL C2DRIV(NHKTA,NHKTB,KHKTA,KHKTB,KHKTCD,KHKTAB,KHABCD,
     *               MHKTA,MHKTB,NORBA,NORBB,NUCA,NUCB,NSETA,NSETB,
     *               NUCAB,NORBCD,NORBAB,ISTEPA,ISTEPB,
     *               TPRIAB,TCONAB,DIAGAB,
     *               PQSYM,IPQ0X,IPQ0Y,IPQ0Z,IAB0X,IAB0Y,IAB0Z,
     *               IPATH,MAXDER,DC2H2,DC2E2,CROSS2,DTEST,SEGMEN,
     *               IPRINT,WORK1,LWORK1,2)
c         IF (TKTIME) TCERI2 = TCERI2 + SECOND() - TIMSTR
      END IF
      RETURN
 1000 FORMAT (//,1X,' WORK SPACE REQUIREMENT ',I9,' EXCEEDS ',
     *        ' CURRENT LIMIT ',I9,' OF WORK1.')
      END
