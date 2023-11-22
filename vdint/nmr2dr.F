      SUBROUTINE NMR2DR(FOCK,JATOM,IODV,NOPV,
     &                  JPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,
     &                  RETUR,WORK1,LWORK1,MSZ,
     &                  NFIELD,DHFA,DHFB,DRELA,DRELB,DNHFA,DNHFB,
     &                  DHFFA,DHFFB,DRELFA,DRELFB,NBAST,NIR)
C
C PROGRAM FOR CALCULATING TWO-ELECTRON INTEGRAL DERIVATIVES
C NEEDED IN GIAO CALCULATIONS. THE ALGORITHM IS BASED ON
C THE MCMURCHIE-DAVIDSON SCHEME.
C
CEND
C
C SEPTEMBER/91 JG 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
      INTEGER DIRPRD
      LOGICAL SHAEQB, SHCEQD, SHAEQC, SHBEQD, SHAAAA, SHABAB,
     &        AACDX, AAAAX, AACDY, AAAAY, AACDZ, AAAAZ,
     &        ABCCX, ABCCY, ABCCZ, DERECC, DERECD,
     &        PRINT1, PRINT2, PRINT3, PRINT4,NOPV,NOPMAT, NODV,
     &        PERTUR, EXPECT, FIRST, NOINT, INTS, UNDIFF,
     &        CAEQCB, CCEQCD, RETUR, DIFX, DIFY, DIFZ, TTIME,
     &        FOCK,GRAD,INTWRIT,IDFGH,JUNK,INTWRIT2
      LOGICAL SHARE,BIGVEC,SEGMEN
      LOGICAL DCORD,DCORGD,NOORBT,DOPERT
      LOGICAL AEQB,CEQD,DIAGAB,DIAGCD,DIACAC,ONECEN,PQSYM,DTEST,
     &        TPRIAB,TPRICD,TCONAB,TCONCD
      LOGICAL SCF,NONHF
      LOGICAL INIT
      LOGICAL IFIRST,IALL,XCOMP,YCOMP,ZCOMP
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "mxcent.par"
#include "baslims.par"
      PARAMETER (AZERO = 0.00 D00, ONE = 1.00 D00, LUPAO = 25)
C
      DIMENSION XAND(0:7),YAND(0:7),ZAND(0:7)
      DIMENSION NUCOA(MXAOVC),  NUCOB(MXAOVC),
     &          NUCOC(MXAOVC),  NUCOD(MXAOVC),
     &          NRCOA(MXCONT),  NRCOB(MXCONT),
     &          NRCOC(MXCONT),  NRCOD(MXCONT),
     &          JSTRA(MXAOVC),  JSTRB(MXAOVC),
     &          JSTRC(MXAOVC),  JSTRD(MXAOVC),
     &          ICORBA(MXAOVC), ICORBB(MXAOVC),
     &          ICORBC(MXAOVC), ICORBD(MXAOVC)
      DIMENSION DHFA(NBAST,NBAST,0:NIR),
     &          DHFB(NBAST,NBAST,0:NIR),
     &          DRELA(NBAST,NBAST,0:NIR),
     &          DRELB(NBAST,NBAST,0:NIR),
     &          DNHFA(NBAST,NBAST,0:NIR),
     &          DNHFB(NBAST,NBAST,0:NIR)
      DIMENSION DHFFA(NBAST*(NBAST+1)/2,0:NIR),
     &          DHFFB(NBAST*(NBAST+1)/2,0:NIR),
     &          DRELFA(NBAST*(NBAST+1)/2,0:NIR),
     &          DRELFB(NBAST*(NBAST+1)/2,0:NIR)
      DIMENSION  WORK1(LWORK1)
      DIMENSION MSZ(3,3)
C
      COMMON/CCOM/THRS,NHTYP,IBX
      COMMON /SHELLSi/ KMAX,
     &                NHKT(MXSHEL),   KHKT(MXSHEL), MHKT(MXSHEL),
     &                ISTBAO(MXSHEL), NUCO(MXSHEL), JSTRT(MXSHEL),
     &                NSTRT(MXSHEL),  MST(MXSHEL),  NCENT(MXSHEL),
     &                NRCO(MXSHEL), NUMCF(MXSHEL),
     &                NBCH(MXSHEL),   KSTRT(MXSHEL)
      COMMON /SHELLS/ CENT(MXSHEL,3), SHARE(MXSHEL)
      COMMON/BLOCKS/CENTSH(MXSHEL,3),
     &              MAXSHL,BIGVEC,SEGMEN,
     &              NHKTSH(MXSHEL),KHKTSH(MXSHEL),MHKTSH(MXSHEL),
     &              ISTBSH(MXSHEL),NUCOSH(MXSHEL),NORBSH(MXSHEL),
     &              NSTRSH(MXSHEL),NCNTSH(MXSHEL),NSETSH(MXSHEL),
     &              JSTRSH(MXSHEL,MXAOVC),
     &              NPRIMS(MXSHEL,MXAOVC),
     &              NCONTS(MXSHEL,MXAOVC),
     &              IORBSH(MXSHEL,MXAOVC),
     &              IORBSB(MXCORB),NRCSH(MXSHEL)
      COMMON/SYMMET/FMULT(0:7),PT(0:7),
     &              MAXLOP,MAXLOT,MULT(0:7),ISYTYP(3),
     &              ITYPE(8,36),NPARSU(8),NPAR(8),NAOS(8),
     &              NPARNU(8,8),IPTSYM(MXCORB,0:7),
     &              IPTCNT(3*MXCENT,0:7),NCRREP(0:7),
     &              IPTCOR(MXCENT*3),NAXREP(0:7),IPTAX(3),
     &              IPTXYZ(3,0:7)
      CHARACTER NAMEX*6
      COMMON /NUCLEIi/ NOORBT(MXCENT),
     &                NUCIND, NUCDEP, NUCPRE(MXCENT), NUCNUM(MXCENT,8),
     &                NUCDEG(MXCENT), ISTBNU(MXCENT), NDCORD,
     &                NDCOOR(MXCOOR), NTRACO, NROTCO, ITRACO(3),
     &                IROTCO(3),
     &                NATOMS, NFLOAT,
     &                IPTGDV(3*MXCENT),
     &                NGDVEC(8), IGDVEC(8)
      COMMON /NUCLEI/ CHARGE(MXCENT), CORD(MXCENT,3),
     &                DCORD(MXCENT,3),DCORGD(MXCENT,3),
     &                DOPERT(0:3*MXCENT)
      COMMON/NUCLEC/NAMEX(MXCOOR)
      COMMON/INTADR/IWKAO,IWKSO,IWKHHS,IWK1HH,IWK1HC,IWKLST
      COMMON/INTINF/THRESH,
     &              NHKTA,NHKTB,NHKTC,NHKTD,
     &              MAXAB,MAXCD,JMAX0,
     &              KHKTA,KHKTB,KHKTC,KHKTD,
     &              KHKTAB,KHKTCD,KHABCD,
     &              MHKTA,MHKTB,MHKTC,MHKTD,
     &              MULA,MULB,MULC,MULD,
     &              NORBA,NORBB,NORBC,NORBD,NORBAB,NORBCD,
     &              NUCA,NUCB,NUCC,NUCD,NUCAB,NUCCD,
     &              NSETA,NSETB,NSETC,NSETD,
     &              ISTEPA,ISTEPB,ISTEPC,ISTEPD,
     &              NSTRA,NSTRB,NSTRC,NSTRD,
     &              AEQB,CEQD,
     &              DIAGAB,IAB0X,IAB0Y,IAB0Z,
     &              DIAGCD,ICD0X,ICD0Y,ICD0Z,
     &              DIACAC,ONECEN,PQSYM,IPQ0X,IPQ0Y,IPQ0Z,
     &              TPRIAB,TPRICD,TCONAB,TCONCD,
     &              MAXDER,DTEST
      COMMON/IPRI/IPRINT
      COMMON/XYZDIF/XDIFF1,YDIFF1,ZDIFF1,XDIFF2,YDIFF2,ZDIFF2
      COMMON/METHOD/IUHF,SCF,NONHF
      COMMON/BASSYM/NBAS(8),NBASIS,NBASSQ,NBASTT
      COMMON/BASSPH/NBAS5(8),NBASIS5,NBASSQ5,NBASTT5
      COMMON/LSYM/NLENQ(8),NLENT(8)
      COMMON/LSYMSPH/NLENQ5(8),NLENT5(8)
      COMMON/OPTION/GRAD,INTWRIT,JUNK
      COMMON/PERT/NTPERT,NPERT(8),IIPERT(8),IXPERT,IYPERT,IZPERT,
     &            IXYPERT(3)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/DFGH/IDFGH 
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/GOFF/IGOFF
      COMMON/TREATP/IFIRST,IALL,XCOMP,YCOMP,ZCOMP,IBF,IBL
C
      DATA ONEM/-1.D0/
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
      PTAND(I,J) = PT(IBTAND(I,J))
      IBTAXO(I,J) = IBTAND(I,IBTXOR(I,J))
C
C SET MAXDIF EQUAL ONE FOR GIAO CALCULATIONS
C
      MAXDIF=1
      FIRST  = .TRUE.
      NOPMAT = NOPV
C
C  DETERMINE RUN TYP
C
      DTEST  = .FALSE.
      PERTUR = .FALSE.
      UNDIFF = .FALSE.
C
C  EXPECTATION VALUES ARE CALCULATED FOR ALL ATOMS
C
      EXPECT = .TRUE.
      MAXDER = MAXDIF
      IATOM = 0
      MULTE = 1
      THRESH = DMAX1(THRS,1.00D-16)
      DO 150 I = 0,MAXLOP
       XAND(I) = PTAND(ISYTYP(1),I)
       YAND(I) = PTAND(ISYTYP(2),I)
       ZAND(I) = PTAND(ISYTYP(3),I)
  150 CONTINUE
C
C  REWIND LUPAO (ONLY NECESSARY FOR EVALUATION OF
C  MAGNETIC SUSZEPRIBIILITIES
C
c      IF((.NOT.SCF).AND. (.NOT.FOCK)) THEN
c       REWIND LUPAO
c      END IF
C
C ALLOCATE CORE MEMORY FOR DERIVATIVES OF FOCK MATRIX 
C
      KSTART=0
      IF(FOCK) THEN
C
       KDFAOA=KSTART+1
       NF=NBAST*(NBAST+1)/2
       KDFAOB=KDFAOA+NF*NFIELD*IUHF*(NIR+1)
       KSTART=KDFAOB+NF*NFIELD*(NIR+1)-1
       CALL ZERO(WORK1(KDFAOA),(IUHF+1)*NF*NFIELD*(NIR+1))
       IF(.NOT.SCF) THEN
C
C THERE ARE A TOTAL OF 3 (NIRREP=1), 9 (NIRREP=2), 12 (NIRREP GT 2)
C BUFFERS
C
        NUMBUF=3
        IF(NIRREP.GT.1) NUMBUF=9
        IF(NIRREP.GT.2) NUMBUF=12
        IRECL=8192
        NBUF=(IRECL-2)/(1+IINTFP)
C
C DETERMINE OPTIMAL RECORD LENGTH IN THE FUTURE (NBUF = ....)
C
        KBUCK=KSTART+1+NF*NFIELD*(NIR+1)
        KIBUCK=KBUCK+NBUF*NUMBUF
        KCOUNT=KIBUCK+NBUF*NUMBUF
        KREORD=KCOUNT+NUMBUF
        KSTART=KREORD+NIRREP*NBASIS-1
        CALL IZERO(WORK1(KCOUNT),NUMBUF)
        CALL IZERO(WORK1(KREORD),NBASIS*NIRREP)
C 
C OPEN FILE FOR GIAO 2E-INTEGRALS
C
        CALL NMROPEN
C
C FILL REORDER VECTOR WHICH GIVES THE CORRECT INDICES FOR SAOS
C
        CALL SAOORD(WORK1(KREORD),NBASIS,NIRREP)
C
       ENDIF
      ENDIF
C
C ZERO ARRAY FOR MAGNETIC SUSCEPTIBILITY TENSOR
C
      CALL ZERO(MSZ,9)
C
C  FIRST SHELL INDEX
C
      DO 100 ISHELA = 1,MAXSHL
       NHKTA = NHKTSH(ISHELA)
       KHKTA = KHKTSH(ISHELA)
       MHKTA = MHKTSH(ISHELA)
       NCENTA = NCNTSH(ISHELA)
       ICENTA = NUCNUM(NCENTA,1)
       CORAX = CENTSH(ISHELA,1)
       CORAY = CENTSH(ISHELA,2)
       CORAZ = CENTSH(ISHELA,3)
       MULA  = ISTBSH(ISHELA)
       MULTA = MULT(MULA)
       NSTRA = IORBSB(IORBSH(ISHELA,1)) - 1
       NSETA = NSETSH(ISHELA)
       DO 104 I = 1, NSETA
        NUCOA(I) = NPRIMS(ISHELA,I)
        NRCOA(I) = NCONTS(ISHELA,I)
        JSTRA(I) = JSTRSH(ISHELA,I)
104    CONTINUE
       NORBA = NORBSH(ISHELA)
       DO 105 I = 1, NORBA
        ICORBA(I) = IORBSH(ISHELA,I) - 1
105    CONTINUE
       PRINT1 = .TRUE.
       SIGNAX = ONE
       SIGNAY = ONE
       SIGNAZ = ONE
C
C  SECOND SHELL INDEX
C
       DO 200 ISHELB = 1,ISHELA
        SHAEQB = ISHELA .EQ. ISHELB
        NHKTB  = NHKTSH(ISHELB)
        KHKTB  = KHKTSH(ISHELB)
        MHKTB  = MHKTSH(ISHELB)
        NCENTB = NCNTSH(ISHELB)
        CORBX0 = CENTSH(ISHELB,1)
        CORBY0 = CENTSH(ISHELB,2)
        CORBZ0 = CENTSH(ISHELB,3)
        NSTRB = IORBSB(IORBSH(ISHELB,1)) - 1
        NSETB = NSETSH(ISHELB)
        DO 204 I = 1, NSETB
         NUCOB(I) = NPRIMS(ISHELB,I)
         NRCOB(I) = NCONTS(ISHELB,I)
         JSTRB(I) = JSTRSH(ISHELB,I)
204     CONTINUE
        NORBB = NORBSH(ISHELB)
        DO 205 I = 1, NORBB
         ICORBB(I) = IORBSH(ISHELB,I) - 1
205     CONTINUE
        MULB   = ISTBSH(ISHELB)
        MULTB  = MULT(MULB)
        MULTAB = MULTA*MULTB
        KHAB   = KHKTA*KHKTB
        KPQ    = IBTAND(MULA,MULB)
        MPQ    = IBTOR(MULA,MULB)
        MAXAB  = NHKTA + NHKTB - 2
C
C CHARGE DISTRIBUTIONS FOR FIRST ELECTRON 
C
        AEQB = SHAEQB .AND. MULTA .EQ. 1
c        TPRIAB = AEQB
        TPRIAB=.FALSE.
        TCONAB = SHAEQB .AND. MAXAB .EQ. 0
        IF (TCONAB) THEN
         NORBAB = NORBA*(NORBA + 1)/2
        ELSE
         NORBAB = NORBA*NORBB
        END IF
        IF (MULTAB .EQ. 1) THEN
         CAEQCB = ICENTA.EQ.NUCNUM(NCENTB,1)
         CALL GIADRV(NHKTA,NHKTB,NSETA,NSETB,ISTEPA,ISTEPB,
     &                NUCA,NUCB,NUCAB,NORBA,NORBB,
     &                NUCOA,NUCOB,NRCOA,NRCOB,JSTRA,JSTRB,
     &                ONE,ONE,ONE,ONE,ONE,ONE,
     &                CORAX,CORAY,CORAZ,CORBX0,CORBY0,CORBZ0,
     &                AACDX,AACDY,AACDZ,IAB0X,IAB0Y,IAB0Z,TPRIAB,
     &                CAEQCB,.TRUE.,.TRUE.,BIGVEC,SEGMEN,DTEST,
     &                12,THRESH,MAXDER,IPRINT)
         XDIFF1=CORAX-CORBX0
         YDIFF1=CORAY-CORBY0
         ZDIFF1=CORAZ-CORBZ0
         IF (NUCAB .EQ. 0) GO TO 200
        END IF
        DIAGAB = SHAEQB .AND. .NOT.BIGVEC
        IF (DIAGAB) THEN
         KHKTAB = KHKTA*(KHKTA + 1)/2
        ELSE
         KHKTAB = KHKTA*KHKTB
        END IF
C
C  THIRD SHELL INDEX
C
        DO 300 ISHELC = 1,ISHELA
         SHAEQC = ISHELA .EQ. ISHELC
         NHKTC  = NHKTSH(ISHELC)
         NCENTC = NCNTSH(ISHELC)
         CORCX0 = CENTSH(ISHELC,1)
         CORCY0 = CENTSH(ISHELC,2)
         CORCZ0 = CENTSH(ISHELC,3)
         NSTRC = IORBSB(IORBSH(ISHELC,1)) - 1
         NSETC = NSETSH(ISHELC)
         DO 304 I = 1, NSETC
          NUCOC(I) = NPRIMS(ISHELC,I)
          NRCOC(I) = NCONTS(ISHELC,I)
          JSTRC(I) = JSTRSH(ISHELC,I)
  304    CONTINUE
         NORBC = NORBSH(ISHELC)
         DO 305 I = 1, NORBC
          ICORBC(I) = IORBSH(ISHELC,I) - 1
  305    CONTINUE
         KHKTC = KHKTSH(ISHELC)
         MHKTC = MHKTSH(ISHELC)
         MULC  = ISTBSH(ISHELC)
         MULTC = MULT(MULC)
         KPQR  = IBTAND(KPQ,MULC)
         MPQR  = IBTOR(MPQ,MULC)
C
C   FOURTH SHELL INDEX
C
         IDMAX = ISHELC
         IF (SHAEQC) IDMAX = ISHELB
         DO 400 ISHELD = 1,IDMAX
          SHBEQD  = ISHELB .EQ. ISHELD
          SHCEQD  = ISHELC .EQ. ISHELD
          NCENTD = NCNTSH(ISHELD)
          CORDX0 = CENTSH(ISHELD,1)
          CORDY0 = CENTSH(ISHELD,2)
          CORDZ0 = CENTSH(ISHELD,3)
          NSTRD = IORBSB(IORBSH(ISHELD,1)) - 1
          NSETD = NSETSH(ISHELD)
          DO 404 I = 1, NSETD
           NUCOD(I) = NPRIMS(ISHELD,I)
           NRCOD(I) = NCONTS(ISHELD,I)
           JSTRD(I) = JSTRSH(ISHELD,I)
404       CONTINUE
          NORBD = NORBSH(ISHELD)
          DO 405 I = 1, NORBD
           ICORBD(I) = IORBSH(ISHELD,I) - 1
405       CONTINUE
          NHKTD = NHKTSH(ISHELD)
          KHKTD = KHKTSH(ISHELD)
          MHKTD = MHKTSH(ISHELD)
          MULD  = ISTBSH(ISHELD)
          MULTD = MULT(MULD)
          MAXCD = NHKTC + NHKTD - 2
          JMAX0 = MAXAB + MAXCD
          KRS    = IBTAND(MULC,MULD)
          SHAAAA = SHAEQB .AND. SHCEQD .AND. SHAEQC
          SHABAB = SHAEQC .AND. ISHELB .EQ. ISHELD
          CEQD   = SHCEQD .AND. MULTC .EQ. 1
C
          TCONCD = SHCEQD .AND. MAXCD .EQ. 0
          IF (TCONCD) THEN
           NORBCD = NORBC*(NORBC + 1)/2
          ELSE
           NORBCD = NORBC*NORBD
          END IF
C
          DIAGCD = SHCEQD .AND. .NOT.BIGVEC
          IF (DIAGCD) THEN
           KHKTCD = KHKTC*(KHKTC + 1)/2
          ELSE
           KHKTCD = KHKTC*KHKTD
          END IF
C
          MRS    = IBTOR(MULC,MULD)
          MPQRS  = IBTOR(KPQ,KRS)
          KPQRS  = IBTAND(KPQR,MULD)
          HKABCD = FMULT(KPQRS)
          KHABCD = KHKTAB*KHKTCD
          NOABCD = NORBAB*NORBCD
C
C FOR SCF, THE CORE REQUIREMENTS MIGHT BE REDUCED, TO BE
C IMPLEMENTED LATER
C
           IF(SCF.OR.FOCK) THEN
            IWKSO=KSTART
           ELSE
            IPSO  = KSTART
            IPAO = IPSO + KHKTA*MULTA*NORBA
     &                   *KHKTB*MULTB*NORBB
     &                   *KHKTC*MULTC*NORBC
     &                   *KHKTD*MULTD*NORBD
            IWKSO=IPAO+NORBA*NORBB*NORBC*NORBD
           ENDIF
           IF(FOCK.AND.(.NOT.SCF)) THEN
C
C ALLOCATE MEMORY FOR SO INTEGRAL DERIVATIVES
C CORRECT FOR SYMMETRY !)
C
            IGOFF=MULTA*MULTB*MULTC*MULTD
            IWKAO = IWKSO + 6*NOABCD*KHABCD*
     &              MULTA*MULTB*MULTC*MULTD
c     &                 (KHABCD*NOABCD*MULTA*MULTB*MULTC*MULTD)/
c     &                        MAX(MULTA,MULTB,MULTC,MULTD)
C
C FOR C1, NO TRANSFORMATIONS FROM AO TO SAO ARE NEEDED
C
            IF(NIRREP.EQ.1) IWKAO=IWKSO
c            IF(NIRREP.EQ.1) CALL ZERO(WORK1(IWKAO+1),6*NOABCD*KHABCD)
           ELSE
C
C CHECK IF THIS IS CORRECT FOR GRADIENT CALCULATIONS ONLY
C
            IWKAO = IWKSO +
     &                 2*(KHABCD*NOABCD*MULTA*MULTB*MULTC*MULTD)/
     &                        MAX(MULTA,MULTB,MULTC,MULTD)
           ENDIF
C
          INTS = .FALSE.
          CALL ZERO(WORK1(IWKSO+1),IWKAO-IWKSO)
C
C        READ TWO-ELECTRON DENSITY (NOT REQUIRED FOR CORRELATED NMR SHIFTS) 
C
c          IF((.NOT.SCF).AND.(.NOT.FOCK)) THEN
c           READ (LUPAO) LENGTH, (WORK1(IPSO + I),I=1,LENGTH)
c          END IF
C
C    FIRST SYMMETRY INDEX R 
C
C        GENERATES DISTINT OVERLAP DISTRIBUTIOONS A*R(B)
C
          INIT=.TRUE.
          INTWRIT2=.FALSE.
          DO 500 ISYMR = 0,MAXLOT
           IF (IBTAND(ISYMR,MPQ) .EQ. 0) THEN
            ICENTB = NUCNUM(NCENTB,IBTAXO(ISYMR,MULB)+1)
            SIGNBX = XAND(ISYMR)
            SIGNBY = YAND(ISYMR)
            SIGNBZ = ZAND(ISYMR)
C
C   CHARGE DISTRIBUTIONS FOR FIRST ELECTRON
C
            IF (MULTAB .GT. 1) THEN
             TPRIAB=.FALSE.
             CAEQCB = ICENTA .EQ. ICENTB
             CALL GIADRV(NHKTA,NHKTB,NSETA,NSETB,ISTEPA,ISTEPB,
     &                    NUCA,NUCB,NUCAB,NORBA,NORBB,
     &                    NUCOA,NUCOB,NRCOA,NRCOB,JSTRA,JSTRB,
     &                    ONE,ONE,ONE,SIGNBX,SIGNBY,SIGNBZ,
     &                    CORAX,CORAY,CORAZ,CORBX0,CORBY0,CORBZ0,
     &                    AACDX,AACDY,AACDZ,IAB0X,IAB0Y,IAB0Z,TPRIAB,
     &                    CAEQCB,.TRUE.,.TRUE.,BIGVEC,SEGMEN,DTEST,
     &                    12,THRESH,MAXDER,IPRINT)
             XDIFF1=CORAX-SIGNBX*CORBX0
             YDIFF1=CORAY-SIGNBY*CORBY0
             ZDIFF1=CORAZ-SIGNBZ*CORBZ0
             IF (NUCAB .EQ. 0) GO TO 500
            END IF
C
C  SECOND SYMMETRY INDEX S 
C
C  GENERATES DISTINCT OVERLAP DITRIBUTIONS C*S(D)
C
            DO 600 ISYMS = 0, MAXLOT
             IF (IBTAND(ISYMS,MRS) .EQ. 0) THEN
              TPRICD=.FALSE.
C
C  THIRD SYMMETRY INDEX T 
C
C  GENERATES DISTINCT QUADRUPLETS A*R(B) * T(C*S(D))
C
              DO 700 ISYMT = 0, MAXLOT
               IF (IBTAND(ISYMT,MPQRS) .EQ. 0) THEN
                ISYMTS = IBTXOR(ISYMT,ISYMS)
C
C  CHARGE DISTRIBUTIONS FOR SECOND ELECTRON
C
                SIGNCX = XAND(ISYMT)
                SIGNCY = YAND(ISYMT)
                SIGNCZ = ZAND(ISYMT)
                SIGNDX = XAND(ISYMTS)
                SIGNDY = YAND(ISYMTS)
                SIGNDZ = ZAND(ISYMTS)
C
C  ICENTC AND ICENTD
C
                ICENTC = NUCNUM(NCENTC,IBTAXO(ISYMT,MULC)+1)
                ICENTD = NUCNUM(NCENTD,IBTAXO(ISYMTS,MULD)+1)
C
                CCEQCD = ICENTC.EQ.ICENTD
                CALL GIADRV(NHKTC,NHKTD,NSETC,NSETD,ISTEPC,ISTEPD,
     &                       NUCC,NUCD,NUCCD,NORBC,NORBD,
     &                       NUCOC,NUCOD,NRCOC,NRCOD,JSTRC,JSTRD,
     &                       SIGNCX,SIGNCY,SIGNCZ,SIGNDX,SIGNDY,SIGNDZ,
     &                       CORCX0,CORCY0,CORCZ0,CORDX0,CORDY0,CORDZ0,
     &                       ABCCX,ABCCY,ABCCZ,ICD0X,ICD0Y,ICD0Z,
     &                       TPRICD,CCEQCD,.TRUE.,.TRUE.,BIGVEC,
     &                       SEGMEN,DTEST,34,THRESH,MAXDER,IPRINT)
                XDIFF2=CORCX0*SIGNCX-SIGNDX*CORDX0
                YDIFF2=CORCY0*SIGNCY-SIGNDY*CORDY0
                ZDIFF2=CORCZ0*SIGNCZ-SIGNDZ*CORDZ0
                IF (NUCCD .EQ. 0) GO TO 600
C
C   CHECK WHETHER THIS INTEGRAL GIVES ZERO CONTRIBUTION
C   THERE ARE CONTRIBUTIONS ONLY IF EITHER CAEQCB OR CCEQCD
C   IS FALSE, OTHERWISE SKIP AND GO TO 700
C
                 IF (CAEQCB.AND.CCEQCD) GO TO 700
                 INTWRIT2=INTWRIT    
C
C  INTEGRAL CONTRIBUTES
C
C   INTEGRAL DIRECTIVES
C
                 NCCINT = NORBAB*NORBCD*KHKTAB*KHKTCD
                 CALL GIADIR(CAEQCB,CCEQCD,NCCINT,MAXDER,
     &                       IPRINT)
C
C  CALCULATE AO INTEGRALS
C
                 ONECEN=.FALSE.
                 CALL GIADER(IPRINT,NOINT,SEGMEN,WORK1,LWORK1)
                 INTS = INTS .OR. .NOT.NOINT
                 IF(NIRREP.EQ.1) INTWRIT2=INTWRIT.AND.(.NOT.NOINT)
                 IF (NOINT) GO TO 700
C
C  CONSTRUCT FOCK MATRIX DERIVATIVES
C
                  IF(IFIRST) THEN
                   CALL GFOCK(WORK1(KDFAOA),WORK1(KDFAOB),NF,ISYMR,
     &                        ISYMS,ISYMT,
     &                        NORBA,NORBB,NORBC,NORBD,NOABCD,
     &                        ICORBA,ICORBB,ICORBC,ICORBD,
     &                        DIAGAB,DIAGCD,SHABAB,SHAEQB,SHCEQD,
     &                        TCONAB,TCONCD,
     &                        KHKTA,KHKTB,KHKTC,KHKTD,
     &                        NHKTA,NHKTB,NHKTC,NHKTD,
     &                        THRESH,HKABCD,IPRINT,
     &                        WORK1,LWORK1,
     &                        DHFFA,DHFFB,NBAST,NFIELD,NIR)
                   ENDIF
                  IF(.NOT.SCF) THEN
C
C
C  TRANSFORM DIFFERENTIATED INTEGRALS TO SYMMETRY BASIS
C         (ONLY FOR NIRREP.GT.1 REQUIRED)
C
                   IF(NIRREP.NE.1) THEN
C
                    LWRKAO = LWORK1 - IWKAO
                    CALL GSYM2(WORK1(IWKSO+1),WORK1(IWKAO+1),ISYMR,
     &                         ISYMT,ISYMTS,SHABAB,NOABCD,
     &                         HKABCD,IPRINT)
C
                   ENDIF
                  END IF
                 END IF
  700          CONTINUE
            END IF
  600       CONTINUE
         END IF
  500    CONTINUE
C
C DUMP GIAO INTEGRAL DERIVATIVES TO AODER FILE
C
         IF(INTWRIT2) THEN
          CALL GPLUNK(WORK1(IWKSO+1),WORK1(KREORD),
     &                WORK1(KBUCK),WORK1(KIBUCK),NBUF,  
     &                WORK1(KCOUNT),
     &                ISYMR,ISYMS,ISYMT,
     &                NORBA,NORBB,NORBC,NORBD,NOABCD,
     &                ICORBA,ICORBB,ICORBC,ICORBD,
     &                DIAGAB,DIAGCD,
     &                SHAEQB,SHCEQD,SHAEQC,SHBEQD,
     &                TCONAB,TCONCD,
     &                KHKTA,KHKTB,KHKTC,KHKTD,
     &                NHKTA,NHKTB,NHKTC,NHKTD,
     &                MULA,MULB,MULC,MULD,
     &                THRESH,HKABCD,IPRINT,
     &                WORK1,LWORK1,
     &                NBAST,NUMBUF,NIR)
         ENDIF
C
  400 CONTINUE
  300 CONTINUE
  200 CONTINUE
  100 CONTINUE
C
C UNFOLD DERIVATIVES OF THE FOCK MATRIX
C
      IF(IFIRST) THEN
       DO 852 IPERT=1,NFIELD
        IRREP=IXYPERT(IPERT) 
        ICOUNT=0
        DO 853 IPERT1=1,IPERT
         IF(IXYPERT(IPERT1).EQ.IRREP) ICOUNT=ICOUNT+1
853     CONTINUE
         IF(IRREP.EQ.1) THEN
c         IF(.NOT.SCF) THEN 
c          CALL FSYM2A(WORK1(KRFAOA+(IPERT-1)*NF*(MAXLOP+1)), 
c     &               WORK1(KRFAOB+(IPERT-1)*NF*(MAXLOP+1)),
c     &               NBAST,IPRINT,WORK1(KSTART+1),
c     &               WORK1(KSTART+1+NLENT(IRREP)),NIR)
c         CALL PUTLST(WORK1(KSTART+1),
c     &               IPERT-IIPERT(IRREP),1,1,IRREP,105)
c         IF(IUHF.NE.0) THEN
c          CALL PUTLST(WORK1(KSTART+1+NLENT(IRREP)),
c     &                IPERT-IIPERT(IRREP),1,1,IRREP,106)
c         ENDIF
c         ENDIF
         CALL FSYM2A(WORK1(KDFAOA+(IPERT-1)*NF*(MAXLOP+1)),
     &               WORK1(KDFAOB+(IPERT-1)*NF*(MAXLOP+1)),
     &               NBAST,IPRINT,WORK1(KSTART+1),
     &               WORK1(KSTART+1+NLENT(IRREP)),NIR)
          ELSE
c          IF(.NOT.SCF) THEN 
c          CALL FSYM3B(IRREP,WORK1(KRFAOA+(IPERT-1)*NF*(MAXLOP+1)), 
c     &               WORK1(KRFAOB+(IPERT-1)*NF*(MAXLOP+1)),
c     &               NBAST,IPRINT,WORK1(KSTART+1),
c     &               WORK1(KSTART+1+NLENT(IRREP)),NIR,ONEM)
c         CALL PUTLST(WORK1(KSTART+1),IPERT-IIPERT(IRREP),
c     &               1,1,IRREP,105)
c         IF(IUHF.NE.0) THEN
c          CALL PUTLST(WORK1(KSTART+1+NLENT(IRREP)),
c     &                IPERT-IIPERT(IRREP),1,1,IRREP,106)
c         ENDIF
c         ENDIF
           CALL FSYM3B(IRREP,WORK1(KDFAOA+(IPERT-1)*NF*(MAXLOP+1)),
     &                 WORK1(KDFAOB+(IPERT-1)*NF*(MAXLOP+1)),
     &                 NBAST,IPRINT,WORK1(KSTART+1),
     &                 WORK1(KSTART+1+NLENT(IRREP)),NIR,ONEM)
          ENDIF
          IF(IDFGH) THEN
           ISCR=KSTART+1+NLENT(IRREP)*(1+IUHF)
           CALL EXPCAR(IRREP,1,WORK1(KSTART+1),WORK1(ISCR),1)
           IF(IUHF.EQ.1) THEN
            CALL EXPCAR(IRREP,1,WORK1(KSTART+1+NLENT(IRREP)),
     &                  WORK1(ISCR),1)
           ENDIF
          ENDIF
 
        CALL GETLST(WORK1(KDFAOA+(IPERT-1)*NF*(MAXLOP+1)),
     &              ICOUNT,1,1,IRREP,102)
c        do 1232 I=1,nlent(irrep)
c         write(*,*)'old ', i, work1(kdfaoa+(ipert-1)*nf*(maxlop+1)
c     &   -1+i)
c1232    continue
c        do 1231 I=1,nlent(irrep)
c         write(*,*)'add ', i, work1(kstart+i)
c1231    continue
        CALL SAXPY(NLENT5(IRREP),ONE,WORK1(KDFAOA+(IPERT-1)*NF*
     &             (MAXLOP+1)),1,WORK1(KSTART+1),1)
c        do 1233 I=1,nlent(irrep)
c         write(*,*)'final ', i, work1(kstart+i)
c1233    continue
        CALL PUTLST(WORK1(KSTART+1),
     &             ICOUNT,1,1,IRREP,102)
        IF(IUHF.NE.0) THEN
         CALL SAXPY(NLENT5(IRREP),ONE,WORK1(KDFAOA+(IPERT-1)*NF*
     &             (MAXLOP+1)),1,WORK1(KSTART+1+NLENT(IRREP)),1)
         CALL PUTLST(WORK1(KSTART+NLENT(IRREP)+1),ICOUNT,
     &               1,1,IRREP,103)
        ENDIF
852    CONTINUE
      ENDIF
C
      IF(INTWRIT)THEN
C
C FLUSH INTEGRAL DERIVATIVE BUFFERS
C
        CALL BFLUSH(WORK1(KBUCK),WORK1(KIBUCK),
     &             WORK1(KCOUNT),NBUF,NUMBUF)
C
2011  FORMAT(I8,' integral derivatives in ',I5,' records for',
     &       'Bx written to file AODER')
2012  FORMAT(I8,' integral derivatives in ',I5,' records for',
     &       'By written to file AODER')
2013  FORMAT(I8,' integral derivatives in ',I5,' records for',
     &       'Bz written to file AODER')
      ENDIF
      RETURN
 2000 FORMAT (/,1X,' ***** Pair of shells (Electron 1) ',2I3,' *****',/)
 2005 FORMAT (/,1X,' ***** Pair of shells (Electron 2) ',2I3,' *****',/)
 2010 FORMAT (/,1X,' ***************************************** ',
     *        /,1X,' ********** Shells ',4I3,     ' ********** ',
     *        /,1X,' ***************************************** ')
 2020 FORMAT (/,1X,' *********************************************** ',
     *        /,1X,' ********** Sym. operations ',3I3,  ' ********** ',
     *        /,1X,' *********************************************** ')
 1000 FORMAT(//,'  ICENT-D   ',4I5)
 1005 FORMAT(//,'  NHKTA-D   ',4I5)
 1010 FORMAT(   '  JMAX0     ',I1,' + ',I1,' = ',I2)
 1020 FORMAT(   '  KHKTA-D   ',4I5)
 1030 FORMAT(   '  KHKTAB-CD ',2I5)
 1040 FORMAT(   '  KHABCD    ',I5)
 1045 FORMAT(   '  MHKTA-D   ',4I5)
 1050 FORMAT(   '  DIAGAB    ',L5,3I5)
 1060 FORMAT(   '  DIAGCD    ',L5,3I5)
 1070 FORMAT(   '  DIACAC    ',L5)
 1075 FORMAT(   '  TPRI/CONAB',2L5)
 1076 FORMAT(   '  TPRI/CONCD',2L5)
 1080 FORMAT(   '  PQSYM     ',L5,3I5)
 1090 FORMAT(   '  NORBA-D   ',4I5)
 1100 FORMAT(   '  NORBAB-CD ',2I5)
 1110 FORMAT(   '  NUCAB-CD  ',2I5)
 1120 FORMAT(   '  ISTEPA-D  ',4I5)
 1130 FORMAT(   '  NSTRA-D   ',4I5)
   50 FORMAT(/,1X,' ***** Timings for TWOINT ***** ',/)
   51 FORMAT(1X, ' Time in ODCVEC: ',F8.2,' seconds')
   52 FORMAT(1X, ' Time in EXCOEF: ',F8.2,' seconds')
   62 FORMAT(1X, ' Time in PATH:   ',F8.2,' seconds')
   59 FORMAT(1X, ' Time in R000:   ',F8.2,' seconds')
   54 FORMAT(1X, ' Time in HERI:   ',F8.2,' seconds')
   55 FORMAT(1X, ' Time in C1DRIV: ',F8.2,' seconds')
   61 FORMAT(1X, ' Time in C2DRIV: ',F8.2,' seconds')
   64 FORMAT(1X, ' Time in SYM2:   ',F8.2,' seconds')
   56 FORMAT(1X, ' Time in DISINT: ',F8.2,' seconds',/)
   57 FORMAT(1X,' Detailed timings for R000 integrals ',
     *       '(JMAX = 1, ',I2,'):',/)
   58 FORMAT(1X,10F6.2)
   60 FORMAT(/,1X,' Detailed timings for Hermitian integrals ',
     *       '(JMAX = 1, ',I2,'):',/)
   63 FORMAT(1X)
   70 FORMAT(/,1X,' ***** Maximum WORK1 space allocations ***** ',/)
   71 FORMAT(1X, ' Undifferentiated SSSS integrals         : ',I7)
   72 FORMAT(1X, ' Hermitian integrals                     : ',I7)
   73 FORMAT(1X, ' Cartesian-Hermitian integrals (special) : ',I7)
   74 FORMAT(1X, ' Cartesian-Hermitian integrals (general) : ',I7)
   75 FORMAT(1X, ' Cartesian integrals (C2HINT)            : ',I7)
   76 FORMAT(1X, ' Cartesian integrals (C2EINT)            : ',I7)
   80 FORMAT(/,1X,' ***** Maximum WORK3 space allocations ***** ',/)
   81 FORMAT(1X, ' DPQX/Y/Z : ',I7)
   82 FORMAT(1X, ' RJ000    : ',I7)
      END
