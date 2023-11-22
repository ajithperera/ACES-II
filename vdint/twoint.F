      SUBROUTINE TWOINT(IRNTYP,MAXDIF,FOCK,JATOM,NODV,NOPV,
     &                  ITIME,JPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,
     &                  RETUR,WORK1,LWORK1,ENEREE,SCF_GRADEE,GRADEE,
     &                  SCF_HESSEE,HESSEE,
     &                  NCOORD,DHFA,DHFB,DRELA,DRELB,DNHFA,DNHFB,
     &                  DHFFA,DHFFB,DRELFA,DRELFB,NBAST,NIR)
C
C     References for calculation of Cartesian integrals:
C
C     L. E. McMurchie & E. R. Davidson, J. Comput. Chem. 26, 218 (1978)
C     V. R. Saunders, in "Methods in Computational Molecular Physics",
C       G. H. F. Diercksen and S. Wilson, eds. (Reidel,Dordrecht,1983)
C     T. U. Helgaker et al., JCP 84, 6266 (1986)
C
C     References for symmetry:
C
C     E. R. Davidson, JCP 62, 400 (1975)
C     P. R. Taylor,   TCA 69, 447 (1986)
C     J. Almloef, MOLECULE Program Description, USIP Report 74 - 29
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
      CHARACTER NAMEX*6
      LOGICAL DCORD, DCORGD, NOORBT, DOPERT
      LOGICAL SHAEQB, SHCEQD, SHAEQC, SHBEQD, SHAAAA, SHABAB,
     *        AACDX, AAAAX, AACDY, AAAAY, AACDZ, AAAAZ,
     *        ABCCX, ABCCY, ABCCZ, DERECC, DERECD,
     *        PRINT1, PRINT2, PRINT3, PRINT4,NOPV,NOPMAT, NODV,
     *        PERTUR, EXPECT, FIRST, NOINT, INTS, UNDIFF,
     *        CAEQCB, CCEQCD, RETUR, DIFX, DIFY, DIFZ, TTIME,
     *        FOCK,GRAD,INTWRIT,TRAINV,IDFGH
      LOGICAL SHARE
      LOGICAL BIGVEC, SEGMEN
      LOGICAL         AEQB,   CEQD,   DIAGAB, DIAGCD, DIACAC,
     *                ONECEN, PQSYM,  DTEST,
     *                TPRIAB, TPRICD, TCONAB, TCONCD
      LOGICAL         DZER
      LOGICAL SCF,NONHF,JUNK
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "mxcent.par"
#include "baslims.par"
      PARAMETER (AZERO = 0.00 D00, ONE = 1.00 D00, LUPAO = 25)
 
#include "icdacc.com"

      DIMENSION XAND(0:7),YAND(0:7),ZAND(0:7)
      DIMENSION NUCOA(MXAOVC),  NUCOB(MXAOVC),
     *          NUCOC(MXAOVC),  NUCOD(MXAOVC),
     *          NRCOA(MXCONT),  NRCOB(MXCONT),
     *          NRCOC(MXCONT),  NRCOD(MXCONT),
     *          JSTRA(MXAOVC),  JSTRB(MXAOVC),
     *          JSTRC(MXAOVC),  JSTRD(MXAOVC),
     *          ICORBA(MXAOVC), ICORBB(MXAOVC),
     *          ICORBC(MXAOVC), ICORBD(MXAOVC)
      DIMENSION DHFA(NBAST,NBAST,0:NIR),
     &          DHFB(NBAST,NBAST,0:NIR),
     &          DRELA(NBAST,NBAST,0:NIR),
     &          DRELB(NBAST,NBAST,0:NIR),
     &          DNHFA(NBAST,NBAST,0:NIR),
     &          DNHFB(NBAST,NBAST,0:NIR)
      dimension ichain(600),ninbck(600)
      DIMENSION DHFFA(NBAST*(NBAST+1)/2,0:NIR),
     &          DHFFB(NBAST*(NBAST+1)/2,0:NIR),
     &          DRELFA(NBAST*(NBAST+1)/2,0:NIR),
     &          DRELFB(NBAST*(NBAST+1)/2,0:NIR)
      DIMENSION JPERT(8),KPERT(8)
      DIMENSION  WORK1(LWORK1)
      DIMENSION GRADEE(NCOORD),HESSEE(NCOORD,NCOORD),
     &          SCF_GRADEE(NCOORD),SCF_HESSEE(NCOORD,NCOORD)
C
      COMMON /CCOM/ THRS, NHTYP, IBX
      COMMON /SHELLSi/ KMAX,
     *                NHKT(MXSHEL),   KHKT(MXSHEL), MHKT(MXSHEL),
     *                ISTBAO(MXSHEL), NUCO(MXSHEL), JSTRT(MXSHEL),
     *                NSTRT(MXSHEL),  MST(MXSHEL),  NCENT(MXSHEL),
     *                NRCO(MXSHEL), NUMCF(MXSHEL),
     *                NBCH(MXSHEL),   KSTRT(MXSHEL)
      COMMON /SHELLS/ CENT(MXSHEL,3),
     *                SHARE(MXSHEL)
      COMMON /BLOCKS/ CENTSH(MXSHEL,3),
     *                MAXSHL, BIGVEC, SEGMEN,
     *                NHKTSH(MXSHEL), KHKTSH(MXSHEL), MHKTSH(MXSHEL),
     *                ISTBSH(MXSHEL), NUCOSH(MXSHEL), NORBSH(MXSHEL),
     *                NSTRSH(MXSHEL), NCNTSH(MXSHEL), NSETSH(MXSHEL),
     *                JSTRSH(MXSHEL,MXAOVC),
     *                NPRIMS(MXSHEL,MXAOVC),
     *                NCONTS(MXSHEL,MXAOVC),
     *                IORBSH(MXSHEL,MXAOVC),
     *                IORBSB(MXCORB), NRCSH(MXSHEL)
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
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
      COMMON /NUCLEC/ NAMEX(MXCOOR)
      COMMON /INTADR/ IWKAO, IWKSO, IWKHHS, IWK1HH, IWK1HC, IWKLST
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
      COMMON /DERZER/ FZERO, DZER, JZERO
      COMMON/METHOD/IUHF,SCF,NONHF
      COMMON/BASSYM/NBAS(8),NBASIS,NBASSQ,NBASTT
      COMMON/BASSPH/NBAS5(8),NBASIS5,NBASSQ5,NBASTT5
      COMMON/LSYM/NLENQ(8),NLENT(8)
      COMMON/LSYMSPH/NLENQ5(8),NLENT5(8)
      COMMON/OPTION/GRAD,INTWRIT,JUNK
      COMMON/PERT/NTPERT,NPERT(8),IIPERT(8),IXPERT,IYPERT,IZPERT,
     &            IYZPERT,IXZPERT,IXYPERT,ITRANSX,ITRANSY,ITRANSZ,
     &            NUCIND1
      COMMON/INVAR/TRAINV
      COMMON/MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD

      COMMON/IPRI/IPRINT
      COMMON/DFGH/IDFGH
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
      PTAND(I,J) = PT(IBTAND(I,J))
      IBTAXO(I,J) = IBTAND(I,IBTXOR(I,J))
C
#ifdef _ASSERT
      if (LWORK1.LT.1) then
         print *, '@TWOINT: Assertion failed.'
         print *, '         LWORK1 = ',LWORK1
         call errex
      end if
#endif
      FIRST  = .TRUE.
      NOPMAT = NOPV
C
C  DETERMINE RUN TYPE
C
      DTEST  = .FALSE.
      PERTUR = .FALSE.
      UNDIFF = .FALSE.
      EXPECT = .FALSE.
C
C     Integrals differentiated with respect to JATOM only
C     (for use in MCSCF response equations)
C
C     Note: if MAXDIF = 0 undifferentiated integrals are written on LU2D
C
c      ELSE IF (IRNTYP .EQ. 1) THEN
c         PERTUR = .TRUE.
c         MAXDER = MAXDIF
c         IATOM  = JATOM
c         MULE   = ISTBNU(IATOM)
c         MULTE  = MULT(MULE)
C
C     Expectation values calculated (all atoms)
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
C ADY:
C This may be run during a Kohn-Sham calculation, but the seamless nature
C of KS and HF has not been exploited yet; therefore, we need to check
C for the HF exchange potential's coefficient in two steps:
C    1) see if the record PCHF exists, then
C    2) if it doesn't, this is a plain old HF calc and set PCHF=1.0
C       otherwise use it.
C
         IONE=IINTFP
         CALL GETREC(0,'JOBARC','PCHF    ',IONE,PCHF)
         IF (IONE.EQ.-1) THEN
            PCHF=1.0
         ELSE
            CALL GETREC(1,'JOBARC','PCHF    ',IINTFP,PCHF)
         END IF
c         write(*,*)
c         write(*,*) "@VDINT: The HF exchange coefficient is ",PCHF
c         write(*,*)
C
C     Rewind LUPAO
C
      IF((.NOT.SCF).AND. (.NOT.FOCK)) THEN
       REWIND LUPAO
      END IF
C
C ALLOCATE CORE MEMORY FOR DERIVATIVES OF FOCK MATRIX 
C
      KSTART=0
      IF(FOCK) THEN
C
       KDFAOA=KSTART+1
       NF=NBAST*(NBAST+1)/2
       KDFAOB=KDFAOA+NF*NCOORD*IUHF*(NIR+1)
       KSTART=KDFAOB+NF*NCOORD*(NIR+1)-1
       if (kstart.gt.lwork1) then
          print *, '@TWOINT: Insufficient memory.'
          print *, '         need ',kstart,' dbls'
          print *, '         have ',lwork1,' dbls'
          call errex
       end if
       CALL ZERO(WORK1(KDFAOA),(IUHF+1)*NF*NCOORD*(NIR+1))
       IF(.NOT.SCF) THEN
        nbuf=600
        irecl=(1+iintfp)*600+2
        KRFAOA=KSTART+1
        KRFAOB=KRFAOA+NF*NCOORD*IUHF*(NIR+1)
        KBUCK=KRFAOB+NF*NCOORD*(NIR+1)
        KIBUCK=KBUCK+NCOORD*NBUF
        KCOUNT=KIBUCK+NCOORD*NBUF
        KCHAIN=KCOUNT+NCOORD*NBUF
        KMAP=KCHAIN+NCOORD*NBUF
        KSTART=KMAP+NCOORD
        if (kstart.gt.lwork1) then
           print *, '@TWOINT: Insufficient memory.'
           print *, '         need ',kstart,' dbls'
           print *, '         have ',lwork1,' dbls'
           call errex
        end if
        CALL MAKMAP(WORK1(KMAP),NCOORD,IXPERT,IYPERT,IZPERT,
     &              ITRANSX,ITRANSY,ITRANSZ,TRAINV)       
        CALL ZERO(WORK1(KRFAOA),(IUHF+1)*NF*NCOORD*(NIR+1))
        CALL IZERO(WORK1(KCOUNT),NCOORD)
        CALL IZERO(WORK1(KCHAIN),NCOORD)
C 
C OPEN SORT FILE FOR INTEGRAL DERIVATIVES
C
c YAU - Why is this IRECL doubled? Many other routines set RECL to
c (600*(1+iintfp)+2)*idaccm
        OPEN(UNIT=30,FILE='AODER',FORM='UNFORMATTED',STATUS='UNKNOWN',
     &       ACCESS='DIRECT',RECL=IRECL*2*idaccm)
C
C RESET NUMBER OF RECORDS TO ZERO
C
        NREC=1
C
       ENDIF
      ENDIF
C
C ZERO ENERGY, GRADIENT AND HESSIAN CONTRIBUTIONS
C
      ENEREE=AZERO
C
      IF(GRAD) THEN
C
       CALL ZERO(SCF_GRADEE,NCOORD)
       CALL ZERO(GRADEE,NCOORD)
C
       IF(MAXDER.EQ.2) THEN
C
        CALL ZERO(SCF_HESSEE,NCOORD*NCOORD)
        CALL ZERO(HESSEE,NCOORD*NCOORD)
C
       ENDIF
C
      ENDIF        
C
C First Shell Index 
C
      DO 100 ISHELA = 1,MAXSHL
CSSS      DO 100 ISHELA = 9,9
CSSS       Write(6,*) "I am HERE1", ISHELA 
       NHKTA = NHKTSH(ISHELA)
       KHKTA = KHKTSH(ISHELA)
       MHKTA = MHKTSH(ISHELA)
       NCENTA = NCNTSH(ISHELA)
       IF (NCENTA .GT. 0) THEN
        ICENTA = NUCNUM(NCENTA,1)
       ELSE
        ICENTA = -1
       END IF
       IF (.NOT.BIGVEC) THEN
        CORAX = CENTSH(ISHELA,1)
        CORAY = CENTSH(ISHELA,2)
        CORAZ = CENTSH(ISHELA,3)
       END IF
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
       IF ((ISHELA .NE. IPRNTA).AND.(IPRNTA .NE. 0)) PRINT1 = .FALSE.
       SIGNAX = ONE
       SIGNAY = ONE
       SIGNAZ = ONE
C
C Second Shell Index 
C
       DO 200 ISHELB = 1,ISHELA
CSSS       DO 200 ISHELB = 1,1
CSSS         Write(6,*) "I AM HERE2", ISHELB
        SHAEQB = ISHELA .EQ. ISHELB
        NHKTB  = NHKTSH(ISHELB)
        KHKTB  = KHKTSH(ISHELB)
        MHKTB  = MHKTSH(ISHELB)
        NCENTB = NCNTSH(ISHELB)
        IF (.NOT.BIGVEC) THEN
         CORBX0 = CENTSH(ISHELB,1)
         CORBY0 = CENTSH(ISHELB,2)
         CORBZ0 = CENTSH(ISHELB,3)
        END IF
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
C Charge Distributions for First Electron 
C
        AEQB = SHAEQB .AND. MULTA .EQ. 1
        TPRIAB = AEQB
        tpriab=.false.
        TCONAB = SHAEQB .AND. MAXAB .EQ. 0
        IF (TCONAB) THEN
         NORBAB = NORBA*(NORBA + 1)/2
        ELSE
         NORBAB = NORBA*NORBB
        END IF
        IF (MULTAB .EQ. 1) THEN
         CAEQCB = .NOT.BIGVEC .AND. ICENTA.EQ.NUCNUM(NCENTB,1)
         CALL ODCDRV(NHKTA,NHKTB,NSETA,NSETB,ISTEPA,ISTEPB,
     *               NUCA,NUCB,NUCAB,NORBA,NORBB,
     *               NUCOA,NUCOB,NRCOA,NRCOB,JSTRA,JSTRB,
     *               ONE,ONE,ONE,ONE,ONE,ONE,
     *               CORAX,CORAY,CORAZ,CORBX0,CORBY0,CORBZ0,
     *               AACDX,AACDY,AACDZ,IAB0X,IAB0Y,IAB0Z,TPRIAB,
     *               CAEQCB,.TRUE.,.TRUE.,BIGVEC,SEGMEN,DTEST,
     *               12,THRESH,MAXDER,IPRINT)
         IF (NUCAB .EQ. 0) GO TO 200
        END IF
        DIAGAB = SHAEQB .AND. .NOT.BIGVEC
        IF (DIAGAB) THEN
         KHKTAB = KHKTA*(KHKTA + 1)/2
        ELSE
         KHKTAB = KHKTA*KHKTB
        END IF
C
C Third Shell Index 
C
        DO 300 ISHELC = 1,ISHELA
CSSS        DO 300 ISHELC = 9,9
CSSS         Write(6,*) "I AM HERE3", ISHELC
         SHAEQC = ISHELA .EQ. ISHELC
         NHKTC  = NHKTSH(ISHELC)
         NCENTC = NCNTSH(ISHELC)
         IF (.NOT.BIGVEC) THEN
          CORCX0 = CENTSH(ISHELC,1)
          CORCY0 = CENTSH(ISHELC,2)
          CORCZ0 = CENTSH(ISHELC,3)
         END IF
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
C Fourth Shell Index 
C
         IDMAX = ISHELC
CSSS         IDMAX = 4

         IF (SHAEQC) IDMAX = ISHELB
         DO 400 ISHELD = 1,IDMAX
CSSS         DO 400 ISHELD = 4,4
CSSS         Write(6,*) "I AM HERE4", ISHELD
          SHBEQD  = ISHELB .EQ. ISHELD
          SHCEQD  = ISHELC .EQ. ISHELD
          NCENTD = NCNTSH(ISHELD)
          IF (.NOT.BIGVEC) THEN
           CORDX0 = CENTSH(ISHELD,1)
           CORDY0 = CENTSH(ISHELD,2)
           CORDZ0 = CENTSH(ISHELD,3)
          END IF
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
c          IF (GRAD) THEN
C
C FOR SCF, THE CORE REQUIREMENTS MIGHT BE REDUCED, TO BE
C IMPLEMENTED LATER
C
           IF(SCF.OR.FOCK) THEN
            IWKSO=KSTART
           ELSE
            IPSO  = KSTART
            IPAO = IPSO + KHKTA*MULTA*NORBA
     *                   *KHKTB*MULTB*NORBB
     *                   *KHKTC*MULTC*NORBC
     *                   *KHKTD*MULTD*NORBD
            IWKSO=IPAO+NORBA*NORBB*NORBC*NORBD

           ENDIF
           IF(FOCK.AND.(.NOT.SCF)) THEN
C
C ALLOCATE MEMORY FOR SO INTEGRAL DERIVATIVES
C CORRECT FOR SYMMETRY !)
C
            IWKAO = IWKSO + 12*NOABCD*KHABCD
c     *                 (KHABCD*NOABCD*MULTA*MULTB*MULTC*MULTD)/
c     *                        MAX(MULTA,MULTB,MULTC,MULTD)
           ELSE
C
C CHECK IF THIS IS CORRECT FOR GRADIENT CALCULATIONS ONLY
C
            IWKAO = IWKSO +
     *                 (KHABCD*NOABCD*MULTA*MULTB*MULTC*MULTD)/
     *                        MAX(MULTA,MULTB,MULTC,MULTD)
           ENDIF
c          ELSE
c           IWKSO = KSTART
c           IWKAO = IWKSO
c     *           + 3*MULTE*(KHABCD*NOABCD*MULTA*MULTB*MULTC*MULTD)/
c     *                      MAX(MULTA,MULTB,MULTC,MULTD)
c          END IF
C
          INTS = .FALSE.
          if (iwkao.gt.lwork1) then
             print *, '@TWOINT: Insufficient memory.'
             print *, '         need ',iwkao,' dbls'
             print *, '         have ',lwork1,' dbls'
             call errex
          end if
          CALL ZERO(WORK1(IWKSO+1),IWKAO-IWKSO)
C
C        Read two-electron density
C
          IF((.NOT.SCF).AND.(.NOT.FOCK)) THEN
           READ (LUPAO) LENGTH, (WORK1(IPSO + I),I=1,LENGTH)
          END IF
CSSS          Write(6,*)  "The lenth of TD block", LENGTH
CSSS          Write(6, "(6(1x,F12.6))") (WORK1(IPSO+i),
CSSS     &                                       i=1,length)
CSSS          Write(6,*)
C
C First Symmetry Index R 
C
C Generates distinct overlap distributions A*R(B)
C
          DO 500 ISYMR = 0,MAXLOT
           IF (IBTAND(ISYMR,MPQ) .EQ. 0) THEN
            IF (NCENTB .GT. 0) THEN
             ICENTB = NUCNUM(NCENTB,IBTAXO(ISYMR,MULB)+1)
            ELSE
             ICENTB = -2
            END IF
            SIGNBX = XAND(ISYMR)
            SIGNBY = YAND(ISYMR)
            SIGNBZ = ZAND(ISYMR)
C
C Charge Distributions for First Electron 
C
            IF (MULTAB .GT. 1) THEN
             TPRIAB = SHAEQB .AND. ISYMR .EQ. 0
             tpriab=.false.
             CAEQCB = .NOT.BIGVEC .AND. ICENTA .EQ. ICENTB
     *                            .AND. IATOM  .NE. NCENTA
             CALL ODCDRV(NHKTA,NHKTB,NSETA,NSETB,ISTEPA,ISTEPB,
     *                   NUCA,NUCB,NUCAB,NORBA,NORBB,
     *                   NUCOA,NUCOB,NRCOA,NRCOB,JSTRA,JSTRB,
     *                   ONE,ONE,ONE,SIGNBX,SIGNBY,SIGNBZ,
     *                   CORAX,CORAY,CORAZ,CORBX0,CORBY0,CORBZ0,
     *                   AACDX,AACDY,AACDZ,IAB0X,IAB0Y,IAB0Z,TPRIAB,
     *                   CAEQCB,.TRUE.,.TRUE.,BIGVEC,SEGMEN,DTEST,
     *                   12,THRESH,MAXDER,IPRINT)
             IF (NUCAB .EQ. 0) GO TO 500
            END IF
C
C Second Symmetry Index S 
C
C Generates distinct overlap distributions C*S(D)
C
            DO 600 ISYMS = 0, MAXLOT
             IF (IBTAND(ISYMS,MRS) .EQ. 0) THEN
              TPRICD = SHCEQD .AND. ISYMS .EQ. 0
              tpricd=.false.
C
C Third Symmetry Index T 
C
C Generates distinct quadruplets A*R(B) * T(C*S(D))
C
              DO 700 ISYMT = 0, MAXLOT
               IF (IBTAND(ISYMT,MPQRS) .EQ. 0) THEN
                ISYMTS = IBTXOR(ISYMT,ISYMS)
C
C Charge Distributions for Second Electron 
C
                SIGNCX = XAND(ISYMT)
                SIGNCY = YAND(ISYMT)
                SIGNCZ = ZAND(ISYMT)
                SIGNDX = XAND(ISYMTS)
                SIGNDY = YAND(ISYMTS)
                SIGNDZ = ZAND(ISYMTS)
c The ALTIX machines equipped with EPIC compilers do not like to execute
c statement functions in logical assignments.  Ajith Perera, 05/2005.
CSSS                CCEQCD = .NOT.BIGVEC .AND.
CSSS     *                   NUCNUM(NCENTC,IBTAXO(ISYMT ,MULC)+1) .EQ.
CSSS     *                   NUCNUM(NCENTD,IBTAXO(ISYMTS,MULD)+1) .AND.
CSSS     *                   NCENTC .NE. IATOM
                ITMPC  = NUCNUM(NCENTC,IBTAXO(ISYMT ,MULC)+1)
                ITMPD  = NUCNUM(NCENTD,IBTAXO(ISYMTS,MULD)+1)
                CCEQCD = .NOT.BIGVEC .AND. ITMPC.EQ.ITMPD .AND.
     *                   NCENTC .NE. IATOM
                CALL ODCDRV(NHKTC,NHKTD,NSETC,NSETD,ISTEPC,ISTEPD,
     *                      NUCC,NUCD,NUCCD,NORBC,NORBD,
     *                      NUCOC,NUCOD,NRCOC,NRCOD,JSTRC,JSTRD,
     *                      SIGNCX,SIGNCY,SIGNCZ,SIGNDX,SIGNDY,SIGNDZ,
     *                      CORCX0,CORCY0,CORCZ0,CORDX0,CORDY0,CORDZ0,
     *                      ABCCX,ABCCY,ABCCZ,ICD0X,ICD0Y,ICD0Z,
     *                      TPRICD,CCEQCD,.TRUE.,.TRUE.,BIGVEC,
     *                      SEGMEN,DTEST,34,THRESH,MAXDER,IPRINT)
                IF (NUCCD .EQ. 0) GO TO 600
C
C                ICENTC and ICENTD
C
                 IF (NCENTC .GT. 0) THEN
                  ICENTC = NUCNUM(NCENTC,IBTAXO(ISYMT,MULC)+1)
                 ELSE
                  ICENTC = -3
                 END IF
                 IF (NCENTD .GT. 0) THEN
                  ICENTD = NUCNUM(NCENTD,IBTAXO(ISYMTS,MULD)+1)
                 ELSE
                  ICENTD = -4
                 END IF
C
C Check whether this integral gives zero contribution
C
                 IF (PERTUR .AND. MAXDER .GT. 0) THEN
                  IF (IATOM .NE. ICENTA .AND. IATOM .NE. ICENTB
     *                                  .AND. IATOM .NE. ICENTC
     *                                  .AND. IATOM .NE. ICENTD)
     *                GO TO 700
                  END IF
                  ONECEN = (ICENTA .EQ. ICENTB) .AND.
     *                     (ICENTA .EQ. ICENTC) .AND.
     *                     (ICENTA .EQ. ICENTD) .AND.
     *                     (ICENTA .NE. 0)
                  IF (ONECEN) THEN
                   IF (MAXDER .EQ. 0) THEN
                    IF (IBTAND(JMAX0,1).EQ.1)  GO TO 700
c                    IF (DTEST .AND. .NOT.DZER) GO TO 700
                   ELSE
                    IF (PERTUR .OR. EXPECT)    GO TO 700
c                    IF (DTEST .AND. .NOT.DZER) GO TO 700
                   END IF
                  END IF
C
C Integral contributes
C
                  IF (BIGVEC) THEN
                   CALL CORDIF(NSETA,NSETC,THRESH,DIFX,DIFY,DIFZ,
     *                         IPRINT,NUCOA,NUCOC,JSTRA,JSTRC)
                  ELSE
                   DIFX  = ABS(CORAX - XAND(ISYMT)*CORCX0) .LT. THRESH
                   DIFY  = ABS(CORAY - YAND(ISYMT)*CORCY0) .LT. THRESH
                   DIFZ  = ABS(CORAZ - ZAND(ISYMT)*CORCZ0) .LT. THRESH
                  END IF
                  AAAAX = AACDX .AND. DIFX .AND. ABCCX
                  AAAAY = AACDY .AND. DIFY .AND. ABCCY
                  AAAAZ = AACDZ .AND. DIFZ .AND. ABCCZ
C
                  PQSYM = AAAAX .OR. AAAAY .OR. AAAAZ
                  IPQ0X = 0
                  IPQ0Y = 0
                  IPQ0Z = 0
                  IF (AAAAX) IPQ0X = 1
                  IF (AAAAY) IPQ0Y = 1
                  IF (AAAAZ) IPQ0Z = 1
                  DIACAC = ((SHAEQC) .AND. (0 .EQ. ISYMT))
     *               .AND. ((ISHELB .EQ. ISHELD).AND.(ISYMR.EQ.ISYMTS))
C
                  IF (IPRINT .GE. 5) THEN
                   WRITE (LUPRI, 2020) ISYMR, ISYMS, ISYMT
                   WRITE (LUPRI, '(A,I5)') ' ISYMTS ', ISYMTS
                  END IF
C
C Integral Directives 
C
                   NCCINT = NORBAB*NORBCD*KHKTAB*KHKTCD
                   CALL DIRECT(BIGVEC,ICENTA,ICENTB,ICENTC,ICENTD,
     *                         NCENTA,NCENTB,NCENTC,NCENTD,
     *                         0,     ISYMR, ISYMT, ISYMTS,
     *                         SIGNAX,SIGNAY,SIGNAZ,
     *                         SIGNBX,SIGNBY,SIGNBZ,
     *                         SIGNCX,SIGNCY,SIGNCZ,
     *                         SIGNDX,SIGNDY,SIGNDZ,
     *                         NCCINT,MAXDER,DERECC,DERECD,
     *                         EXPECT,IATOM,MULTE,IPRINT)
                    DERECC = .TRUE.
                    DERECD = .TRUE.
C
C Calculate AO Integrals 
CSSS                   Write(6,"(4(1x,i4),1x,i8))") KHKTA*MULTA*NORBA, 
CSSS     *                                    KHKTB*MULTB*NORBB,
CSSS     *                                    KHKTC*MULTC*NORBC,
CSSS     *                                    KHKTD*MULTD*NORBD,IWKAO-IWKSO
C
CSSS                  Write(6,*) "Entering intder"
                  CALL INTDER(IPRINT,NOINT,SEGMEN,WORK1,LWORK1)
   
CSSS                  Write(6,*)
CSSS                  Write(6, "(6(1x,F12.6))") (WORK1(IWKAO+i),
CSSS     &                                       i=1,3*(IWKAO-IWKSO))
CSSS                  Write(6,*)

                  INTS = INTS .OR. .NOT.NOINT
                  IF (NOINT) GO TO 700

C
C Transform to symmetry basis 
C
C
C Expectation value of differentiated integrals
C
                   IF(SCF.OR.(.NOT.FOCK)) THEN 
                     CALL INTEXP(WORK1(IPSO+1),WORK1(IPAO+1),
     *                           ISYMR,ISYMS,ISYMT,
     *                           NORBA,NORBB,NORBC,NORBD,NOABCD,
     *                           ICORBA,ICORBB,ICORBC,ICORBD,
     *                           DIAGAB,DIAGCD,SHABAB,SHAEQB,SHCEQD,
     *                           TCONAB,TCONCD,
     *                           KHKTA,KHKTB,KHKTC,KHKTD,
     *                           NHKTA,NHKTB,NHKTC,NHKTD,
     *                           THRESH,HKABCD,IPRINT,
     *                           BIGVEC,NOPMAT,NODV,WORK1,LWORK1,
     &                           ENEREE,SCF_GRADEE,GRADEE,SCF_HESSEE,
     &                           HESSEE,NCOORD,
     &                           DHFA,DHFB,DRELA,DRELB,DNHFA,
     &                           DNHFB,NBAST,NIR,PCHF)
                   ENDIF
                   IF(FOCK) THEN
                    CALL DFOCK(WORK1(KDFAOA),WORK1(KDFAOB),NF,ISYMR,
     &                         ISYMS,ISYMT,
     &                         NORBA,NORBB,NORBC,NORBD,NOABCD,
     &                         ICORBA,ICORBB,ICORBC,ICORBD,
     &                         DIAGAB,DIAGCD,SHABAB,SHAEQB,SHCEQD,
     &                         TCONAB,TCONCD,
     &                         KHKTA,KHKTB,KHKTC,KHKTD,
     &                         NHKTA,NHKTB,NHKTC,NHKTD,
     &                         THRESH,HKABCD,IPRINT,
     &                         WORK1,LWORK1,
     &                         DHFFA,DHFFB,NBAST,NCOORD,NIR)
                    IF(.NOT.SCF) THEN
                     CALL DFOCK(WORK1(KRFAOA),WORK1(KRFAOB),NF,
     &                          ISYMR,ISYMS,ISYMT,
     &                          NORBA,NORBB,NORBC,NORBD,NOABCD,
     &                          ICORBA,ICORBB,ICORBC,ICORBD,
     &                          DIAGAB,DIAGCD,SHABAB,SHAEQB,SHCEQD,
     &                          TCONAB,TCONCD,
     &                          KHKTA,KHKTB,KHKTC,KHKTD,
     &                          NHKTA,NHKTB,NHKTC,NHKTD,
     &                          THRESH,HKABCD,IPRINT,
     &                          WORK1,LWORK1,
     &                          DRELFA,DRELFB,NBAST,NCOORD,NIR)
C
C Differentiated integrals
C
                     LWRKAO = LWORK1-IWKAO
                     CALL DRSYM2(WORK1(IWKSO+1),WORK1(IWKAO+1),HKABCD,
     &                           ISYMR,ISYMT,ISYMTS,MULTE,MULE,SHABAB,
     &                           NOABCD,LWRKAO,IPRINT)
C
C ADD HERE INTEGRAL DERIVATIVES TO SO INTEGRAL DERIVATIVES
C (HAS TO BE MODIFIED FOR SYMMETRY LATER)
C
c                     CALL DRSYM3(WORK1(IWKSO+1),WORK1(IWKAO+1),HKABCD,
c     *                           ISYMR,ISYMT,ISYMTS)
                    ENDIF
                  END IF
               END IF
  700          CONTINUE
            END IF
  600       CONTINUE
         END IF
  500    CONTINUE
C
         IF(INTWRIT.AND..NOT.ONECEN) THEN
c         IF(.NOT.ONECEN) THEN
          CALL DPLUNK(WORK1(iwkao+1),ISYMR,ISYMS,ISYMT,
     &                NORBA,NORBB,NORBC,NORBD,NOABCD,
     &                ICORBA,ICORBB,ICORBC,ICORBD,
     &                DIAGAB,DIAGCD,TCONAB,TCONCD,
     &                KHKTA,KHKTB,KHKTC,KHKTD,
     &                NHKTA,NHKTB,NHKTC,NHKTD,
     &                THRESH,HKABCD,IPRINT,
     &                WORK1,LWORK1,
     &                NBAST,NCOORD,NIR,
     &                WORK1(KBUCK),WORK1(KIBUCK),NBUF,  
     &                WORK1(KCOUNT),WORK1(KCHAIN),NREC,
     &                WORK1(KMAP))
         ENDIF
C
C Process Integrals 
C
         IF (.NOT. EXPECT .AND. INTS) THEN
C
C Differentiated integrals
C
C x direction
C
               CALL DSOUT2(WORK1(IWKSO+1),NOABCD,IATOM,1,MULTE,MULE,
     *                     SHAEQB,SHCEQD,SHABAB,FIRST,.FALSE.,IPRINT)
               FIRST = .FALSE.
C
C y direction
C
               CALL DSOUT2(WORK1(IWKSO+1),NOABCD,IATOM,2,MULTE,MULE,
     *                     SHAEQB,SHCEQD,SHABAB,FIRST,.FALSE.,IPRINT)
C
C z direction
C
               CALL DSOUT2(WORK1(IWKSO+1),NOABCD,IATOM,3,MULTE,MULE,
     *                     SHAEQB,SHCEQD,SHABAB,FIRST,.FALSE.,IPRINT)
         END IF
         IF (RETUR) THEN
            IF (ISHELA .EQ. IPRNTA .AND.
     *          ISHELB .EQ. IPRNTB .AND.
     *          ISHELC .EQ. IPRNTC .AND.
     *          ISHELD .EQ. IPRNTD) RETURN
         END IF
  400 CONTINUE
  300 CONTINUE
  200 CONTINUE
  100 CONTINUE
C
C UNFOLD DERIVATIVES OF THE FOCK MATRIX
C
      IF(FOCK) THEN
C
       CALL IZERO(JPERT,8)
       CALL IZERO(KPERT,8)
C
C LOOP OVER ALL PERTURBATIONS
C
       DO 852 IPERT=1,NCOORD
C
C DETERMINE THE CORRESPONDING IRREP
C
        DO 853 IREP=1,MAXLOP+1
C
         IF(IPERT.GT.IIPERT(IREP)) THEN
C
          IRREP=IREP
C
C RESET JPERT
C
         ENDIF
C
853     CONTINUE
C
C INCREMENT JPERT
C
        JPERT(IRREP)=JPERT(IRREP)+1
        KPERT(IRREP)=KPERT(IRREP)+1
C
C ELIMINATE TRANSLATIONAL REDUNDANCIES
C
        IF(TRAINV) THEN
C
         IF(IXPERT.EQ.IRREP.AND.KPERT(IRREP).EQ.ITRANSX) THEN
          JPERT(IRREP)=JPERT(IRREP)-1
          GO TO 852
         ENDIF
         IF(IYPERT.EQ.IRREP.AND.KPERT(IRREP).EQ.ITRANSY) THEN
          JPERT(IRREP)=JPERT(IRREP)-1
          GO TO 852
         ENDIF
         IF(IZPERT.EQ.IRREP.AND.KPERT(IRREP).EQ.ITRANSZ) THEN
          JPERT(IRREP)=JPERT(IRREP)-1
          GO TO 852
         ENDIF
        ENDIF
C
C ADDRESS FOR SCRATCH USED IN EXPCAR
C
        ISCR=KSTART+1+NLENT(IRREP)*(1+IUHF)
C
C DIFFER BETWEEN TOTAL SYMMETRIC AND OTHER PERTURBATIONS
C
        IF(IRREP.EQ.1) THEN
C
C TOTAL SYMMETRIC PERTURBATIONS
C
         IF(.NOT.SCF) THEN
C
C DEAL FIRST WITH THE D(REL)*dI/d chi TERMS
C 
          CALL FSYM2A(WORK1(KRFAOA+(IPERT-1)*NF*(MAXLOP+1)), 
     &                WORK1(KRFAOB+(IPERT-1)*NF*(MAXLOP+1)),
     &                NBAST,IPRINT,WORK1(KSTART+1),
     &                WORK1(KSTART+1+NLENT(IRREP)),NIR)
          IF(IDFGH) THEN
           CALL EXPCAR(IRREP,1,WORK1(KSTART+1),WORK1(ISCR),0)
          ENDIF
          CALL PUTLST(WORK1(KSTART+1),
     &                JPERT(IRREP),1,1,IRREP,105)
          IF(IUHF.NE.0) THEN
           IF(IDFGH) THEN
            CALL EXPCAR(IRREP,1,WORK1(KSTART+1+NLENT(IRREP)),
     &                  WORK1(ISCR),0)
           ENDIF
           CALL PUTLST(WORK1(KSTART+1+NLENT(IRREP)),
     &                 JPERT(IRREP),1,1,IRREP,106)
          ENDIF
         ENDIF
C
C UNFOLD FOCK MATRIX DERIVATIVES
C
         CALL FSYM2A(WORK1(KDFAOA+(IPERT-1)*NF*(MAXLOP+1)),
     &               WORK1(KDFAOB+(IPERT-1)*NF*(MAXLOP+1)),
     &               NBAST,IPRINT,WORK1(KSTART+1),
     &               WORK1(KSTART+1+NLENT(IRREP)),NIR)
         IF(IDFGH) THEN
          CALL EXPCAR(IRREP,1,WORK1(KSTART+1),WORK1(ISCR),0)
          IF(IUHF.EQ.1) THEN
           CALL EXPCAR(IRREP,1,WORK1(KSTART+1+NLENT(IRREP)),
     &                 WORK1(ISCR),0)
          ENDIF
         ENDIF
        ELSE
C
C OTHER PERTURBATIONS
C
         IF(.NOT.SCF) THEN 
          CALL FSYM2B(IRREP,WORK1(KRFAOA+(IPERT-1)*NF*(MAXLOP+1)), 
     &                WORK1(KRFAOB+(IPERT-1)*NF*(MAXLOP+1)),
     &                NBAST,IPRINT,WORK1(KSTART+1),
     &                WORK1(KSTART+1+NLENT(IRREP)),NIR)
          ISCR=KSTART+1+NLENT(IRREP)*(1+IUHF)
          IF(IDFGH) THEN
           CALL EXPCAR(IRREP,1,WORK1(KSTART+1+NLENT(IRREP)),
     &                 WORK1(ISCR),0)
          ENDIF
          CALL PUTLST(WORK1(KSTART+1),JPERT(IRREP),
     &                1,1,IRREP,105)
          IF(IUHF.NE.0) THEN
           IF(IDFGH) THEN
            CALL EXPCAR(IRREP,1,WORK1(KSTART+1+NLENT(IRREP)),
     &                  WORK1(ISCR),0)
           ENDIF
           CALL PUTLST(WORK1(KSTART+1+NLENT(IRREP)),
     &                 JPERT(IRREP),1,1,IRREP,106)
          ENDIF
         ENDIF
         CALL FSYM2B(IRREP,WORK1(KDFAOA+(IPERT-1)*NF*(MAXLOP+1)),
     &               WORK1(KDFAOB+(IPERT-1)*NF*(MAXLOP+1)),
     &               NBAST,IPRINT,WORK1(KSTART+1),
     &               WORK1(KSTART+1+NLENT(IRREP)),NIR)
         IF(IDFGH) THEN
          CALL EXPCAR(IRREP,1,WORK1(KSTART+1),WORK1(ISCR),0)
          IF(IUHF.EQ.1) THEN
           CALL EXPCAR(IRREP,1,WORK1(KSTART+1+NLENT(IRREP)),
     &                 WORK1(ISCR),0)
          ENDIF
         ENDIF
        ENDIF
C
C ADD ONE AND TWO-ELECTRON PIECES OF FOCK-MATRIX DERIVATIVES
C
        CALL GETLST(WORK1(KDFAOA+(IPERT-1)*NF*(MAXLOP+1)),
     &              JPERT(IRREP),1,1,IRREP,102)
        CALL SAXPY(NLENT5(IRREP),ONE,WORK1(KDFAOA+(IPERT-1)*NF*
     &             (MAXLOP+1)),1,WORK1(KSTART+1),1)
        CALL PUTLST(WORK1(KSTART+1),
     &              JPERT(IRREP),1,1,IRREP,102)
        IF(IUHF.NE.0) THEN
         CALL SAXPY(NLENT5(IRREP),ONE,WORK1(KDFAOA+(IPERT-1)*NF*
     &              (MAXLOP+1)),1,WORK1(KSTART+1+NLENT(IRREP)),1)
         CALL PUTLST(WORK1(KSTART+NLENT(IRREP)+1),JPERT(IRREP),
     &               1,1,IRREP,103)
        ENDIF
852    CONTINUE
      ENDIF
C
C  Write Final Buffer 
C
      IF (.NOT.EXPECT) THEN
         IF (UNDIFF) THEN
            CALL SYMOUT(WORK1(IWKSO+1),NOABCD,SHAEQB,SHCEQD,SHABAB,
     *                  FIRST,.TRUE.,IPRINT)
         ELSE
            CALL DSOUT2(WORK1(IWKSO+1),NOABCD,IATOM,1,MULTE,MULE,
     *                  SHAEQB,SHCEQD,SHABAB,FIRST,.TRUE.,IPRINT)
         END IF
      END IF
C
C     <<<<< Print Section - Gradient and Hessian Elements >>>>>>
C
      IF (GRAD) THEN
       IF(IPRINT.NE.0) THEN
       CALL HEADER('Two-electron integral gradient',-1)
       CALL PRIGRD(GRADEE,0)
       ENDIF
       IF (MAXDER.EQ.2) THEN
        DO 800 I = 1,NCOORD
         DO 810 J = 1,I
          IF (I .NE. J) THEN
           HESSEE(I,J) = HESSEE(I,J) + HESSEE(J,I)
          END IF
          HESSEE(J,I) = HESSEE(I,J)
810      CONTINUE
800     CONTINUE
        IF(IPRINT.NE.0) THEN
         CALL HEADER('Two-electron integral Hessian',-1)
c         CALL PRIHES(SCF_HESSEE,NCOORD,'CENTERS',0)
         CALL PRIHES(HESSEE,NCOORD,'CENTERS',0)
        ENDIF
       ENDIF
      ENDIF
C
      IF(INTWRIT)THEN
        CALL ACESFLUSH(30,WORK1(KBUCK),WORK1(KIBUCK),WORK1(KCHAIN),
     &             WORK1(KCOUNT),NBUF,NREC,NCOORD)
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
