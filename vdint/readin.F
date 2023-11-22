      SUBROUTINE READIN(WRTONE,WORK1,LWORK1)
C
C     INPUT PROCESSING ROUTINE FOR  V D I N T.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
      INTEGER PQM
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
CSSS#include "ecplims.par"
#include "ecp.par"
#include "mxcent.par"
#include "baslims.par"
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D4 = 4.0D0,
     &           DP25 = 0.25D0, DP5 = 0.5D0, DP75 = 0.75D0,
     &           PI = 3.14159 26535 89793 D0,
     &           BOHRAD = 0.52917 706 D0,
     &           THRMIN = 1.D-15)
C
      LOGICAL SHARE
      LOGICAL BIGVEC,SEGMEN
      LOGICAL BIGVC,GENCON
      LOGICAL ELECT,GEOM,MAGNET,MAGNET2,GIAO
      LOGICAL JSO,JFC,JSD
      LOGICAL LPRINT,WRTONE,LAST
      LOGICAL SHARER(MXQN,MXAOVC)
      LOGICAL DCORD,DCORGD,NOORBT,DOPERT
      CHARACTER NAMEX*6
      CHARACTER*80 MOLBUF
      CHARACTER*1 ISPD(8), IBLANK, KASYM(3,3), ID3, KO(3)
      CHARACTER*2 KOKO(3)
      CHARACTER*4 KWO(MXQN*(MXQN+1)*(MXQN+2)/6), JPRX(MXCORB),
     &            NAMN(MXCENT), LAMN(MXSHEL), MAMN(MXCORB)
      CHARACTER*6 NAVN
C
CCH
      COMMON /FLAGS/IFLAGS(100)
      COMMON /CBIREA/ IPRINT, LUMLCL, MAXPRI, BIGVC, GENCON
      COMMON/PROP/ELECT,MAGNET,GEOM,MAGNET2,GIAO
      COMMON/JNMR/JSO,JFC,JSD
      COMMON /HRUNIT/ LUONEL, LUINTA, LUSUPM, LUDASP, LUSOL
      COMMON /CCOM/ THRS, NHTYP, IBX
C
      DIMENSION DEGEN(MXCENT), CORR(MXCENT,3), IFBA(MXCENT)
      DIMENSION A(MXPRIM), NCCF(MXPRIM), C(MXPRIM,MXCONT), JTRAN(MXCORB)
      DIMENSION KA(8,3), NB(8),
     &          JCO(MXQN), KHK(MXQN), MHK(MXQN), NUCMAX(MXQN),
     &          NUCSUM(MXQN),
     &          NUC(MXQN,MXAOVC), NRC(MXQN,MXAOVC)
      DIMENSION WORK1(LWORK1)
C
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
      COMMON /PRIMIT/ PRIEXP(MXPRIM), PRICCF(MXPRIM,MXCONT),
     *                PRICRX(MXPRIM), PRICRY(MXPRIM), PRICRZ(MXPRIM)
      COMMON /CTRORB/ ICNCTR(MXPRIM)
      COMMON /SHELLSi/ KMAX,
     &                NHKT(MXSHEL),   KHKT(MXSHEL), MHKT(MXSHEL),
     &                ISTBAO(MXSHEL), NUCO(MXSHEL), JSTRT(MXSHEL),
     &                NSTRT(MXSHEL),  MST(MXSHEL),  NCENT(MXSHEL),
     &                NRCO(MXSHEL), NUMCF(MXSHEL),
     &                NBCH(MXSHEL),   KSTRT(MXSHEL)
      COMMON /SHELLS/ CENT(MXSHEL,3), SHARE(MXSHEL)
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
      COMMON /PINCOM/ IPIND(MXCORB), IBLOCK(MXCORB), INDGEN(MXCORB)
      COMMON /AOSOTR/ CTRAN(MXCORB,8), ITRAN(MXCORB,8), IAOAO(MXCORB)
      CHARACTER*1 CHRSGN(-1:1)
      DATA CHRSGN /'-',' ','+'/
      DATA ISPD /'s','p','d','f','g','h','i'/, IBLANK/' '/
      DATA KO   / 'X', 'Y', 'Z'/
      DATA KOKO /'YZ','XZ','XY'/
      DATA KWO /'s   ','x   ','y   ','z   ',
     *          'xx  ','xy  ','xz  ','yy  ','yz  ','zz  ',
     *          'xxx ','xxy ','xxz ','xyy ','xyz ','xzz ',
     *          'yyy ','yyz ','yzz ','zzz ',
     &          'xxxx','xxxy','xxxz','xxyy','xxyz','xxzz',
     &          'xyyy','xyyz','xyzz','xzzz','yyyy','yyyz',
     &          'yyzz','yzzz','zzzz','h500','h410','h401',
     &          'h320','h311','h302','h230','h221','h212',
     &          'h203','h140','h131','h122','h113','h104',
     &          'h050','h041','h032','h023','h014','h005',
     &          'i600','i510','i501','i420','i411','i402',
     &          'i330','i321','i312','i303','i240','i231',
     &          'i222','i213','i204','i150','i141','i132',
     &          'i123','i114','i105','i060','i051','i042',
     &          'i033','i024','i015','i006'/

CCH---------------------------------------------------------------------
C     PARAMETER FOR ECP
C-----------------------------------------------------------------------
CSS      PARAMETER (NDI10=200)
CSS      PARAMETER (NDILMX = 5)
CSS      PARAMETER (MAXJCO = 10)
CCH---------------------------------------------------------------------
C     COMMON FOR ECPGRAD
CCH---------------------------------------------------------------------
      COMMON /ECP   / IECP
CSS      COMMON /PSEUD / NELECP(NDI10),IPSEUD
C-----------------------------------------------------------------------
C     VARIABLES FOR THE ECPSTUFF 
C-----------------------------------------------------------------------
      INTEGER IPR,IUNIT,IERR
CSSS      DIMENSION IPSEUX(NDI10)
CSSS      CHARACTER*80 INFILE
CSSS      CHARACTER*80 CNICK(NDI10)
CSSS      CHARACTER*2 NAMAT(NDI10)
CSSS      CHARACTER MTYPE(NDI10)*(2)
CSSS      DIMENSION IQMSTR(NDI10),JCOSTR(NDI10,NDILMX)
CSSS      DIMENSION NUCSTR(NDI10,NDILMX,MAXJCO),NRCSTR(NDI10,NDILMX,MAXJCO)
C--
      dimension iqmstr(Mxatms),jcostr(Mxatms,Maxang)
      dimension nucstr(Mxatms,Maxang,Mxprim),
     &          nrcstr(Mxatms,Maxang,Mxprim)
      character*4 namat(Mxatms)
      character*80 infile
      character*80 cnick(Mxatms)
      character mtype(Mxatms)*(2)
C-----------------------------------------------------------------------
C     MEANING
C-----------------------------------------------------------------------
C     NDI10   : max. number of atoms
C     NDILMX  : max. l-quantumnumber
C     MAXJCO  : max. number of sets of contractionblocks
C-----------------------------------------------------------------------
C     IPR     : debug level, here: ipr=0
C     IUNIT   : unit number for infile
C     IERR    : errorflag; ierr=0 means no error during the process
C     IPSEUX  : pseudo potential type of atom
C               0  means : no pseudo potential or not yet assigned
C     IPSEUD  : number of types of ecps
C     CNICK   : pseudo potential nicknames indexed by pseudo potential
C               types
C     INFILE  : name of the inputfile for ecp
C     MTYPE   : atomic identifier (starting with element symbol)
C     IQMSTR  : max. l-quantum number of a given shell
C     NAMAT   : name of the atom
C     JCOSTR  : number of different sets of blocks of contraction
C     NUCSTR  : number of primitive exponents of a given shell
C     NRCSTR  : number of blocks of contraction coefficients of a
C               given shell
CCH---------------------------------------------------------------------

      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
C     Statement function to "convert" 0 - 7 to a printable
C     "binary" form
C
      IFAKBN(I) = 25*IBTAND(I,4) + 5*IBTAND(I,2) + IBTAND(I,1)
CCH---------------------------------------------------------------------
C     ECP: GET ECPFLAG (IFLAGS(71)=1: ECP ARE USED!)
C          GET SPHERICALFLAG (IFLAGS(62)=1: SPHERICAL HARMONICS ARE 
C                                           USED!)
CCH---------------------------------------------------------------------
      IECP=IFLAGS(71)
cch      NHARM=IFLAGS(62)
      NHARM=0
CCH---------------------------------------------------------------------
c      IF (IPRINT .GE. 11) CALL TIMER('START ',TIMSTR,TIMEND)
C     ... transfer /CBIREA/ BIGVC to /BLOCKS/ BIGVEC
      BIGVEC = BIGVC
      GENCON = .NOT.SEGMEN
      LUIN = LUMLCL 
C
C  CHANGED TO INTERFACE WITH VMOL, READS NOW MOL FILE
C
      IF (LUIN .NE. LUCMD) THEN
         OPEN (LUIN,FILE='MOL',STATUS='OLD',
     *         FORM='FORMATTED',ERR=3000)
      END IF
C
C     Read two title cards
C
cMOL:1:INTGRL        1    0    1    0    0    0    0    0    0
cMOL:2:         *** ACES2 Program System (Release V0.1) ***
cMOL:3:Cyclopentane, 40 electrons
      READ(LUIN,'(A)') MOLBUF
      IF (MOLBUF(1:6).EQ.'INTGRL'.OR.MOLBUF(1:6).EQ.'SUPMAT') THEN
         READ(LUIN,'(A)') MOLBUF
      END IF
      READ(LUIN,'(A)') MOLBUF
C
C NONTYP NUMBER OF NON-SYMMETRY-EQUIVALENT NUCLEI
C NSYMOP NUMBER oF NON-REDUNDANT SYMMETRY OPERATIONS
C
cMOL:4:   15    0            0.10E-13              0    0
c------12345678901234567890123456789012345678901234567890
      READ(LUIN,'(A)') MOLBUF
      READ(MOLBUF,1000) ITE,NONTYP,NSYMOP,
     *       ((KASYM(I,J),I=1,3),J=1,3),ID3,THRS
 1000 FORMAT (BN,I1,I4,I5,10A1,D10.2,6I5)
      NPPYP = NONTYP
      IF (IPRINT.EQ.0) THEN
         LPRINT = .FALSE.
      ELSE
         LPRINT = .TRUE.
      END IF
CJDW 3/31/97.
      IF(IFLAGS(1).GE.10)THEN
         LPRINT = .TRUE.
      ELSE
         LPRINT = .FALSE.
      ENDIF
1010  FORMAT(//T2,78('*'),/T2,'*',76X,'*',
     *        /T2,'*',2X,A72,2X,'*',/T2,'*',76X,'*',
     *        /T2,'*',2X,A72,2X,'*',/T2,'*',76X,'*',
     *        /T2,78('*'),//)
C
      IF (THRS .LE. D0) THEN
         WRITE (6,'(/A,2(/T6,A,1P,D10.2),/)')
     &      ' *** WARNING FROM READIN ***',
     &      'Threshold for discarding integrals was',THRS,
     &      'Threshold is reset to                 ',THRMIN
         THRS = THRMIN
      ENDIF
C
C  SET THRESHOLD TO 10.D-14 AS DEFAULT
C
C THRS THRESHOLD FOR NEGLECTING INTEGRALS
C
      THRS=1.D-14
C
C     Print threshold for integrals
C
      IF (LPRINT) WRITE(6,1020) THRS
 1020 FORMAT (/T11,'Integrals less than ',1P,D9.2,' are neglected')
C
C PROCESS SYMMETRY INPUT 
C
      IF (LPRINT) WRITE (LUPRI,'(//I11,A,/I11,A/)')
     &     NONTYP, ' types of atoms',
     &     NSYMOP, ' symmetry operations'
C
C     DETERMINE ISYTYP(I)
C
      ISYTYP(1) = 0
      ISYTYP(2) = 0
      ISYTYP(3) = 0
C
C MAXLOP + 1 = NUMBER OF IRREDUCIBLE REPRESENTATIONS,
C THE LOOPS ARE ALWAYS FROM 0 TO MAXLOP
C
C  ISYTYP DETERMINES THE BEHAVIOUR OF X,Y,Z UNDER THE SYMMETRY
C  OPERATIONS : BIT 1,2,3 ARE SET TO ONE IF X, OR Y OR Z ARE
C  CHANGING UNDER THIS PARTICULAR OPERATIONS. HOWEVER, ISYTYP
C  ALSO GIVES THE IRREP OF X,Y,Z WHEN COUNTING FROM 0 TO MAXLOP
C
      MAXLOP = 2**NSYMOP - 1
      IF (NSYMOP.GT.0) THEN
       DO 111 J = 1,NSYMOP
        NBL = 1
        DO 112 I = 1,3
         IF (KASYM(I,J).EQ.IBLANK) NBL = NBL + 1
         IF (KASYM(I,J).EQ.KO(1)) ISYTYP(1) = ISYTYP(1) + 2**(J-1)
         IF (KASYM(I,J).EQ.KO(2)) ISYTYP(2) = ISYTYP(2) + 2**(J-1)
         IF (KASYM(I,J).EQ.KO(3)) ISYTYP(3) = ISYTYP(3) + 2**(J-1)
  112   CONTINUE
        GO TO (21,22,23),NBL
   21   IF (LPRINT) WRITE (LUPRI,'(9X,A)') ' Inversion centre'
        GO TO 111
   22   DO 42 K = 1,3
         DO 41 I = 1,3
          IF (KASYM(I,J) .EQ. KO(K)) GO TO 42
   41    CONTINUE
         IF (LPRINT) WRITE (LUPRI,'(9X,3A)')
     &               ' Rotation about the ',KO(K),'-axis'
   42   CONTINUE
        GO TO 111
   23   DO 52 K = 1,3
         DO 51 I = 1,3
          IF (KASYM(I,J) .EQ. KO(K)) GO TO 53
   51    CONTINUE
   52   CONTINUE
        GO TO 111
   53   IF (LPRINT) WRITE (6,'(9X,3A)')
     *       ' Reflection in the ', KOKO(K), '-plane'
  111  CONTINUE
      END IF
C
C     Determine KA(I,J), PT(I), MULT(I) and FMULT(I)
C
C     Note new ranges for PT, MULT and FMULT
C
      DO 101 I = 1,8
       IS = 0
       DO 107 K = 1,3
        KA(I,K) = (I-1)/(2**(K-1))-2*((I-1)/(2**K))
        IS = IS + KA(I,K)
  107  CONTINUE
       PT(I-1)    = (-1.0)**IS
       MULT(I-1)  = 2**MAX(0,NSYMOP-IS)
       FMULT(I-1) = FLOAT(MULT(I-1))
  101 CONTINUE
C
C     Determine ITYPE(I,J), note that these are renumbered
C     according to new FMULT etc definitions
C
      ITYPE(1,1) = 0
      NDO = 0
      NFO = 0
      NGO = 0
      NHO = 0
      NIO = 0
     `
      DO 115 JQM = 1,3
       ITYPE(2,JQM) = ISYTYP(JQM)
       DO 116 KQM = JQM,3
        NDO   = NDO+1
        IVARA = ISYTYP(JQM)
        IVARB = ISYTYP(KQM)
        IVARC = IBTXOR(IVARA,IVARB)
        ITYPE(3,NDO) = IVARC
        DO 117 LQM = KQM,3
         NFO   = NFO+1
         IVARD = ISYTYP(LQM)
         IVARE=IBTXOR(IVARC,IVARD)
         ITYPE(4,NFO) = IVARE
         DO 118 MQM = LQM, 3
          NGO  = NGO+1
          IVARF=ISYTYP(MQM)
          IVARG=IBTXOR(IVARE,IVARF)
          ITYPE(5,NGO)=IVARG
          DO 119 NQM = MQM, 3
           NHO = NHO+1
           IVARH=ISYTYP(NQM)
           IVARI=IBTXOR(IVARG,IVARH)
           ITYPE(6,NHO)=IVARI
           DO 120 PQM = NQM, 3
            NIO = NIO+1
            IVARP=ISYTYP(PQM)
            IVARQ=IBTXOR(IVARI,IVARP)
            ITYPE(7,NIO)=IVARQ
  120       CONTINUE
  129      CONTINUE
  119     CONTINUE
  118    CONTINUE
  117   CONTINUE
  116  CONTINUE
  115 CONTINUE
C
C  Process orbital input data 
C
C     Run over types of atoms in the molecule
C
      PIPPI = (DP5/PI)**DP75
      KMAX  = 0
      KBCH  = 0
      LMAX  = 0
      NMAX  = 0
      MAXLO = 7
      IND   = 0
      MMAX  = 0
      NHTYP = 2
      JATOM = 0
cMOL:5:   9999.00      3.00
      READ(LUIN,'(A)') MOLBUF
      DO 5 I = 1, NONTYP
cMOL:6:          6.00000000    1    3    1    1    1
c------123456789012345678901234567890123456789012345
         READ(LUIN,'(A)') MOLBUF
         READ(MOLBUF,1030) MAXNUC,Q,NONT,IQM,(JCO(J),J=1,IQM)
         Do J = 1, IQM
           MHK(J) = 0
         ENDDO
 1030    FORMAT (BN,I5,1X,F14.8,12I5)
CCH---------------------------------------------------------------------
C        ECP: GET IQM AND JCO
CCH---------------------------------------------------------------------
         IF (IECP.NE.0) THEN
           IQMSTR(I)=IQM
           DO IZAHL=1,IQM
             JCOSTR(I,IZAHL)=JCO(IZAHL)
           ENDDO
         ENDIF
CCH---------------------------------------------------------------------
         IF (MAXNUC .EQ. 0) MAXNUC = MAXPRI
         IF (MAXNUC .LT. 0) Q = ABS(FLOAT(MAXNUC))*Q
         IF (LPRINT) THEN
            IF (IQM .GT. 0) THEN
               WRITE (6,1040)  I,NINT(Q),NONT,ISPD(IQM)
            ELSE
               WRITE (6,1050)  I,NINT(Q),NONT
            END IF
            IF (IQM.GT.0) WRITE (6,1060) (JCO(J),ISPD(J),J=1,IQM)
         END IF
 1040    FORMAT(////10X,'Atomic type number',I2,
     *             /10X,'--------------------',
     *            //10X,'Nuclear charge:                      ',I3,
     *             /10X,'Number of symmetry independent atoms:',I3,
     *             /10X,'Highest orbital type:                ',2X,A1,
     *             /)
 1050    FORMAT (////10X,'Atomic type number',I2,
     *             /10X,'--------------------',
     *            //10X,'Nuclear charge:                      ',I3,
     *             /10X,'Number of symmetry independent atoms:',I3,
     *             /10X,'No orbitals attached',/)
 1060    FORMAT(11X,I3,' CGTO''s of ',A1,' type')
C
C        ***** KHK() and MHK() *****
C
         DO 1342 J = 1, IQM
            KHKTA = J*(J+1)/2
            IF (MHK(J) .EQ. 0) THEN
               MHK(J) = 268435455
               KHK(J) = KHKTA
            ELSE
               IF (NSYMOP .GT. 0) THEN
                  WRITE (6, '(/A/)')
     *              ' No symmetry allowed when components are deleted.'
                  STOP
               END IF
               MHKTA = MHK(J)
               ICMP = 0
               DO 1323 K = 1, KHKTA
                  IF (IBTAND(IBTSHR(MHKTA,K-1),1) .EQ. 1) THEN
                     ICMP = ICMP + 1
                  END IF
 1323          CONTINUE
               KHK(J) = ICMP
               WRITE (6,'(/10X,3A,15I1)') ' The following ',
     *              ISPD(J),' components are used: ',
     *              (IBTAND(IBTSHR(MHKTA,K-1),1),K=1,KHKTA)
            END IF
 1342    CONTINUE
         NMIN = NMAX+1
C
C        Run over symmetry independent atoms of this type
C
         DO 113 N = 1,NONT
C
C           Read coordinates
C
cMOL:7:C #1     -1.445961800726     -1.946295242451     -0.236376233553
            READ(LUIN,'(A)') MOLBUF
            IF (MOLBUF(72:72).EQ.'*') THEN
               READ(MOLBUF,1212) NAVN(1:4),(CORR(N,J), J = 1,3)
            ELSE
               READ(MOLBUF,12) NAVN(1:4),(CORR(N,J), J = 1,3)
            END IF
   12       FORMAT(BN,A4,3F20.12)
 1212       FORMAT(BN,A4,3F20.15)
            JATOM = JATOM + 1
            DO 435 II = 1, 3
               IF (NAVN(1:1) .EQ. ' ') NAVN(1:4) = NAVN(2:4)//' '
  435       CONTINUE
            NAMEX(3*JATOM)     = NAVN(1:4)//' z'
            NAMEX(3*JATOM - 1) = NAVN(1:4)//' y'
            NAMEX(3*JATOM - 2) = NAVN(1:4)//' x'
            IF(I.LE.NPPYP) MMAX=MMAX+1
            NMAX = NMAX + 1
            IFBA(NMAX) = IQM
            READ (NAVN,'(A4)') NAMN(NMAX)
CCH---------------------------------------------------------------------
C           ECP: GET ATOMSYMBOL
CCH---------------------------------------------------------------------
            NAMAT(NMAX)=NAMN(NMAX)(1:2)
CCH---------------------------------------------------------------------
            DO 114 J = 1,3
               IF (ID3.NE.IBLANK) CORR(N,J) = CORR(N,J)/BOHRAD
               CORD(NMAX,J) = CORR(N,J)
  114       CONTINUE
  113    CONTINUE
         IF (IPRINT.GE.2) THEN
            WRITE(6,1293)
            WRITE(6,1211) (NAMN(N),
     *         (CORD(N,J),J = 1,3),N = NMIN,NMAX)
         END IF
 1293    FORMAT (//T11,'Nuclear coordinates:',/)
 1211    FORMAT(10X,A4,3X,3F12.6)
         IMAX=0
         IF (IQM.GT.0) THEN
C
C AP: IQM is the maximum angular momentum of a symmetry unique atom of a given
C type (s=1,p=2,d=3.....)
C
C           Run over angular quantum number
C
            if (iqm.gt.mxqn) then
               print *, '@READIN: Assertion failed.'
               print *, '         mxqn = ',mxqn
               print *, '         iqm  = ',iqm
               call errex
            end if
            DO 9002 J=1,IQM
               NUCMAX(J) = 0
               NUCSUM(J) = 0
               KM = JCO(J)
C
C AP: KM is the number of shells (here "shell" is loosely used and entirely
C depends on how the user chooses to input the basis set)
C
               IF (KM .GT. 0) THEN
C
C                 Run over sets of primitives to be contracted
C
                  if (km.gt.mxaovc) then
                     print *, '@READIN: Assertion failed.'
                     print *, '         mxaovc = ',mxaovc
                     print *, '         km     = ',km
                     call errex
                  end if
                  DO 9004 K = 1, KM
C
C                    Number of primitive and contracted functions
C
                     READ(LUIN,'(A)') MOLBUF
                     READ(MOLBUF,'(BN,10I5)') NUC(J,K), NRC(J,K)
                     SHARER(J,K) = NUC(J,K) .LT. 0
                     NUC(J,K)    = ABS(NUC(J,K))
C
C AP: The current maxima allowed for J and K are 5 (MXQN) and 32(MXAOVC).
C 5 is the maximum angular momentum and 32 is the number of "shells". The 32
C seems too high at this point, but we need to go further and see how MXAOVC
C is being used. The maxima of NUCMAX, NUSUM, NUC, and NRC are set correctly.
C
CCH---------------------------------------------------------------------
C     ECP: GET NUC AND NRC
CCH---------------------------------------------------------------------
c     I : counter of atoms
c     J : counter of l-quantumnumbers of the atom
c     K : counter of jco
C
C    j ---> ....
C
                     IF (IECP.NE.0) THEN
                       IF (K.LE.JCOSTR(I,J)) THEN
                         NUCSTR(I,J,K)=NUC(J,K)
                         NRCSTR(I,J,K)=NRC(J,K)
                       ENDIF
                     ENDIF
C
C AP: CH has defined NUCSTR and NRCSTR to be 50, 5, and 10. The 10 in principle
C has to be 32 to be compatible. It seems to be a counter for symmetry unique
C atoms and that has to tightened to maxatms.
CCH---------------------------------------------------------------------
                     IMIN = IMAX + 1
                     IMAX = IMAX + NUC(J,K)
C
C AP: IMIN is a counter for the total primitives prior to the processing of the
C current "shell". IMAX is a counter for the total number of primitives.
C
                     NUCMAX(J) = MAX(NUCMAX(J),NUC(J,K))
                     NUCSUM(J) = NUCSUM(J) + NUC(J,K)
C
C AP: NUCMAX is the maximum number of primitives allowed for a given angular
C momentum. NUCSUM is the total number of primitives for a given angular mom.
C
C                    Exponents and contraction matrix
C
                     if (imax.gt.mxprim.or.nrc(j,k).gt.mxcont) then
                        print *, '@READIN: Assertion failed.'
                        print *, '         imax = ',imax,' > ',mxprim
                        print *, '         nrc = ',nrc(j,k),' > ',mxcont
                        call errex
                     end if
                     DO 9222 L = IMIN,IMAX
cMOL:9:    6665.000000000       0.000692000      -0.000146000       0.000000000
c------123456789012345678901234567890123456789012345678901234567890123456789012
CSSS                        READ (LUIN,'(4F18.4)')
CSSS     *                        A(L), (C(L,M), M=1,NRC(J,K))

                        READ (LUIN,*)
     *                        A(L), (C(L,M), M=1,NRC(J,K))
                        NCCF(L) = NRC(J,K)
 9222                CONTINUE
C
C AP: The array A contains the exponents and C contains the contraction coefs.
C Important to note that they contain the exponents and coefficients
C for all the symmetry unique atoms. The current maximum for A is 5184 (MXPRIM)
C and for C it is 5184 (MXPRIM) by 32 (MXCONT=MXAOVC). So, currently we can
C process 165888 primitives. NCCF (MAXPRIM) has the number of contracted
C functions for a set of primitives.
C
C                    Reorder primitive orbitals
C
                     IF (IMIN.LT.IMAX) THEN
                        IMAM = IMAX-1
                        DO 9069 L=IMIN,IMAM
                           LP = L + 1
                           DO 9068 M=LP,IMAX
                           IF(A(L) .GT. A(M)) GO TO 9068
                              HOLD = A(L)
                              A(L) = A(M)
                              A(M) = HOLD
                              DO 9071 N = 1, NRC(J,K)
                                 HOLD   = C(L,N)
                                 C(L,N) = C(M,N)
                                 C(M,N) = HOLD
 9071                         CONTINUE
 9068                      CONTINUE
 9069                   CONTINUE
                     END IF
C
C                    Normalize orbitals
C
                     DO 9070 N = 1, NRC(J,K)
                        SUM = D0
                        DO 9001 L=IMIN,IMAX
                        DO 9001 M=IMIN,L
                           T = D2*SQRT(A(L)*A(M)) / (A(L)+A(M))
                           T = C(L,N)*C(M,N)*(T**(J + DP5))
                           IF (L .NE. M) SUM = SUM + T
                           SUM = SUM + T
 9001                   CONTINUE
                        SUM=D1/SQRT(SUM)
                        DO 9016 L=IMIN,IMAX
                           C(L,N)  = C(L,N)
     *                         *SUM*(D4*A(L))**(DP5*J+DP25)*PIPPI
 9016                   CONTINUE
 9070                CONTINUE
 9004             CONTINUE
C
C                 Determine NUCMAX
C
                  IF (MAXNUC .EQ. 0) THEN
                     NUCMAX(J) = NUCSUM(J)
                  ELSE
                     NUCMAX(J) = MAX(NUCMAX(J),MAXNUC)
                     NUCMAX(J) = MIN(NUCMAX(J),NUCSUM(J))
                  END IF
               END IF
 9002       CONTINUE
         END IF
C
C AP: NUCMAX is the maximum number of primitives for an angular momentum.
C
C        End of orbital processing for this atom
C
C        Run over symmetry independent atoms
C
         DO 6 N = NMIN,NMAX
C
C AP: Actually NMIN (as far as I can tell) and NMAX are the total number of
C symmetry unique atoms of given type.
C
C           PRIEXP(), PRICCF(), PRICRX(), PRICRY(), PRICRZ()
C
            IF (IMAX.GT.0) THEN
               if (lmax+imax.gt.mxprim) then
                  print *, '@READIN: Assertion failed.'
                  print *, '         mxprim = ',mxprim
                  print *, '         lmax   = ',lmax+imax
                  call errex
               end if
               DO 7 L = 1,IMAX
                  LMAX = LMAX + 1
                  PRIEXP(LMAX) = A(L)
C
C AP: NCCF (lim to MXPRIM) has the number of contracted functions for a set of
C primitives. MXCONT (=MXAOVC) is the current maximum for an element of
C this array. The maximum for PRIEXP (MXPRIM) is correctly set and
C the assertions are correct.
C
                  if (nccf(l).gt.mxcont) then
                     print *, '@READIN: Assertion failed.'
                     print *, '         mxcont = ',mxcont
                     print *, '         nccf   = ',nccf(l)
                     call errex
                  end if
                  DO 9007  J = 1, NCCF(L)
                     PRICCF(LMAX,J) = C (L,J)
 9007             CONTINUE
                  PRICRX(LMAX) = CORD(N,1)
                  PRICRY(LMAX) = CORD(N,2)
                  PRICRZ(LMAX) = CORD(N,3)
    7          CONTINUE
            END IF
C
C AP: The maxima for PRICCF (MXPRIM and MXCONT) are correct. Also,
C the maximum for PRICR{X,Y,Z} is correctly set to MXPRIM.
C
C           MULK and MAXLO, new definitions
C
            MULK = 1
            LL   = 1
            DO 9 L = 1, NSYMOP
               DO 8 M = 1, 3
                  IVARA = ISYTYP(M)
                  IF (IBTAND(LL,IVARA) .NE. 0) THEN
                     IF (ABS(CORD(N,M)) . GE . 1.D-06) GO TO 9
                  END IF
    8          CONTINUE
               MULK = MULK + LL
    9       LL = 2*LL
C
            MULK = MULK - 1
            MAXLO     = IBTAND(MULK,MAXLO)
            ISTBNU(N) = MULK
            CHARGE(N) = Q
            NOORBT(N) = IQM .EQ. 0
            NHTYP     = MAX(NHTYP,IQM)
            IND       = IND + 1
            IF ( IQM.GT.0 ) THEN
               DO 1274  J=1,IQM
                  DO 1275 K = 1, JCO(J)
                     KBCH = KBCH + 1
c                     DO 1277 KK = 1, NRC(J,K)
                        KMAX = KMAX + 1
                        NCENT (KMAX) = IND
                        NUCO  (KMAX) = NUC(J,K)
                        NRCO  (KMAX) = NRC(J,K)
c                        NUMCF (KMAX) = KK
                        NUMCF (KMAX) =  1
                        NBCH  (KMAX) = KBCH
                        SHARE (KMAX) = SHARER(J,K)
                        LAMN  (KMAX) = NAMN(N)
                        ISTBAO(KMAX) = MULK
                        NHKT  (KMAX) = J
                        KHKT  (KMAX) = KHK(J)
                        MHKT  (KMAX) = MHK(J)
                        DO 1276 M = 1, 3
                           CENT(KMAX,M) = CORD(N,M)
 1276                   CONTINUE
 1277                CONTINUE
 1275             CONTINUE
 1274          CONTINUE
            END IF
    6    CONTINUE
    5 CONTINUE
     
C
C AP: The variable KMAX is a counter for the total number of "shells". Note
C that DO 1277 is commented (in the current version). All of the arrays indexed
C to KMAX have the maximum dimension based on MXSHEL=MXTNPR.
C
C     **********************************************
C     ***** Tabulate Incomplete Gamma Function *****
C     **********************************************
C
      JMAX = MAX(4,4*(NHTYP - 1) + 2)
      CALL GAMTAB(JMAX)
C
C     ***** MAXLOT *****
C
      MAXLOT = 1
      IF (NSYMOP .GT. 0) THEN
         II = 1
         DO 17 I = 1, NSYMOP
            IF(IBTAND(II,MAXLO).EQ.0) MAXLOT = 2*II
            II = 2*II
   17    CONTINUE
      END IF
      MAXLOT = MAXLOT - 1
C
C     ***** NSTRT() AND KSTRT() *****
C
      NSTRT(1) = 0
      KSTRT(1) = 0
      DO 9030 IA = 1, KMAX - 1
         NSTRT(IA + 1) = NSTRT(IA) + KHKT(IA)*NRCO(IA)
         KSTRT(IA + 1) = KSTRT(IA) + KHKT(IA)*MULT(ISTBAO(IA))*NRCO(IA)
 9030 CONTINUE
C
C     ***** JSTRT() *****
C
      NOFF = 0
      DO 9031 IA = 1, KMAX - 1
         JSTRT(IA) = NOFF
         IF (NUMCF(IA + 1). EQ. 1) NOFF = NOFF + NUCO(IA)
 9031 CONTINUE
      JSTRT(KMAX) = NOFF
C
C     ******************************************
C     ***** Orbital Contraction Processing *****
C     ******************************************
C
C AP: KMAX corresponds to the number of "shells" (maxed at MXSHEL). In the IA
C loop, what we do is to take into account the symmetry non-unique atoms. The
C degeneracy of p(=3), d(=6), ... is also taken care of. So, all the arrays
C that are indexed with IABAS should have maximum dims set to the maximum
C number of basis functions MXCORB.
C
      IABAS = 0
      IBBAS = 0
      DO 30 IA = 1, KMAX
         MULA  = ISTBAO(IA)
         NUCA  = NUCO(IA)
         NRCA  = NRCO(IA)
         NUMA  = NUMCF(IA)
         JSTA  = JSTRT(IA)
         NHKTA = NHKT(IA)
         KHKTA = KHKT(IA)
         MHKTA = MHKT(IA)
         IF (IA .LT. KMAX) THEN
            LAST  = NUMCF(IA + 1) .EQ. 1
         ELSE
            LAST = .TRUE.
         END IF
         IF (NUCA.GT.MXAOVC) THEN
            print *, '@READIN: Assertion failed.'
            print *, '         # AO/sh =',NUCA
            print *, '         allowed =',MXAOVC
            call errex
         END IF
         IF (NRCA.GT.MXCONT) THEN
            print *, '@READIN: Assertion failed.'
            print *, '         # AO/cn =',NRCA
            print *, '         allowed =',MXCONT
            call errex
         END IF
         DO 300 NA = 1, NHKTA*(NHKTA+1)/2
            IF (IBTAND(IBTSHR(MHKTA,NA-1),1) .EQ. 1) THEN
               MXYZ = ((NHKTA+1)*NHKTA*(NHKTA-1))/6 + NA
               DO 1516 KB = 0, MAXLOP
                  IF (IBTAND(KB,MULA) .EQ. 0) THEN
                     IF (IABAS+nrca.GT.MXCORB) THEN
                        print *, '@READIN: Assertion failed.'
                        print *, '         # cont. AOs =',IABAS+nrca
                        print *, '         allowed max =',MXCORB
                        call errex
                     END IF
                     IF (IBBAS+NUCA.GT.MXPRIM) THEN
                        print *, '@READIN: Assertion failed.'
                        print *, '         # prim. AOs =',IBBAS+NUCA
                        print *, '         allowed max =',MXPRIM
                        call errex
                     END IF
                     do 1517 irc=1,nrca
                     IABAS = IABAS + 1
                     ICNCTR(IABAS) = NCENT(IA)
1517                 continue
                     IF (LAST) IBBAS = IBBAS + NUCA
                  END IF
 1516          CONTINUE
            END IF
  300    CONTINUE
   30 CONTINUE
C
C     ***** MST(I) *****
C
      MST(1) = 0
      DO 3 IA = 1, KMAX - 1
         MAD = 0
         DO 15 IB = 1, IA
            MAD = MAD + NUCO(IA)*NUCO(IB)
   15    CONTINUE
         MST(IA + 1) = MST(IA) + MAD*(MAXLOT+1)
    3 CONTINUE
C
C   Orbital Symmetry Processing 
C
C     Determine the following variables:
C
C     CTRAN,IPIND,IPTSYM,ITRAN,JPRX,JTRAN,MAMN,
C     NAOS,NB,NPAR,NPARNU,NPARSU,NORB,NSYM
C
      CALL SYMPRO(CTRAN,ITRAN,JPRX,JTRAN,KWO,MAMN,NB,NORB,
     *            NSYM,LAMN,IPRINT)
C
C     IPTAX(3) and NAXREP(0:7)
C
      IF(.NOT.MAGNET) THEN
       IPTAXI = 0
       DO 1001 IREP = 0, MAXLOP
        NAXIS = 0
        DO 1002 ICOOR = 1, 3
         IF (IBTXOR(IREP,ISYTYP(ICOOR)) .EQ. 0) THEN
          NAXIS  = NAXIS + 1
          IPTAXI = IPTAXI + 1
          IPTAX(ICOOR) = IPTAXI
          IPTXYZ(NAXIS,IREP) = ICOOR
         END IF
1002    CONTINUE
        NAXREP(IREP) = NAXIS
1001   CONTINUE
      ELSE 
       IPTAXI=0
       DO 1003 IREP=0,MAXLOP
        NAXIS=0
        IF(IBTXOR(IREP,IBTXOR(ISYTYP(2),ISYTYP(3))).EQ.0) THEN
         NAXIS=NAXIS+1
         IPTAXI=IPTAXI+1
         IPTAX(1)=IPTAXI
         IPTXYZ(NAXIS,IREP)=1
        ENDIF
        IF(IBTXOR(IREP,IBTXOR(ISYTYP(1),ISYTYP(3))).EQ.0) THEN
         NAXIS=NAXIS+1
         IPTAXI=IPTAXI+1
         IPTAX(2)=IPTAXI
         IPTXYZ(NAXIS,IREP)=2
        ENDIF
        IF(IBTXOR(IREP,IBTXOR(ISYTYP(2),ISYTYP(1))).EQ.0) THEN
         NAXIS=NAXIS+1
         IPTAXI=IPTAXI+1
         IPTAX(3)=IPTAXI
         IPTXYZ(NAXIS,IREP)=3
        ENDIF
        NAXREP(IREP)=NAXIS
1003   CONTINUE
      ENDIF
      IF (IPRINT.GT. 20) THEN
       WRITE (6,'(1X,A,3I5)') ' IPTAX  ', (IPTAX(I),I=1,3)
       WRITE (6,'(1X,A,8I5)') ' NAXREP ', (NAXREP(I),I=0,MAXLOP)
      END IF
C
C Calculate symmetry-adapted nuclear coordinates, in the case
C of magnetic pertubations, calculate symmetry adapted nuclear
C magnetic moments. Note while (x,y,z) coordinates transform
C as a vector, the nuclear magnetic moments transform as
C pseudo-vectors in an analogous way as the magnetic field.
C
      DO 7030 IREP = 0,MAXLOP
         DO 7040 ICENTA = 1, 3*NMAX
            IPTCNT(ICENTA,IREP) = 0
7040     CONTINUE
7030  CONTINUE
      IOFF = 0
      IF(.NOT.MAGNET) THEN
C
C NUCLEAR COORDINATES
C
       IF(.NOT.JSD.AND..NOT.JFC) THEN
C
        DO 7000 IREP = 0,MAXLOP
         ICENTA = 0
         DO 7010 ICENT = 1,NMAX
          MULC = ISTBNU(ICENT)
          DO 7020 IDIRN = 1,3
           ICENTA = ICENTA + 1
           IF (IBTAND(MULC,IBTXOR(IREP,ISYTYP(IDIRN))) .EQ. 0) THEN
            IOFF = IOFF + 1
            IPTCNT(ICENTA,IREP) = IOFF
           ENDIF
7020      CONTINUE
7010     CONTINUE
7000    CONTINUE
C
       ELSE IF(JFC) THEN
C
        DO 7002 IREP=0,MAXLOP
         ICENTA=0
         DO 7012 ICENT=1,NMAX
          MULC=ISTBNU(ICENT)
          ICENTA=ICENTA+1
          IF(IBTAND(MULC,IREP).EQ.0) THEN
           IOFF=IOFF+1
           IPTCNT(ICENTA,IREP)=IOFF
          ENDIF
7012     CONTINUE
7002    CONTINUE
       ELSE IF(JSD) THEN
C
        DO 7003 IREP = 0,MAXLOP
         ICENTA = 0
         DO 7013 ICENT = 1,NMAX
          MULC = ISTBNU(ICENT)
          DO 7023 IDIRN = 1,3
           DO 7024 JDIRN=IDIRN,3
            ICENTA = ICENTA + 1
            ISYTIJ=IBTXOR(ISYTYP(IDIRN),ISYTYP(JDIRN))
            IF (IBTAND(MULC,IBTXOR(IREP,ISYTIJ)) .EQ. 0) THEN
            IOFF = IOFF + 1
            IPTCNT(ICENTA,IREP) = IOFF
           ENDIF
7024       CONTINUE
7023      CONTINUE
7013     CONTINUE
7003    CONTINUE
       ENDIF
C
      ELSE
C
C NUCLEAR MAGNETIC MOMENTS
C
       DO 7001 IREP = 0,MAXLOP
        ICENTA = 0
        DO 7011 ICENT = 1,NMAX
         MULC = ISTBNU(ICENT)
         DO 7021 IDIRN = 1,3
          IF(IDIRN.EQ.1) THEN
           ISYMTYP=IBTXOR(ISYTYP(2),ISYTYP(3))
          ELSE IF(IDIRN.EQ.2) THEN
           ISYMTYP=IBTXOR(ISYTYP(1),ISYTYP(3))
          ELSE IF(IDIRN.EQ.3) THEN
           ISYMTYP=IBTXOR(ISYTYP(1),ISYTYP(2))
          ENDIF
          ICENTA = ICENTA + 1
          IF (IBTAND(MULC,IBTXOR(IREP,ISYMTYP)) .EQ. 0) THEN
           IOFF = IOFF + 1
           IPTCNT(ICENTA,IREP) = IOFF
          ENDIF
7021     CONTINUE
7011    CONTINUE
7001   CONTINUE
      ENDIF
      IF(.NOT.MAGNET) THEN
       IF(JFC) THEN
        CALL HEADER('Symmetry-adapted nuclei',-1)
       ELSE IF(JSD) THEN
        CALL HEADER('Symmetry-adapted spin-dipolar integrals',-1) 
       ELSE
CJDW 3/31/97
        IF(LPRINT)THEN
         CALL HEADER('Symmetry-adapted nuclear coordinates',-1)
        ENDIF
       ENDIF
      ELSE
       CALL HEADER('Symmetry-adapted nuclear magnetic moments',-1)
      ENDIF
      ICENTA = 0
      IF(JFC) THEN
       WRITE(6,'(11X,8I4)') (IREP, IREP = 0,MAXLOP)
       DO 7048 ICENT=1,NMAX
        WRITE(6,'(4X,A4,3X,8I4)')
     &  NAMN(ICENT),(IPTCNT(ICENTA+1,IREP),IREP = 0,MAXLOP)
        ICENTA=ICENTA+1
7048   CONTINUE
      ELSE IF(JSD) THEN
       WRITE(6,'(12X,8I4)') (IREP, IREP = 0,MAXLOP)
       DO 7049 ICENT=1,NMAX
        WRITE(6,'(4X,A4,2X,A2,8I4)')
     &  NAMN(ICENT),KWO(5),(IPTCNT(ICENTA+1,IREP),IREP=0,MAXLOP)
        WRITE(6,'(4X,A4,2X,A2,8I4)')
     &  NAMN(ICENT),KWO(6),(IPTCNT(ICENTA+2,IREP),IREP=0,MAXLOP)
        WRITE(6,'(4X,A4,2X,A2,8I4)')
     &  NAMN(ICENT),KWO(7),(IPTCNT(ICENTA+3,IREP),IREP=0,MAXLOP)
        WRITE(6,'(4X,A4,2X,A2,8I4)')
     &  NAMN(ICENT),KWO(8),(IPTCNT(ICENTA+4,IREP),IREP=0,MAXLOP)
        WRITE(6,'(4X,A4,2X,A2,8I4)')
     &  NAMN(ICENT),KWO(9),(IPTCNT(ICENTA+5,IREP),IREP=0,MAXLOP)
        WRITE(6,'(4X,A4,2X,A2,8I4)')
     &  NAMN(ICENT),KWO(10),(IPTCNT(ICENTA+6,IREP),IREP=0,MAXLOP)
        ICENTA=ICENTA+6
7049   CONTINUE

      ELSE
CJDW 3/31/97
       IF(LPRINT)THEN
        WRITE(6,'(11X,8I4)') (IREP, IREP = 0,MAXLOP)
        DO 7050 ICENT = 1,NMAX
         WRITE(6,'(4X,A4,2X,A1,8I4)')
     &   NAMN(ICENT), KWO(2), (IPTCNT(ICENTA+1,IREP),IREP = 0,MAXLOP)
         WRITE(6,'(4X,A4,2X,A1,8I4)')
     &   NAMN(ICENT), KWO(3), (IPTCNT(ICENTA+2,IREP),IREP = 0,MAXLOP)
         WRITE(6,'(4X,A4,2X,A1,8I4)')
     &   NAMN(ICENT), KWO(4), (IPTCNT(ICENTA+3,IREP),IREP = 0,MAXLOP)
         ICENTA = ICENTA + 3
7050    CONTINUE
       ENDIF
      ENDIF
C
C     ***** Calculate Distances and Nuclear Potential Energy *****
C
      POTNUC = D0
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('Internuclear distances in angstroms',-1)
      END IF
      DO 8724 N = 1, NMAX
         IVARB = ISTBNU(N)
         NDEG  = 0
         DO 70 LA = 0,MAXLOP
            IF (IBTAND(IVARB,LA) .EQ. 0) NDEG = NDEG + 1
70       CONTINUE
         DEGEN(N) = NDEG
         IF (IPRINT .GE. 5)
     *      WRITE (6,127) NAMN(N),(CORD(N,I),I=1,3),
     *                        IFAKBN(ISTBNU(N)), NDEG - 1
127      FORMAT(//,2X,A4,5X,3F10.5,6X,/,
     *          16X, 'Stabilizer ',I3,',   with ',I2,
     *          ' symmetry equivalent atoms',/)
         DO 25 M = N, NMAX
            IVARA = ISTBNU(M)
            DO 26 KB = 0, MAXLOP
               IF (IBTAND(KB,IVARA) .EQ. 0) THEN
                  IF(M .EQ. N .AND. KB .EQ. 0) GO TO 26
                  IVARC = IBTAND(IVARA,IVARB)
                  HKAB  = FMULT(IVARC)
                  IF (M .EQ. N) HKAB = 0.5*HKAB
                  DIST2 = D0
                  DO 27 I = 1, 3
                     IVARC = ISYTYP(I)
                     TERM  = (CORD(N,I) - CORD(M,I)*
     *                           PT(IBTAND(IVARC,KB)))**2
                     DIST2 = DIST2 + TERM
   27             CONTINUE
                  DIST = SQRT(DIST2)
                  BDIS = BOHRAD*DIST
                  IF (DIST.LT.0.1.AND.CHARGE(M)*CHARGE(N).NE.D0) THEN
                     WRITE (6,29) N, KB, M, SQRT(DIST2)
   29                FORMAT (' The nucleus',I3,' is too close to the',
     *                      I3,'''th. transformation of nucleus No.',I3,
     *                      ' DIST=',D15.5)
                     CALL ERREX
                  END IF
                  IF (IPRINT .GE. 5 .AND. M .LE. MMAX)
     *               WRITE (6,128) NAMN(M),KB,BDIS
  128             FORMAT (15X,A4,I4,5X,'Dist=',F6.3)
                  IF (IBTAND(KB,IVARB) .EQ. 0) THEN
C                    IF(IFBA(M)+IFBA(N).EQ.0) GO TO 26
                     IF((ABS(CHARGE(M)) .GT. 100.0D0) .AND.
     *                  (ABS(CHARGE(N)) .GT. 100.0D0)) GOTO 26
                     IF(DIST.EQ.D0) GOTO 26
                     POTNUC=POTNUC+CHARGE(M)*CHARGE(N)*HKAB/DIST
                  END IF
               END IF
   26       CONTINUE
   25    CONTINUE
 8724 CONTINUE
C
C     Set parameters in COMMON /NUCLEI/
C
      NUCDEP    = 0
      NUCPRE(1) = 0
      NUCIND    = NMAX
      NFLOAT    = 0
      DO 510 JATOM = 1,NUCIND
         IDEG          = DEGEN(JATOM)
         NUCDEG(JATOM) = IDEG
         NUCPRE(JATOM + 1) = NUCPRE(JATOM) + IDEG
         IF (NINT(CHARGE(JATOM)) .EQ. 0) NFLOAT = NFLOAT + 1
         ICENT = NUCPRE(JATOM)
         MULX  = ISTBNU(JATOM)
         DO 512 ISYMOP = 0,MAXLOP
            IF (IBTAND(MULX,ISYMOP) .EQ. 0) THEN
               NUCDEP = NUCDEP + 1
               ICENT  = ICENT  + 1
               NUCNUM(JATOM,ISYMOP+1) = ICENT
            ELSE
               NUCNUM(JATOM,ISYMOP+1) = 0
            END IF
  512    CONTINUE
  510 CONTINUE
      NATOMS = NUCDEP - NFLOAT
C
C     ***** Determine IAOAO *****
C
      CALL AOTOAO(IAOAO,IPRINT)
C
C     **** Transformation matrix Cartesian to symmetry coordinates ****
C
      NCOOR=3*NUCDEP
      IF(JFC) NCOOR=NUCDEP
      IF(JSD) NCOOR=6*NUCDEP
      KCSTRA=1
      KSCTRA=KCSTRA+NCOOR*NCOOR
      KDWORK=KSCTRA+NCOOR*NCOOR
      KIWORK=KDWORK+NCOOR
      CALL TRACOR(WORK1(KCSTRA),WORK1(KSCTRA),WORK1(KDWORK),
     &            WORK1(KIWORK),NCOOR,IPRINT)
C
C PRINT SECTION
C
C GROUP MULTIPLICATION TABLE
C
      IF (MAXLOP .GT. 0 .AND. IPRINT .GE. 2) THEN
         WRITE (6,'(/)')
         CALL HEADER('Group Multiplication Table',1)
         DO 1113 I=0,MAXLOP
            WRITE (6,'(1X,8I4)') (IBTXOR(I,J)+1,J=0,MAXLOP)
 1113    CONTINUE
      END IF
C
C     ***** Printing of atomic coordinates *****
C
CJDW 4/1/97. Always call PRICAR --- otherwise the GRD file is not created
C            and we have troubles.
C
c     IF (LPRINT) CALL PRICAR
      CALL PRICAR
C
      IF (IPRINT .GE. 10 .OR. WRTONE) THEN
         WRITE (6,'(////10X,A,1P,D13.7////)')
     *       'Nuclear repulsion energy : ', POTNUC
      END IF
C
      IBX = 0
      DO 1127 I = 1, MAXLOP+1
         IBX = IBX + NAOS(I)
 1127 CONTINUE
#ifdef _OLDECP
C
CCH---------------------------------------------------------------------
C     ECP: PREPARE INFORMATION FOR THE ECPSTUFF
CCH---------------------------------------------------------------------
      IF (IECP.GT.0) THEN
        IERR=0
        CALL CHRGES(CHARGE,MTYPE,NDI10,IERR)
        IF (IERR.EQ.0) THEN
          IPR=0
          IUNIT=55
          INFILE='ECPDATA'
          CALL ECPDEF(NMAX,IPSEUX,MTYPE,CNICK)
          IF (IERR.EQ.0) THEN
            IPR=0
            NATOM=NMAX
            CALL ECPPAR(IPR,IPSEUX,CNICK,CHARGE,NATOM,.FALSE.,0)
          ELSE
            WRITE(*,*) 'E R R O R in subroutine ecppar!'
            WRITE(*,*) 'Please check your ecp-inputfile ECPDATA.'
            CALL ERREX
          ENDIF
        ELSE
          WRITE(*,*) 'E R R O R in subroutine chrges!'
          WRITE(*,*)
          CALL ERREX
        ENDIF
        CALL DOITGR(NAMAT,IQMSTR,JCOSTR,NUCSTR,NRCSTR,NHARM)
      ENDIF
CCH---------------------------------------------------------------------
C----- New ecplibray
#else
      Call Ecp_set4_grdint(Iecp,Nmax,namat,iqmstr,jcostr,nucstr,nrcstr)
#endif

      RETURN
C
 3000 WRITE (6,'(//A)') ' *** ERROR *** unable to open MOLECULE.INP'
      CALL ERREX
      END
