      SUBROUTINE READIN(I2,ITRAN,CUN,CONTUN,CTEMP,NU,NR,AT,KABS,
     &                  I1,NXYZ,NUC,NRC,NBAS,CORR,KA,C,JCO,PT,
     &                  TITLE,JTRAN,A,CTRAN,ICDS,NDEG,CHARGT,X,
CSSS     &                  IAD,NFUNMX,KWD,NH4,JtranC,ITranC,
     &                  IAD,NFUNMX,KWD,JtranC,ITranC,
     &                  CTranC,NSYMOP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cch---------------------------------------------------------------------
c     parameter for ecp (modified for Aces2!)
c-----------------------------------------------------------------------
C#include "ecplims.par"
c-----------------------------------------------------------------------
c     ndi10 : max. number of atoms
c     ndilmx: max. l-quantum number
cch---------------------------------------------------------------------
C-----------------------------------------------------------------------
C     (Other) Parameters
C-----------------------------------------------------------------------
C
#include "ecp.par"
#include "baslims.par"
#include "flags.h"
C
      PARAMETER (NPRIMX = 10)
      PARAMETER (NCONTS = 15, NCONTP = 8, NCONTD = 4)
      PARAMETER (NCONTF = 2, NCONTG = 1)
      PARAMETER (NCONTH = 1, NCONTI = 1)
C-----------------------------------------------------------------------
      INTEGER  AIQM
      INTEGER  IPRX(MXCBF)
      INTEGER  IAMN(MXCBF)
      INTEGER STARS, FORCIN, EODATA
      INTEGER  IAD(NHT)
      INTEGER  IFBA(1000)
      INTEGER  NCONTR(7)
C SG 8/25/98 atom name variables changed to be characters instead of integers
#ifdef _NAMN_THROUGH_CHRTMP
      CHARACTER*4 NAMN(MXCBF),NAME 
      INTEGER DIM_CHRTMP
      PARAMETER (DIM_CHRTMP=MXCBF*4/10)
      CHARACTER*(DIM_CHRTMP) CHRTMP
      CHARACTER*1 CHRNUM
#else
      INTEGER NAMN(MXCBF)
#endif /* _NAMN_THROUGH_CHRTMP */
      CHARACTER*4 LAMN(MXTNSH),MAMN(MXCBF),LAM,LAM1
C
      character*1 ioverp
      character*6 kaset
      character*6 slask
      character*8 hblnk
      character*8 popinp,labinp,polinp,syminp,CSymIN(2)
      character*8 blankf
      character*8 nanum(8)
      character*1 ispd(7)
      character*1 ko(3)
      character*1 kka(8,3) 
      character*1 iblank
      character*2 koko(3)
      character*4 kwo(84)
      character*4 jprx(MXCBF)
      character*4 jjprx(MXCBF)
      character*4 nako 
      character*4 kktyp(84,2)
      character*6 iblnk
cch
      character*80 fname
C
C....    INPUT PROCESSING ROUTINE FOR  M O L E C U L E.
C....    JAN ALML0F, STOCKHOLM, DEC. 1971.
C
      INTEGER   AND,OR,EOR,DSTRT
      LOGICAL IFIF,LPRT,LHARM,DELE
      DIMENSION POLINP(2), SYMINP(2), BLANKF(6), POPINP(2), LABINP(2) 
      DIMENSION NXBAS(8),JCNTAM(28,7),KDELB(200)
      DIMENSION I2(21), ITRAN(MXCBF,224), CUN(MXTNCC),
     1            CONTUN(MXTNCC), CTEMP(NFUNMX,NFUNMX), 
     2            NU(NFUNMX), NR(NFUNMX), AT(NFUNMX),
     3            KABS(MXSHEL,28,30,8)
      dimension iadq(8),charm(1596)
C SG 6/25/97 Changed to agree with Stanton's code to use 500 functions
      dimension last(MXCBF),iast(MXCBF)
      dimension nast(MXCBF),jaast(MXCBF)
C
C     DIMENSION I1(50),NXYZ(KWD,3),NUC(NHT,25),NRC(NHT,25),
      DIMENSION I1(50),NXYZ(KWD,3),NUC(NHT,NFUNMX),NRC(NHT,NFUNMX),
     1 NBAS(8),
     2 CORR(100*3),KA(8,3),C(MXTNCC),JCO(NHT),PT(8),TITLE(18),
     3 JTRAN(MxCBF),A(50),CTRAN(MxCBF,224),
     4 ICDS(MXTNPR), NDEG(100), CHARGT(100), X(3)
C
C     To hold SOAO matrix pointers for cartesian basis
C
      Integer JTranC(MxCBF), ITranC(MxCBF, 224)
      INTEGER FILEOUT,FILEOUT2 
      Double precision CTranC(MxCBF, 224)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD 
      COMMON /REP/ NEWIND(MXCBF) , MSTOLD(8)
      COMMON /INDX/ PC(512),DSTRT(8,MXCBF),NTAP,LU2,NRSS,NUCZ,ITAG,
     1 MAXLOP,MAXLOT,KMAX,NMAX,KHKT(7),MULT(8),ISYTYP(3),ITYPE(7,28),
     2 AND(8,8),OR(8,8),EOR(8,8),NPARSU(8),NPAR(8),MULNUC(MXatms),
     3 NHKT(MXTNSH),MUL(MXTNSH),NUCO(MXTNSH),NRCO(MXTNSH),JSTRT(MXTNSH),
     4 NSTRT(MXTNSH),MST(MXTNSH),JRS(MXTNSH)
      COMMON /DAT/  ALPHA(MXTNPR),CONT(MXTNCC),CENT(3,MXTNSH),
     1              CORD(Mxatms,3),CHARGE(Mxatms),FMULT(8),TLA, TLC
      COMMON /HARMON/ LHARM
cch      COMMON /ECPOPT/ IECP
cold  COMMON /SYMIND/ LABORB(400)
      COMMON /SYMIND/ LABORB(MXCBF)
      Common /VMTask/ ITask
#include "flags.com"
#include "flags2.com"
c
cch---------------------------------------------------------------------
c     variables for ecp(input: ecpdef,doit):
c-----------------------------------------------------------------------
      integer ipr,iunit,ierr
      dimension iatdef(Mxatms),ipseux(Mxatms),ipseud(Mxatms)
      dimension iqmstr(Mxatms),jcostr(Mxatms,Maxang)
      dimension nucstr(Mxatms,Maxang,Mxprim),
     &          nrcstr(Mxatms,Maxang,Mxprim)
      character*4 namat(Mxatms)
      character*80 infile
      character*80 cnick(Mxatms)
c     character*80 cnick(Mxatms),infile
      character mtype(Mxatms)*(2)
c-----------------------------------------------------------------------
c     meaning of variables :
c-----------------------------------------------------------------------
c     ipr     : debug level, here: ipr=0
c     iunit   : unit number for infile
c     ierr    : errorflag; ierr=0 means no error during the process
c     iatdef  : lowest index within a group of atoms with the same
c               attributes (attention : is set to -1 during attribute
c               list processing !)
c     ipseux  : pseudo potential type of atom
c               0  means : no pseudo potential or not yet assigned
c     ipseud  : number of types of ecps
c     iqmstr  : max. l-quantum number of a given shell
c     namat   : name of the atom
c     nucstr  : number of primitive exponents of a given shell
c     nrcstr  : number of blocks of contraction coefficients of a
c               given shell
c     cnick   : pseudo potential nicknames indexed by pseudo potential
c               types
c     infile  : name of the inputfile for ecp
c     mtype   : atomic identifier (starting with element symbol)
cch---------------------------------------------------------------------
      common /ecpsts/ iecp
cch---------------------------------------------------------------------
C
      DATA IZILCH /0/
      DATA IECPU /35/
      DATA INU /5/, IF10 /0/, IF18 /0/
      DATA HBLNK/'        '/
      DATA JCNTAM /
     & 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,1,0,0,0,0,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     & 0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1/
      DATA POLINP/'********','POLYAINP'/
      DATA SYMINP/'********','SYMTRANS'/
      DATA CSymIn/'********','CSYMTRAN'/
      DATA LABINP/'********','LABBASIS'/
      DATA POPINP/'********','MULPOPDA'/
      DATA BLANKF/6*'        '/
      DATA LUGRAD /36/, LUT /6/, LUAUX /62/, LUSCR /13/
C
      DATA NANUM /'       1','       2','       3','       4',
     1            '       5','       6','       7','       8'/
      DATA ISPD /'s', 'p','d','f','g','h','i'/
      DATA KO /'X','Y','Z'/   , IBLANK /' '/
      DATA KOKO /'YZ','XZ','XY'/
      DATA (KWO(I),I=1,56)
     O         /'S   ', 'X   ', 'Y   ', 'Z   ',
     1          'XX  ', 'XY  ', 'XZ  ', 'YY  ',
     2          'YZ  ', 'ZZ  ', 'F300', 'F210',
     3          'F201', 'F120', 'F111', 'F102',
     4          'F030', 'F021', 'F012', 'F003',
     5          'G400', 'G310', 'G301', 'G220',
     6          'G211', 'G202', 'G130', 'G121',
     7          'G112', 'G103', 'G040', 'G031',
     8          'G022', 'G013', 'G004', 'H500',
     9          'H410', 'H401', 'H320', 'H311',
     O          'H302', 'H230', 'H221', 'H212',
     1          'H203', 'H140', 'H131', 'H122',
     2          'H113', 'H104', 'H050', 'H041',
     3          'H032', 'H023', 'H014', 'H005'/
      data (kwo(i),i=57,84) /
     $ 'I600', 'I510', 'I501', 'I420', 'I411', 'I402',
     $ 'I330', 'I321', 'I312', 'I303', 'I240', 'I231',
     $ 'I222', 'I213', 'I204', 'I150', 'I141', 'I132',
     $ 'I123', 'I114', 'I105', 'I060', 'I051', 'I042',
     $ 'I033', 'I024', 'I015', 'I006'/
      DATA NCONTR /NCONTS, NCONTP, NCONTD, NCONTF, NCONTG,
     $ NCONTH, NCONTI/
      data iadq/0,1,10,46,146,371,812,1596/
      data (charm(i),i = 1,146) /1,  1,0,0, 0,1,0, 0,0,1,
     D       -1., 0 , 0 ,-1., 0 , 2.,
     D        0 , 1., 0 , 0 , 0 , 0 ,
     D        0 , 0 , 1., 0 , 0 , 0 ,
     D        1., 0 , 0 ,-1., 0 , 0 ,
     D        0 , 0 , 0 , 0 , 1., 0 ,
     D        1., 0 , 0 , 1., 0 , 1.,
     F       -1., 0 , 0 ,-1., 0 , 4 , 0 , 0 , 0 , 0 ,
     F        0 ,-1., 0 , 0 , 0 , 0 ,-1., 0 , 4 , 0 ,
     F        0 , 0 ,-3., 0 , 0 , 0 , 0 ,-3., 0 , 2.,
     F        1., 0 , 0 ,-3., 0 , 0 , 0 , 0 , 0 , 0 ,
     F        0 , 0 , 0 , 0 , 1., 0 , 0 , 0 , 0 , 0 ,
     F        1., 0 , 0 , 1., 0 , 1., 0 , 0 , 0 , 0 ,
     F        0 , 3., 0 , 0 , 0 , 0 ,-1., 0 , 0 , 0 ,
     F        0 , 0 , 1., 0 , 0 , 0 , 0 ,-1., 0 , 0 ,
     F        0 , 1., 0 , 0 , 0 , 0 , 1., 0 , 1., 0 ,
     F        0 , 0,  1., 0 , 0 , 0 , 0 , 1., 0 , 1./
      data (charm(i),i = 147,371) /
     G  3., 0 , 0 , 6 , 0 ,-24, 0 , 0 , 0 , 0 , 3., 0 ,-24, 0 , 8 ,
     G  0 ,-1., 0 , 0 , 0 , 0 ,-1., 0 , 6 , 0 , 0 , 0 , 0 , 0 , 0 ,
     G  0 , 0 ,-3., 0 , 0 , 0 , 0 ,-3., 0 , 4 , 0 , 0 , 0 , 0 , 0 ,
     G  1., 0 , 0 ,-6 , 0 , 0 , 0 , 0 , 0 , 0 , 1., 0 , 0 , 0 , 0 ,
     G  0 , 0 , 0 , 0 , 3., 0 , 0 , 0 , 0 , 0 , 0 ,-1., 0 , 0 , 0 ,
     G -1., 0 , 0 , 0 , 0 , 6 , 0 , 0 , 0 , 0 , 1., 0 ,-6 , 0 , 0 ,
     G  0 , 1., 0 , 0 , 0 , 0 ,-1., 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
     G  0 , 0 , 1., 0 , 0 , 0 , 0 ,-3., 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
     G  0 , 1., 0 , 0 , 0 , 0 , 1., 0 , 1., 0 , 0 , 0 , 0 , 0 , 0 ,
     G  0 , 0 , 1., 0 , 0 , 0 , 0 , 1., 0 , 1., 0 , 0 , 0 , 0 , 0 ,
     G  1., 0 , 0 , 2., 0 , 2., 0 , 0 , 0 , 0 , 1., 0 , 2., 0 , 1.,
     G  0 , 0 , 0 , 0 ,-3., 0 , 0 , 0 , 0 , 0 , 0 ,-3., 0 , 4 , 0 ,
     G  1., 0 , 0 , 0 , 0 , 1., 0 , 0 , 0 , 0 ,-1., 0 ,-1., 0 , 0 ,
     G  0 , 0 , 0 , 0 , 1., 0 , 0 , 0 , 0 , 0 , 0 , 1., 0 , 1., 0 ,
     G -1., 0 , 0 ,-2., 0 , 1., 0 , 0 , 0 , 0 ,-1., 0 , 1., 0 , 2./
      DATA (charm(I),I = 372,539) /
     H  1.,0.,0.,2.,0.,-12.,0.,0.,0.,0.,1.,0.,
     H    -12.,0.,8.,0.,0.,0.,0.,0.,0.,
     H  0.,1.,0.,0.,0.,0.,2.,0.,-12.,0.,0.,0.,
     H     0.,0.,0.,1.,0.,-12.,0.,8.,0.,
     H  0.,0.,-1.,0.,0.,0.,0.,0.,0.,2.,0.,0.,
     H     0.,0.,0.,0.,1.,0.,-2.,0.,0.,
     H -1.,0.,0.,2.,0.,8.,0.,0.,0.,0.,3.,0.,
     H    -24.,0.,0.,0.,0.,0.,0.,0.,0.,
     H  0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,-1.,
     H     0.,0.,0.,0.,0.,0.,0.,0.,0.,
     H  1.,0.,0.,2.,0.,2.,0.,0.,0.,0.,1.,0.,
     H     2.,0.,1.,0.,0.,0.,0.,0.,0.,
     H  0.,5.,0.,0.,0.,0.,-10.,0.,0.,0.,0.,0.,
     H     0.,0.,0.,1.,0.,0.,0.,0.,0.,
     H  0.,0.,1.,0.,0.,0.,0.,-6.,0.,0.,0.,0.,
     H     0.,0.,0.,0.,1.,0.,0.,0.,0. /
      DATA (charm(I),I = 540,707) /
     H  0.,-3.,0.,0.,0.,0.,-2.,0.,24.,0.,0.,0.,
     H     0.,0.,0.,1.,0.,-8.,0.,0.,0.,
     H  0.,0.,15.,0.,0.,0.,0.,30.,0.,-40.,0.,0.,
     H     0.,0.,0.,0.,15.,0.,-40.,0.,8.,
     H  1.,0.,0,-10.,0.,0.,0.,0.,0.,0.,5.,0.,
     H     0.,0.,0.,0.,0.,0.,0.,0.,0.,
     H  0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,1.,
     H     0.,1.,0.,0.,0.,0.,0.,0.,0.,
     H  1.,0.,0.,-2.,0.,1.,0.,0.,0.,0.,-3.,0.,
     H    -3.,0.,0.,0.,0.,0.,0.,0.,0.,
     H  0.,0.,0.,0.,-1.,0.,0.,0.,0.,0.,0.,-1.,
     H     0.,2.,0.,0.,0.,0.,0.,0.,0.,
     H -1.,0.,0.,-2.,0.,3.,0.,0.,0.,0.,-1.,0.,
     H     3.,0.,4.,0.,0.,0.,0.,0.,0.,
     H  0.,1.,0.,0.,0.,0.,2.,0.,2.,0.,0.,0.,
     H     0.,0.,0.,1.,0.,2.,0.,1.,0. /
      DATA (charm(I),I = 708,812) /
     H  0.,0.,1.,0.,0.,0.,0.,2.,0.,2.,0.,0.,
     H     0.,0.,0.,0.,1.,0.,2.,0.,1.,
     H  0.,3.,0.,0.,0.,0.,2.,0.,3.,0.,0.,0.,
     H     0.,0.,0.,-1.,0.,-1.,0.,0.,0.,
     H  0.,0.,1.,0.,0.,0.,0.,0.,0.,1.,0.,0.,
     H     0.,0.,0.,0.,-1.,0.,-1.,0.,0.,
     H  0.,-1.,0.,0.,0.,0.,-2.,0.,3.,0.,0.,0.,
     H     0.,0.,0.,-1.,0.,3.,0.,4.,0.,
     H  0.,0.,-3.,0.,0.,0.,0.,-6.,0.,-1.,0.,0.,
     H     0.,0.,0.,0.,-3.,0.,-1.,0.,2. /
      DATA (charm(I),I = 813,980) /
     I   1.,   0.,   0., -15.,   0.,   0.,   0.,   0.,   0.,   0.,
     I  15.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  -1.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   1.,   0.,   0.,   0.,   0.,   2.,   0., -16.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   1.,   0., -16.,   0.,  16.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,   1.,   0.,   0.,   0.,   0., -10.,   0.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   5.,   0.,   0.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I  -1.,   0.,   0.,   5.,   0.,  10.,   0.,   0.,   0.,   0.,
     I   5.,   0., -60.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  -1.,   0.,  10.,   0.,   0.,   0.,   0.,
     I   0.,   0.,   0.,   0.,   5.,   0.,   0.,   0.,   0.,   0.,
     I   0., -10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,   1.,   0.,   0.,   0.,   0.,   0.,
     I   1.,   0.,   0.,   1.,   0., -16.,   0.,   0.,   0.,   0.,
     I  -1.,   0.,   0.,   0.,  16.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  -1.,   0.,  16.,   0., -16.,   0.,   0. /
      DATA (charm(I),I = 981,1148) /
     I   0.,   6.,   0.,   0.,   0.,   0., -20.,   0.,   0.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   6.,   0.,   0.,   0.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,  -3.,   0.,   0.,   0.,   0.,   6.,   0.,   8.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   9.,   0., -24.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  -4.,   0.,   0.,   0.,   0.,   0.,   0.,  40.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   4.,   0., -40.,   0.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,   1.,   0.,   0.,   0.,   0.,  -2.,   0.,   1.,
     I   0.,   0.,   0.,   0.,   0.,   0.,  -3.,   0.,  -3.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I  -5.,   0.,   0., -15.,   0.,  90.,   0.,   0.,   0.,   0.,
     I -15.,   0., 180.,   0.,-120.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  -5.,   0.,  90.,   0.,-120.,   0.,  16.,
     I   0.,   0.,   0.,   0.,  -9.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  -6.,   0.,  24.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,   3.,   0.,  -8.,   0.,   0.,   0. /
      DATA (charm(I),I = 1149,1316) /
     I   1.,   0.,   0.,  -5.,   0.,   1.,   0.,   0.,   0.,   0.,
     I  -5.,   0.,  -6.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   1.,   0.,   1.,   0.,   0.,   0.,   0.,
     I   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  20.,   0., -40.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,  10.,   0., -40.,   0.,  16.,   0.,
     I  -1.,   0.,   0.,  -1.,   0.,   5.,   0.,   0.,   0.,   0.,
     I   1.,   0.,   0.,   0.,   6.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   1.,   0.,  -5.,   0.,  -6.,   0.,   0.,
     I   0.,  -1.,   0.,   0.,   0.,   0.,  -2.,   0.,   5.,   0.,
     I   0.,   0.,   0.,   0.,   0.,  -1.,   0.,   5.,   0.,   6.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,  10.,   0.,   0.,   0.,   0.,  20.,   0., -40.,
     I   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0., -40.,   0.,
     I  16.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   1.,   0.,   0.,   0.,   0.,   0.,   0.,   1.,   0.,
     I   0.,   0.,   0.,   0.,   0.,  -1.,   0.,  -1.,   0.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0. /
      DATA (charm(I),I = 1317,1484) /
     I   0.,   0.,  -3.,   0.,   0.,   0.,   0.,  -6.,   0.,   1.,
     I   0.,   0.,   0.,   0.,   0.,   0.,  -3.,   0.,   1.,   0.,
     I   4.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   1.,   0.,   0.,   0.,   0.,   2.,   0.,   2.,   0.,
     I   0.,   0.,   0.,   0.,   0.,   1.,   0.,   2.,   0.,   1.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,   1.,   0.,   0.,   0.,   0.,   2.,   0.,   2.,
     I   0.,   0.,   0.,   0.,   0.,   0.,   1.,   0.,   2.,   0.,
     I   1.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     I  -1.,   0.,   0.,  -3.,   0.,   0.,   0.,   0.,   0.,   0.,
     I  -3.,   0.,   0.,   0.,   3.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  -1.,   0.,   0.,   0.,   3.,   0.,   2.,
     I   0.,   0.,   0.,   0.,   3.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   2.,   0.,   3.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,  -1.,   0.,  -1.,   0.,   0.,   0.,
     I   1.,   0.,   0.,   1.,   0.,   2.,   0.,   0.,   0.,   0.,
     I  -1.,   0.,   0.,   0.,   1.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  -1.,   0.,  -2.,   0.,  -1.,   0.,   0. /
      DATA (charm(I),I = 1485,1596) /
     I   0.,   0.,   0.,   0.,  -3.,   0.,   0.,   0.,   0.,   0.,
     I   0.,  -6.,   0.,   1.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,  -3.,   0.,   1.,   0.,   4.,   0.,
     I   3.,   0.,   0.,   9.,   0., -21.,   0.,   0.,   0.,   0.,
     I   9.,   0., -42.,   0., -16.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   3.,   0., -21.,   0., -16.,   0.,   8.,
     I   0.,   0.,   0.,   0.,   1.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   2.,   0.,   2.,   0.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   0.,   1.,   0.,   2.,   0.,   1.,   0.,
     I   1.,   0.,   0.,   3.,   0.,   3.,   0.,   0.,   0.,   0.,
     I   3.,   0.,   6.,   0.,   3.,   0.,   0.,   0.,   0.,   0.,
     I   0.,   1.,   0.,   3.,   0.,   3.,   0.,   1. /
      DATA (KKTYP(I,1),I = 1,84) /'S000', 'P100', 'P010', 'P001',
     1            'D200', 'D110', 'D101', 'D020', 'D011', 'D002',
     2            'F300', 'F210', 'F201', 'F120', 'F111', 'F102',
     3            'F030', 'F021', 'F012', 'F003',
     4            'G400', 'G310', 'G301', 'G220', 'G211', 'G202',
     5            'G130', 'G121', 'G112', 'G130', 'G040', 'G031',
     6            'G022', 'G031', 'G004',
     7            'H500', 'H410', 'H401', 'H320', 'H311', 'H302',
     8            'H230', 'H221', 'H212', 'H203', 'H140', 'H131',
     9            'H122', 'H113', 'H104', 'H050', 'H041', 'H032',
     O            'H023', 'H014', 'H005',
     1            'I600', 'I510', 'I501', 'I420', 'I411', 'I402',
     1            'I330', 'I321', 'I312', 'I303', 'I240', 'I231',
     1            'I222', 'I213', 'I204', 'I150', 'I141', 'I132',
     1            'I123', 'I114', 'I105', 'I060', 'I051', 'I042',
     1            'I033', 'I024', 'I015', 'I006'/
      DATA (KKTYP(I,2),I = 1,84) /'1s  ', '2px ', '2py ', '2pz ',
     1            '3d0 ', '3d2-', '3d1+', '3d2+', '3d1-', '3s  ',
     2            '4f1+', '4f1-', '4f0 ', '4f3+', '4f2-', '4px ',
     3            '4f3-', '4f2+', '4py ', '4pz ',
     4            '5g0 ', '5g2-', '5g1+', '5g4+', '5g3-', '5g2+',
     5            '5g4-', '5g3+', '5d2-', '5d1+', '5s  ', '5g1-',
     6            '5d2+', '5d1-', '5d0' ,
     H '6h1+','6h1-','6h2+','6h3+','6h4-','6px ','6h5-','6h4+','6h3-',
     H   '6h0 ','6h5+','6f2-','6f3+','6h2-','6f1+','6py ','6pz ',
     H   '6f3-','6f2+','6f1-','6f0 ',
     I            '7i6+', '7i2-', '7i5+', '7i4+', '7i5-', '7i2+',
     I            '7i6-', '7i3+', '7i4-', '7g3+', '7i0 ', '7i3-',
     I            '7g4+', '7i1-', '7g2+', '7g2-', '7i1+', '7g4-',
     I            '7g1+', '7d2-', '7d1+', '7d0 ', '7g3-', '7d2+',
     I            '7g1-', '7g0 ', '7d1-', '7s  '/
      DATA IBLNK/'      '/
      DATA NCMAX /100/
C
C...  STATEMENT FUNCTIONS
C
      PTAND(I,J) = PC(AND(I,J))

#ifdef _TRACE_VMOL
      write(*,*) '@VMOL: Entered file ',__FILE__
#endif

      IPRINT = IFLAGS(1)
cch---------------------------------------------------------------------
c     init for ecp
c-----------------------------------------------------------------------
      call izero(iqmstr,ndi10)
      call izero(jcostr,ndi10*ndilmx)
      call izero(nucstr,ndi10*ndilmx)
      call izero(nrcstr,ndi10*ndilmx)
cch---------------------------------------------------------------------
C
C  OPEN INPUT FILE MOL
C
      OPEN(UNIT=INU,FILE='MOL',FORM='FORMATTED',STATUS='OLD')
      OPEN(UNIT=LUSCR,FILE='VMSCR1',FORM='UNFORMATTED',STATUS='UNKNOWN')
      FILEOUT=51
      FILEOUT2=52
      LUTR1=53
      OPEN(UNIT=FILEOUT,FILE='ACES2_BASINFO',status='UNKNOWN')
      OPEN(UNIT=FILEOUT2,FILE='ACES2_BASINFO2',status='UNKNOWN')
      OPEN(UNIT=LUTR1,FILE='ACES2_TR1',status='UNKNOWN')
C
CCH
C      IPRINT=1
CJDW 4/20/94. Let's go back to letting the user control how much vmol
C             printing he/she will have --- this might start making people
C             check their basis sets again.
C      IPRINT=0
CCH
      IECP=0
c      ISTATS=0
      ID3=0 
      TLA=1.D-10 
      SCALE=0.0
      JUNK=1
      NONT=1 
      TLC=0.D0
      NTAP=2
      DO 3335 I=1,3
       ISYTYP(I)=1
3335  CONTINUE
C
      NM=0
      DO 6 I=0,NHT
      DO 6 IX=I,0,-1
      DO 6 IY=I-IX,0,-1
      NM=NM+1
      NXYZ(NM,1)=IX
      NXYZ(NM,2)=IY
      NXYZ(NM,3)=I-IX-IY
    6 CONTINUE
C
      PI=2.D0*DACOS(0.0D00)
      PIPPI=(0.5D0/PI)**0.75D0
      write(*,*) 'One- and two-electron integrals over ',
     &           'symmetry-adapted AOs are calculated.\n'
C
C-----------------------------------------------------------------------
C     Begin reading the MOL file.
C-----------------------------------------------------------------------
C
      READ(INU,83383) SLASK,IGENER,IF10,IF18,IECP,IECPU,IJUNK,
     1                INUTMP,ISTATS,NHARM,IDELA,IDELB,IBORTH,ITask
83383 FORMAT(A6,4X,13I5)
83384 FORMAT(4X,A6,4X,12I5)
      IF (SLASK .EQ. 'SUPMAT') THEN
        WRITE(*,*) '\n\n'
        write(*,*) '* *  THIS CODE DOES NOT GENERATE SUPERMATRICES  * *'
        CALL ERREX
      ENDIF
c      IF(IECP.NE.0)IF18=1
      LHARM = NHARM.eq.0
      IF (LHARM) THEN
        IDOSPH = 1
        write(*,*) '@READIN: Spherical harmonics are used.'
      ELSE
        IDOSPH = 0
        write(*,*) '@READIN: Cartesian gaussians are used.'
        idela=1
      ENDIF
c      if (idela.eq.0) then
c         write(*,*) '  Contaminants will be deleted'
c      else
c         write(*,*) '  Contaminants will be retained'
c      end if
      if(iprint.ge.1)write(6,176)
      if (iecp.ne.0) then
         write(*,*) '@READIN: Effective core potentials are used.'
      end if
      IF (INUTMP .NE. 0) INU = INUTMP
      LPRT=.FALSE.
      IF(IPRINT.GE.1)LPRT = .TRUE.
      READ(INU,833) TITLE
  833 FORMAT (9A8)
C
C-----------------------------------------------------------------------
C     Read number of nuclear types (NONTYP), symmetry information
C     (NSYMOP, KKA), tolerances (TLA, TLC), and coordinate units (ID3).
C     IRSTRT is not used anywhere.
C-----------------------------------------------------------------------
      READ(INU,8800) NONTYP, NSYMOP, ((KKA(J,I),I=1,3),J=1,3),
     1              TLA,  TLC, ID3, IRSTRT
CJDW  We now set TLA in MKVMOL in xjoda (keyword INTGRL_TOL), which is
C     written to the MOL file..
C     TLA=1.D-14
8800  FORMAT(2I5,9A1,1X,E10.2,E10.2,I5,I5)
C
C-----------------------------------------------------------------------
C     Open integral files.
C-----------------------------------------------------------------------
C
      CALL GFNAME('IIII    ',FNAME,ILENGTH)
      OPEN(UNIT=2,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
      IF(NSYMOP.NE.0)THEN
       CALL GFNAME('IJIJ    ',FNAME,ILENGTH)
       OPEN(UNIT=81,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
       CALL GFNAME('IIJJ    ',FNAME,ILENGTH)
       OPEN(UNIT=82,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &      STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
       IF(NSYMOP.GT.1)THEN
        CALL GFNAME('IJKL    ',FNAME,ILENGTH)
        OPEN(UNIT=83,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &       STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
       ENDIF
      ENDIF
C
      IF (TLC.eq.0.0) TLC=0.1*TLA
      READ(INU,180)TIML,SAFL
180   FORMAT (2F10.5)
c      IF (IECP .NE. 0) INU = IECPU
C
C-----------------------------------------------------------------------
C     Process symmetry information and construct symmetry arrays.
C-----------------------------------------------------------------------
C
C  MAXLOP = NUMBER OF IRREDUCIBLE REPRESENTATIONS,
C  THE LOOPS ARE ALWAYS FROM 1 TO MAXLOP
C
C  ISYTYP DETERMINES THE BEHAVIOUR OF X,Y,Z UNDER THE SYMMETRY
C  OPERATIONS: BIT 1,2,3 ARE SET TO ONE IF X, OR Y OR Z ARE
C  CHANGING UNDER THIS PARTICULAR OPERATIONS.
C  !ATTENTION! ISYTYP IS INITIALIZED TO 1 AND NOT TO 0 IN VMOL.
C  THIS HAS TO BE CONSIDERED BY USING THE DEFINITION OF ISYTYP
C  GIVEN ABOVE!
C  EXAMPLE: ISYTYP(2)=2 MEANS THE VALUE OF Y CHANGES UNDER
C           1 SYMMETRYOPERATION.
C
      MAXLOP=2**NSYMOP
      if(iprint.ge.1)write(6,1010) nontyp
1010  format(' ',I3,' types of atoms')
      if(iprint.ge.1)write(6,1011) nsymop
1011  format(' ',I3,' symmetry operations',/)
CCH
CCH   GENERATORS:
CCH   NSYMOP: NUMBER OF GENERATORS,
CCH   MAXLOP: NUMBER OF SYMMETRYOPERATIONS (NUMBER OF IRREPS).
CCH
      IF(NSYMOP.GT.0) THEN
        DO 111 J=1,NSYMOP
          NBL=1
          DO 112 I=1,3
            IF (KKA(J,I).EQ.IBLANK) NBL=NBL+1
            DO 112 K=1,3
              IF (KKA(J,I).EQ.KO(K)) ISYTYP(K)=ISYTYP(K)+2**(J-1)
  112     CONTINUE
          IF (NBL.EQ.1) THEN
            if(iprint.ge.1)write(6,31)
   31       format('   Inversion centre')
          ELSE IF (NBL.EQ.2) THEN
            DO 42 K=1,3
              DO 41 I=1,3
                IF (KKA(J,I).EQ.KO(K)) GO TO 42
   41         CONTINUE
              if(iprint.ge.1)write(6,67) ko(k)
67            format('   Rotation about the ',a1,'-axis')
   42       CONTINUE
          ELSE
            DO 52 K=1,3
              DO 51 I=1,3
                IF (KKA(J,I).EQ.KO(K)) THEN
                  if(iprint.ge.1)write(6,77) koko(k)
77                format('   Reflection in the ',a2,'-plane')
                  GO TO 111
                ENDIF
   51         CONTINUE
   52       CONTINUE
          ENDIF
  111   CONTINUE
      ENDIF
      if(iprint.ge.1)then
      write(6,176)
176   format('  ')
      write(6,177) tla
177   format('  Integrals less than ',e9.2,' are neglected.',/)
      endif
C
C     Determine KA(I,J), PT(I), MULT(I) and FMULT(I)
CCH   KA, PT: fixed
CCH   MULT (FMULT): depend only on the value of NSYMOP
C
      DO 101 I=1,8
        IS=0
        DO 107 K=1,3
          KA(I,K)=(I-1)/(2**(K-1))-2*((I-1)/(2**K))
  107   IS=IS+KA(I,K)
        PT(I)=(-1.0)**IS
        MULT(I)=2**MAX0(0,NSYMOP-IS)
  101 FMULT(I)=DFLOAT(MULT(I))
C
CCH   CREATE THE AND, OR AND EOR-ARRAYS
C
      DO 108 I=1,8
        DO 108 J=1,8
        IAN=1
        IOR=1
        DO 102 K=1,3
          IAN=IAN+MIN0(KA(I,K),KA(J,K))*2**(K-1)
  102   IOR=IOR+MAX0(KA(I,K),KA(J,K))*2**(K-1)
        AND(I,J)=IAN
        OR(I,J)=IOR
        EOR(I,J)=IOR-IAN+1
  108 CONTINUE
C
CCH   CREATE ITYPE-ARRAY
CCH   ITYPE(1,1)=1: BEHAVIOUR OF R^0 (NO CHANGE OF SIGN AT ALL!)
CCH   ITYPE(2,1..3)=ISYTYP(1..3): BEHAVIOUR OF R^1 (X, Y, Z)
CCH   ...AND SO ON UP TO R^6.
C
      ITYPE(1,1)=1
      NDO=0
      NFO=0
      NGO=0
      NHO=0
      NIO=0
      DO 105 JQM=1,3
        ITYPE(2,JQM)=ISYTYP(JQM)
        DO 105 KQM=JQM,3
          NDO=NDO+1
          ITYPE(3,NDO) = EOR(ISYTYP(JQM),ISYTYP(KQM))
          DO 105 LQM=KQM,3
            NFO=NFO+1
            ITYPE(4,NFO)=EOR(ITYPE(3,NDO),ISYTYP(LQM))
            DO 105 M2M=LQM,3
              NGO=NGO+1
              ITYPE(5,NGO)=EOR(ITYPE(4,NFO),ISYTYP(M2M))
              DO 105 M3M=M2M,3
                NHO=NHO+1
                ITYPE(6,NHO)=EOR(ITYPE(5,NGO),ISYTYP(M3M))
                DO 105 M4M=M3M,3
                  NIO=NIO+1
                  ITYPE(7,NIO)=EOR(ITYPE(6,NHO),ISYTYP(M4M))
  105 CONTINUE
C
CCH   KHKT: ENTARTUNG OF L-QUANTUM NUMBERS 1-7
      DO 106 JJ=1,7
  106 KHKT(JJ)=JJ*(JJ+1)/2
C
C     RUN OVER TYPES OF ATOMS IN THE MOLECULE
C
      KMAX=0
      LMAX=0
      LLMAX = 0
      LLMAX = 0
      NMAX=0
      MAXLO=8
      NHTYP=2
C-----------------------------------------------------------------------
C     Begin loop over center (nuclear) types.
C-----------------------------------------------------------------------
      DO 5 I=1,NONTYP
      DO 900 J = 1,5
        IAD(J) = 0
900   CONTINUE
      READ(INU,11) KASET, Q,NONT,IQM,(JCO(J),J=1,IQM), IGENES
C
C IQM is the maximum angular momentum for a given center. NO analytical
C gradinets for > g functions. Lets trap it here since there is
C no point continuing until vdint find this out. Ajith Perera, 06/2006.
C This restriction is eliminated, 04/2013, Ajith Perera

CSSS      If (Iflags2(h_IFLAGS2_grad_calc) .NE.0 .AND. Iqm .GE. 5 .AND.
CSSS     &    Iflags2(h_IFLAGS2_geom_opt)  .NE.0 .AND. 
CSSS     &    Iflags(h_IFLAGS_vibration)  .NE.0) Then
CSSS          Print*, "Analytical gradients are not yet implemented for"
CSSS          Print*, "h or higher angular momentum basis functons "
CSSS          Call Errex
CSSS      Endif
c
c      WRITE(LUT,12861) I
c12861 FORMAT(////20X,'ATOMIC TYPE NUMBERu
c      WRITE(LUT,12861) I
c12861 FORMAT(////20X,'ATOMIC TYPE NUMBER ',I2,/)
      if(iprint.ge.1)write(6,12861) I
12861  format('  Atomic type number ',I2,' :',/
     &        '  -----------------------',/)
      LUBAS = INU
11    FORMAT (A6,F14.8,12I5)
      if(iprint.ge.1)then
       write(6,12862) INT(Q)
12862  format('    Nuclear charge:                       ',I2)
       write(6,12863) nont
12863  format('    number of symmetry independent atoms: ',I2)
       write(6,12864) ispd(iqm)
12864  format('    highest orbital type:                  ',A1,/)
      endif
      IF(IQM .EQ. 0) GOTO 920
      if(iprint.ge.1)then
      IF (IGENER .EQ. 0 .AND. IGENES .EQ. 0) THEN
        IF (KASET .EQ. '      ') THEN
          WRITE(LUT,1287) (JCO(J),ISPD(J),J=1,IQM)
        ELSE
          WRITE(LUT,1278) (JCO(J), IAD(J), ISPD(J),J = 1,IQM)
          IF (SCALE .NE. 0. .AND. LPRT) WRITE(LUT,12781) SCALE
        ENDIF
      ELSE
        IF (KASET .EQ. '      ') THEN
          write(6,81287) (jco(j),ispd(j),j=1,iqm)
81287     format(4x,I3,' groups of CGTOs of ',a1,' type')
        ELSE
          WRITE(LUT,81278) (JCO(J), IAD(J), ISPD(J),J = 1,IQM)
          IF (SCALE .NE. 0. .AND. LPRT) WRITE(LUT,12781) SCALE
        ENDIF
      ENDIF
      endif
1287  FORMAT(15X,I3,' CGTOS OF  ',A1,' TYPE')
1278  FORMAT(15X,I3,' ORIGINAL AND ',I2,'  ADDED CGTOS OF ',A1,' TYPE')
81278 FORMAT(15X,I3,' ORIGINAL AND ',I2,'  ADDED GROUPS OF CGTOS OF ',
     1           A1,' TYPE')
12781 FORMAT(/,' BASIS LIBRARY S FUNCTION EXPONENTS WILL BE SCALED',
     1         ' USING A SLATER EXPONENT OF ',F12.8,/)
920   CONTINUE
      IQM4 = MAX0(IQM,1)
C      WRITE(LUGRAD) Q, NONT, IQM, (JCO(J),J = 1,IQM4)
C
C....    SET TO READ GEOMETRY TO HIGHER PRECISION IF REQUIRED 
C
      if(iprint.ge.1)write(6,1293)
1293  format(/,'    nuclear coordinates (in a.u.) are :',/)
      NMIN=NMAX+1
C
CCH   LOOP OVER SYMMETRY INDEPENDENT ATOMS IN THE MOLECULE
CCH   N=1,NONT
C-----------------------------------------------------------------------
C     Begin loop over symmetry inequivalent centers of each type,
C     reading coordinates.
C-----------------------------------------------------------------------
C
      DO 114 N=1,NONT
CCH   NMAX IS THE TOTAL NUMBER OF SYMMETRY INDEPENDENT ATOMS
      NMAX = NMAX + 1
      IF (NMAX .LE. NCMAX) GO TO 113
      WRITE(6,1305) NCMAX
1305  FORMAT ('1',5X,'SPACE LIMITATION ENCOUNTERED IN READIN'
     1/11X,'NUMBER OF CENTERS EXCEEDS',I4)
      CALL ERREX
113   IFBA(NMAX) = IQM

      READ(INU,'(A4,3F20.12,I5)') NAMN(NMAX),(CORD(NMAX,J),J=1,3),IDX
#ifdef _DEBUG
      write(*,'(A6,a4,3f13.9)') 'NAMN: ',namn(nmax),(cord(nmax,j),j=1,3)
#endif /* _DEBUG */
      IF (IDX.NE.0) READ(INU,'(3F20.5)') (CORD(NMAX,J),J=1,3)

cch---------------------------------------------------------------------
c     ecp: get IQM for later calculations in doit
c-----------------------------------------------------------------------
c     save information for ecpcalculation
c     iqmstr: array of the highest l-quantum numbers of each atom
c-----------------------------------------------------------------------
        if (iecp.ne.0) then
          iqmstr(nmax)=IQM
          do 9234 izahl=1,iqm
            jcostr(nmax,izahl)=jco(izahl)
 9234     continue
c-----------------------------------------------------------------------
c     ecp: get atomsymbol
c-----------------------------------------------------------------------
C    warning, machine dependent code
C
C    namu(nmax)    Character string    nmax   integer
C                   --->   sometimes there are problems
C
          namat(nmax)=namn(nmax)(1:2)
        endif
cch---------------------------------------------------------------------
       if(iprint.ge.1)write(6,176)
C
      IF (ID3 .EQ. 1) THEN
        DO 1148 J=1,3
        CORD(NMAX,J)=CORD(NMAX,J)/0.5291771
1148    CONTINUE
      ENDIF
114   CONTINUE
C-----------------------------------------------------------------------
C     End loop over symmetry inequivalent centers of each type.
C-----------------------------------------------------------------------
      IMAX=0
      IIMAX = 0
      IF(IQM.EQ.0) GO TO 1272
C-----------------------------------------------------------------------
C     Loop over angular momentum for each nuclear type.
C-----------------------------------------------------------------------
      DO 2 J=1,IQM
      KM=JCO(J) + IAD(J)
13    FORMAT (10I5)
211   FORMAT(8F10.4)
213   FORMAT(4F18.4)
14    FORMAT (8F9.4)
C
C
      IF (JCO(J) .EQ. 0 .AND. IAD(J) .EQ. 0) GOTO 2
C-----------------------------------------------------------------------
C     Loop over shells of each angular momentum. Read number of primit-
C     ives and contracted functions.
C-----------------------------------------------------------------------
      DO 88130 K = 1, KM
      IF (IGENER .EQ. 0 .AND. IGENES .EQ. 0) THEN
        READ(LUSCR) NUC(J,K), NRC(J,K)
      ELSE
        IF (K .LE. JCO(J)) THEN
          READ(LUBAS,13) NUC(J,K), NRC(J,K)
        ELSE
          READ(INU,13) NUC(J,K), NRC(J,K)
        ENDIF
      ENDIF
C
       IF(NUC(J,K) .GT. NFUNMX)THEN
        WRITE(6,*) '  @READIN-F, Too many primitives in this shell '
        WRITE(6,*) '             Value is ',NUC(J,K)
        WRITE(6,*) '             Limit is ',NFUNMX
        CALL ERREX
       ENDIF
C
       IF(NRC(J,K) .GT. NFUNMX)THEN
        WRITE(6,*) '  @READIN-F, Too many cont. fxns. in this shell '
        WRITE(6,*) '             Value is ',NRC(J,K)
        WRITE(6,*) '             Limit is ',NFUNMX
        CALL ERREX
       ENDIF
cch---------------------------------------------------------------------
c     ecp: get NUC and NRC
c-----------------------------------------------------------------------
c     i : counter of atoms
c     j : counter of l-quantumnumbers of the atom
c     k : counter of jco
C
C    j ---> ....
C
            if (iecp.ne.0) then
              if (k.le.jcostr(i,j)) then
                nucstr(i,j,k)=NUC(J,K)
                nrcstr(i,j,k)=NRC(J,K)
              endif
            endif
cch---------------------------------------------------------------------
      IMIN=IMAX+1
      IMAX=IMAX+NUC(J,K)
C-----------------------------------------------------------------------
C     Loop over primitives in each shell, reading exponents and cont-
C     raction coefficients.
C-----------------------------------------------------------------------
      DO 88140 LL = 1, NUC(J,K)
      IF (IGENER .EQ. 0 .AND. IGENES .EQ. 0) THEN
        READ(LUSCR) A(IMIN-1+LL), (CTEMP(LL,MC), MC = 1, NRC(J,K))
      ELSE
        IF (KASET .NE. '      ') THEN
          IF (K .LE. JCO(J)) THEN
            IF (IF10 .EQ. 1) THEN
CSSS              READ(LUBAS,211) A(IMIN-1+LL), (CTEMP(LL,MC),MC=1,NRC(J,K))
              READ(LUBAS,*) A(IMIN-1+LL), (CTEMP(LL,MC),MC=1,NRC(J,K))
            ELSE
CSSS              READ(LUBAS,14) A(IMIN-1+LL), (CTEMP(LL,MC),MC=1,NRC(J,K))
              READ(LUBAS,*) A(IMIN-1+LL), (CTEMP(LL,MC),MC=1,NRC(J,K))
            ENDIF
C
C....    SCALE S FUNCTIONS IF DESIRED
C
            IF (J .EQ. 1 .AND. SCALE .NE. 0.) THEN
              DO 88960 L = IMIN, IMAX
                A(L) = A(L)*SCALE*SCALE
88960         CONTINUE
            ENDIF
          ELSE
            IF (IF10 .EQ. 1) THEN
CSSS              READ(INU,211) A(IMIN-1+LL), (CTEMP(LL,MC),MC=1,NRC(J,K))
              READ(INU,*) A(IMIN-1+LL), (CTEMP(LL,MC),MC=1,NRC(J,K))
            ELSE
CSSS              READ(INU,14) A(IMIN-1+LL), (CTEMP(LL,MC),MC=1,NRC(J,K))
              READ(INU,*) A(IMIN-1+LL), (CTEMP(LL,MC),MC=1,NRC(J,K))
            ENDIF
          ENDIF
        ELSE
          IF(IF18.EQ.1)GO TO 88212
          IF(IF10.EQ.1)GO TO 88210
CSSS          READ(INU,14) A(IMIN-1+LL), (CTEMP(LL,MC),MC = 1,NRC(J,K))
          READ(INU,*) A(IMIN-1+LL), (CTEMP(LL,MC),MC = 1,NRC(J,K))
          GO TO 88215
88210     CONTINUE
CSSS          READ(INU,211) A(IMIN-1+LL), (CTEMP(LL,MC),MC = 1,NRC(J,K))
          READ(INU,*) A(IMIN-1+LL), (CTEMP(LL,MC),MC = 1,NRC(J,K))
          GO TO 88215
88212     CONTINUE
CSSS          READ(INU,213) A(IMIN-1+LL), (CTEMP(LL,MC),MC = 1,NRC(J,K))
          READ(INU,*) A(IMIN-1+LL), (CTEMP(LL,MC),MC = 1,NRC(J,K))
C
C This is added for IP/EA work of Monika (using a very small exponent 
C one can use a exciattion energy code to do IP/EA). 03/2010, 
C Ajith Perera
C
          If (A(IMIN-1+LL) .EQ. 0.0D0) A(IMIN-1+LL) = 1.0D-57

88215     CONTINUE
        ENDIF
      ENDIF
C-----------------------------------------------------------------------
C     End loop over primitives.
C-----------------------------------------------------------------------
88140 CONTINUE
C
C....    NOW RENORMALIZE COEFFICIENT ARRAY AND MOVE FROM CTEMP TO C.
C....    INPUT COEFFICIENTS (POSSIBLY UNNORMALIZED) ARE HELD IN CUN.
C
      IF (IBORTH .EQ. 0) THEN
C-----------------------------------------------------------------------
C     Loop over contracted functions in each shell.
C-----------------------------------------------------------------------
        DO 89150 MC = 1, NRC(J,K)
        IIMIN = IIMAX
        SUM = 0.
        DO LL = 1, NUC(J,K)
           AL = A(IMIN-1+LL)
           DO LM = 1, LL
              AM = A(IMIN-1+LM)
              dtmp = (2.0*DSQRT(AL*AM)/(AL+AM))
              T =     sqrt(dtmp)
              T = T * ( dtmp**J )
              T = T * CTEMP(LL,MC) * CTEMP(LM,MC)
              SUM = SUM + T
              IF (LL .NE. LM) SUM = SUM + T
           END DO
        END DO
        XNORM = 1./DSQRT(SUM)
        IF(DABS(XNORM-1.D0).GT.1.D-3)THEN
         WRITE(6,22222)I,ISPD(J),MC,XNORM
22222    FORMAT(T3,'@READIN-W, Basis function not normalized: ',/,
     &          T3,'           Center (as ordered in MOL file): ',I3,/,
     &          T3,'           Angular Momentum type          : ',A,/,
     &          T3,'           Offset within angular momentum : ',I2,/,
     &          T3,'           <chi|chi>                      : ',
     &          F15.10)
        ENDIF
        DO 89180 LL = 1, NUC(J,K)
        AL = A(IMIN-1+LL)
        CUN(IIMIN+LL) = CTEMP(LL,MC)
        C(IIMIN+LL) = CTEMP(LL,MC)*PIPPI*XNORM*(4.*AL)**(0.5*J+0.25)
89180   CONTINUE
        IIMAX = IIMAX + NUC(J,K)
C-----------------------------------------------------------------------
C     End of loop over contracted functions.
C-----------------------------------------------------------------------
89150   CONTINUE
      ELSE
C-----------------------------------------------------------------------
C     Loop over contracted functions in each shell.
C-----------------------------------------------------------------------
      DO 88150 MB = 1, NRC(J,K)
      IIMIN = IIMAX
      DO 88150 MC = 1, MB
      SUM = 0.
      DO 88170 LL = 1, NUC(J,K)
      AL = A(IMIN-1+LL)
      DO 88170 LM = 1, NUC(J,K)
      AM = A(IMIN-1+LM)
      T = CTEMP(LL,MB)*CTEMP(LM,MC)*(2.0*DSQRT(AL*AM)/(AL+AM))**(J+0.5)
88170 SUM = SUM+T
      IF (MC .eq. MB) THEN
      do 88201 ll=nuc(j,k),1,-1
      if(abs(ctemp(ll,mc)).gt.0.1) goto 88202
88201 continue
      print *,' Normalization failed. largest C<0.1'
      stop
88202 ac=ctemp(ll,mc)/abs(ctemp(ll,mc))
      DO 88180 LL = 1, NUC(J,K)
      AL = A(IMIN-1+LL)
      ctemp(ll,mc) = ac*ctemp(ll,mc)/sqrt(SUM)
      CUN(IIMIN+LL) = CTEMP(LL,MC)
      C(IIMIN+LL)=CTEMP(LL,MC)*PIPPI*(4.*AL)**(0.5*J+0.25)
88180 CONTINUE
      IIMAX = IIMAX + NUC(J,K)
      ELSE
      do 88190 LL=1, NUC(J,K)
      CTEMP(LL,MB) = CTEMP(LL,MB) - SUM*CTEMP(LL,MC)
88190 CONTINUE
      ENDIF
88150 CONTINUE
C-----------------------------------------------------------------------
C     End of loop over contracted functions.
C-----------------------------------------------------------------------
      ENDIF
88130 CONTINUE
C-----------------------------------------------------------------------
C     End of loop over shells of each angular momentum.
C-----------------------------------------------------------------------
      IF (IGENER .NE. 0 .OR. IGENES .NE. 0) JCO(J) = JCO(J) + IAD(J)
2     CONTINUE
C-----------------------------------------------------------------------
C     End of loop over angular momentum for each nuclear type.
C-----------------------------------------------------------------------
C
 1272 CONTINUE
C
C-----------------------------------------------------------------------
C     Loop over centers of current type, fill arrays with exponents,
C     contraction coefficients, symmetry information, and data for each
C     shell.
C-----------------------------------------------------------------------
      DO 5 N=NMIN,NMAX
C
C-----------------------------------------------------------------------
C     Load up exponents of primitives.
C-----------------------------------------------------------------------
        IF(LMAX+IMAX .GT. MXTNPR)THEN
         WRITE(6,*) ' @READIN-F, Too many primitives. Limit is ',MXTNPR
         CALL ERREX
        ENDIF
C
        DO 7 L=1,IMAX
    7   ALPHA(L+LMAX) = A(L)
        LMAX = LMAX + IMAX
C
C-----------------------------------------------------------------------
C     Load up contraction coefficients.
C-----------------------------------------------------------------------
        IF(LLMAX+IIMAX .GT. MXTNCC)THEN
         WRITE(6,*)
     &   ' @READIN-F, Too many contraction coefficients. Limit is ',
     &                                                       MXTNCC
         CALL ERREX
        ENDIF
C
        DO 1273 L=1,IIMAX
          CONTUN(L+LLMAX) = CUN(L)
 1273   CONT(L+LLMAX) = C(L)
        LLMAX = LLMAX+ IIMAX
        MULK=1
        LL=1
        DO 9 L=1,NSYMOP
          DO 8 M=1,3
            IF(AND(LL+1,ISYTYP(M)).NE.1 .AND. ABS(CORD(N,M)).GE.1.E-6)
     &      GO TO 9
    8     CONTINUE
          MULK=MULK+LL
    9   LL=2*LL
        MAXLO=AND(MULK,MAXLO)
        MULNUC(N)=MULK
        CHARGE(N)=Q
        NHTYP=MAX0(NHTYP,IQM)
        IF(IQM.EQ.0) GO TO 5
        DO 1274  J=1,IQM
          JCOJ=JCO(J)
          DO 1274 K=1,JCOJ
            KMAX=KMAX+1
C
            IF(KMAX .GT. MXTNSH)THEN
             WRITE(6,*) ' @READIN-F, Too many shells. Limit is ',MXTNSH
             CALL ERREX
            ENDIF
C
            NUCO(KMAX)=NUC(J,K)
            NRCO(KMAX) = NRC(J,K)
            LAMN(KMAX)=NAMN(N)
            MUL(KMAX)=MULK
            NHKT(KMAX)=J
            DO 1274 M=1,3
 1274   CENT(M,KMAX)=CORD(N,M)
    5 CONTINUE
C-----------------------------------------------------------------------
C     End of loop over nuclear types. All of MOL file has now been read.
C-----------------------------------------------------------------------
C
      CALL VPTAG
      CALL RFINIT
      if(idelb .gt. 0) read 828, (kdelb(i), i=1,idelb)
      NPARSU(1)=0
      MST(1) = 0
      NORB=1
      norxb=0
C
C     Tracks the number of deleted orbitals and where
C     their symmetry transformations are stored
C
      NOrbDel = MxCBF
C
      IBBAS=0
      DO 19 LA=1,MAXLOP
         IF(LA.EQ.1) IABAS=0
         KABAS=0
         NBAS(LA)=0
         JSTRT(1)=0
         JRS(1) = 0
         NSTRT(1)=0
         DO 20 IA=1,KMAX
            NHKTA=NHKT(IA)
            KHKTA=KHKT(NHKTA)
            IF(IA.NE.KMAX) THEN
               JSTRT(IA+1)=JSTRT(IA)+NUCO(IA)
               JRS(IA+1) = JRS(IA) + NRCO(IA) * NUCO(IA)
               NSTRT(IA+1)=NSTRT(IA)+KHKTA
            ENDIF
            DO 200 NA=1,KHKTA
               MXYZ=((NHKT(IA)+1)*NHKT(IA)*(NHKT(IA)-1))/6+NA
               IFIF = AND(MUL(IA),EOR(LA,ITYPE(NHKT(IA),NA))) .EQ. 1
               DO 199 JAA = 1, NRCO(IA)
                  JKB=0
                  DO 1515 KB=1,MAXLOP
                     IF(AND(KB,MUL(IA)) . NE . 1) GO TO 1515
                     JKB=JKB+1
                     KABAS=KABAS+1
C
                     IF(NORB.GT.MXCBF)THEN
                      WRITE(6,*) ' @READIN-F, NORB too large ',NORB
                      CALL ERREX
                     ELSEIF(JKB.GT.224)THEN
                      WRITE(6,*) ' @READIN-F, JKB too large ',JKB
                      CALL ERREX
                     ENDIF
C
                     CTRAN(NORB,JKB)=
     $                  PC(AND(KB,EOR(LA,ITYPE(NHKT(IA),NA))))
                     ITRAN(NORB,JKB)=KABAS
                     JTranC(NOrb) = JKB
                     CTranC(NOrb, JKB) = CTran(NOrb, JKB)
                     ITranC(NOrb, JKB) = ITran(NOrb, JKB)
 1515             CONTINUE
                  IF(IFIF) THEN
                     dele = (idela.eq.0) .and. (jcntam(na,nhkta).eq.1)
                     if(idelb .gt. 0) then
                        do 777 i=1,idelb
 777                    dele = dele .or. kdelb(i).eq.norb
                     endif
                     if (dele) then
                        newind(norb)= -NOrbDel
C
                        IF(NORBDEL.GT.MXCBF)THEN
                         WRITE(6,*) ' @READIN-F, NORBDEL too large ',
     &                                                        NORBDEL
                         CALL ERREX
                        ENDIF
C
                        JTRAN(NOrbDel)=JKB
                        JPRX(NOrbDel)=KWO(MXYZ)
                        MAMN(NOrbDel)=LAMN(IA)
C     
C     NOrbDel starts at a useful value, so the update is done _after_ the
C     value is used.
C
                         NOrbDel = NOrbDel - 1
                      else
C
                        IF(NORXB.GT.MXCBF)THEN
                         WRITE(6,*) ' @READIN-F, NORXB too large ',
     &                                                       NORXB
                         CALL ERREX
                        ENDIF
C
                         norxb=norxb+1
                         newind(norb)=norxb
                         JTRAN(NORXB)=JKB
                         JPRX(NORXB)=KWO(MXYZ)
                         MAMN(NORXB)=LAMN(IA)
                      endif
                      NORB=NORB+1
                      NBAS(LA)=NBAS(LA)+1
                    ENDIF
 199             CONTINUE
 200          CONTINUE
  20      CONTINUE
          IF(LA.GT.1) THEN
             NPARSU(LA)=NPARSU(LA-1)+NPAR(LA-1)
             MST(LA) = MST(LA-1) +  NBAS(LA-1)
          ENDIF
 19    NPAR(LA)=NBAS(LA)*(NBAS(LA)+1)/2
c
      if (lharm) then
c     map ia,nb,jaa and kb onto 1-dim. seqno kabs(*,*,*,*)
      kab=0
      do 1001 ia=1,kmax
      khkta=khkt(nhkt(ia))
      khkta=nhkt(ia)*(nhkt(ia)+1)/2
      do 1001 nb=1,khkta
      do 1001 jaa=1,nrco(ia)
      do 1001 kb=1,maxlop
      if (and(kb,mul(ia))  .eq.  1) then
      kab=kab+1
      kabs(ia,nb,jaa,kb)=kab
      endif
1001  continue
c
c     map 1-dim. seqno onto la,ia,na and jaa
      noryb=0
      do 1002 la=1,maxlop
      do 1002 ia=1,kmax
      khkta=nhkt(ia)*(nhkt(ia)+1)/2
      do 1002 na=1,khkta
      do 1002 jaa=1,nrco(ia)
      if(and(mul(ia),eor(la,itype(nhkt(ia),na))) .ne. 1) go to 1002
      noryb=noryb+1
C
      IF(NORYB.GT.MXCBF)THEN
       WRITE(6,*) ' @READIN-F, NORYB too large ',NORYB
       CALL ERREX
      ENDIF
C
      last(noryb)=la
      iast(noryb)=ia
      nast(noryb)=na
      jaast(noryb)=jaa
1002  continue
c     
c     find connection between SO's and CGTO's
C     
C     NOTE: NOrb is one too _large_ from the way it is updated.  Hence
C     the limit norb-1.  The value of NOrb is corrected later on.
C
      WRITE(FILEOUT2,'(A)') '* start'
      do 1003 i=1,norb-1
      inew=newind(i)
C
      IF(INEW.GT.MXCBF)THEN
       WRITE(6,*) ' @READIN-F, INEW too large ',INEW
       CALL ERREX
      ENDIF
C
C     Instead of skiping deleted functions, construct their
C     transformations.  Remember they live at the very end of
C     the JTran/CTran/ITran arrays, so they won't harm anything.
C
      if(inew .le. 0) inew = -inew
      jkb=0
      la=last(i)
      ia=iast(i)
      khkta=khkt(nhkt(ia))
      na=nast(i)
      jaa=jaast(i)
      iaq=iadq(nhkt(ia))+(na-1)*khkta
      inx=(nhkt(ia)+1)*(nhkt(ia)-1)*nhkt(ia)/6 + na
      do 1004 kb=1,maxlop
      if(and(kb,mul(ia)) .ne. 1) go to 1004
      p=pc(and(kb,eor(la,itype(nhkt(ia),na))))
      do 1005 nb=1,khkta
      ch=charm(iaq+nb)
      if (abs(ch) .lt. 1.e-03) go to 1005
      kab=kabs(ia,nb,jaa,kb)
      jkb=jkb+1
      ctran(inew,jkb)=p*ch
      WRITE(FILEOUT2,'(7I5,2F5.1)') i,kb,khkta,nb,inew,jkb,kab,p,ch
      itran(inew,jkb)=kab
 1005 continue
      WRITE(FILEOUT2,'(5I5,F5.1,I5,50F5.1)') i,kb,inew,jkb,kab,
     &     p,khkta,(charm(iaq+nb),nb=1,khkta)
 1004 continue
      jtran(inew)=jkb
      jjprx(inew) = kktyp(inx,1+idosph)
 1003 continue
      WRITE(FILEOUT2,'(A)') '* end'
      endif
C NEWIND: spherical index
      WRITE(LUTR1,'(A)')'Trans1'
      DO I=1,NORB-1
          WRITE(LUTR1,'(2I4)') I,NEWIND(I)
      ENDDO
#ifdef _OLDECP
cch---------------------------------------------------------------------
c     ecp-part
c-----------------------------------------------------------------------
      if (iecp.ne.0) then
c-----------------------------------------------------------------------
c     grep atomsymbols and put them to mtype-array
c-----------------------------------------------------------------------
        ierr=0
        call chrges(CHARGE,mtype,nmax,ierr)
c-----------------------------------------------------------------------
c     read ecp data and prepare data for evulation of the
c     integrals
c-----------------------------------------------------------------------
        if (ierr.eq.0) then
          ipr=0
          iunit=55
          infile='ECPDATA'
c         call ecpdef(ipr,nmax,iatdef,ipseux,ipseud,mtype,cnick,
c    &                iunit,infile,ierr)
          call ecpdef(nmax,ipseux,mtype,cnick)
          if (ierr.eq.0) then
            ipr=0
            natoms=NMAX
            call ecppar(ipr,ipseux,cnick,CHARGE,natoms,.FALSE.,0)
c-----------------------------------------------------------------------
c           Summation over the charges and storage of this value
c-----------------------------------------------------------------------
            chgsum=0
            do 780 i=1,natoms
C
C account for symmetry !
C
              chgsum=chgsum+charge(i)*fmult(mulnuc(i))
  780       continue
            nproton=idint(chgsum)
            call putrec(20,'JOBARC','NMPROTON',1,nproton)
c-----------------------------------------------------------------------
c           ecp: prepare data for ecpintegralevaluation
c-----------------------------------------------------------------------
            call doit(namat,iqmstr,jcostr,nucstr,nrcstr,NHARM)
c-----------------------------------------------------------------------
          else
            write(*,*) 'E R R O R in subroutine ecppar!'
            write(*,*) 'Please check your ecp-inputfile ECPDATA.'
            stop
          endif
        else
          write(*,*) 'E R R O R in subroutine chrges!'
          write(*,*)
          stop
        endif
c-----------------------------------------------------------------------
c     omega
c-----------------------------------------------------------------------
      endif
c-----------------------------------------------------------------------
C----- New ecplibray
#else
       If (iecp.ne.0) then
          Call Getrec(20, 'JOBARC', 'NREALATM', 1, NATOMS)
          Call Ecp_main(Iecp,Nmax,namat,iqmstr,jcostr,nucstr,nrcstr)
       Endif 
#endif
C
C
CJDW 6/6/95. Add Ajith's block of PUTRECs.
C
C  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  &                                                               &
C  & Added by ajith 03/94. Information required to symmetry adapt  &
C  & nuclear perturbations.                                        &
C  &                                                               &
C  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C      
      CALL PUTREC (20, 'JOBARC', 'MULCINP ', 100, MULNUC)
      CALL PUTREC (20, 'JOBARC', 'SYMTYPE ', 3 , ISYTYP)
      CALL PUTREC (20, 'JOBARC', 'NOSYMMOP', 1 , NSYMOP)
      CALL PUTREC (20, 'JOBARC', 'NOSYATOM', 1 , NMAX) 

#ifdef _NAMN_THROUGH_CHRTMP
c   o NAMN is divided into 10 'columns'
      DO J = 0, 9
c      o CHRTMP will hold column 1+J
         WRITE (CHRNUM,'(I1)') J
         DO I = 0, DIM_CHRTMP-4, 4
            CHRTMP(1+I:4+I) = NAMN(1+(DIM_CHRTMP*J+I)/4)
         END DO
         CALL PUTCREC (20,'JOBARC','ATMSYMB'//CHRNUM,DIM_CHRTMP,CHRTMP)
      END DO
#else
      CALL PUTREC (20, 'JOBARC', 'ATOMSYMB', 300, NAMN)
#endif /* _NAMN_THROUGH_CHRTMP */
C
CEND OF ADDITION
C This block of code should get executed only when IRC searches
C are used. For those cases, the geometries are read from the 
C JOBARC file and the MOL file should correspond to the IRC reference
C point (usually the TS), 04/2019. Ajith Perera.

      Call GETREC(0,"JOBARC","COORD   ",Length,CORR)
      Ipes_scan =  Iflags2(h_IFLAGS2_pes_scan)
      If (Length .GT. 0 .AND. Ipes_scan .EQ. 1) Then 
         Call GETREC(20,"JOBARC","COORD   ",Length,CORR)
         Icoff = 0
         Do Icent = 1, Nmax
            Do Ixyz = 1,3
               Cord(Icent,Ixyz) = Corr(icoff+Ixyz)
            Enddo 
            Icoff=Icoff+ 3
         Enddo 
      Endif 

      WRITE(FILEOUT,'(A)') '* geometry'
      NNUC = 0
      REWIND LUSCR
      POTNUC=0.0D0
      if(iprint.ge.1) write(6,126)
126   format('   Internuclear distances (A) :')
      DO 25 N=1,NMAX
      if(iprint.ge.1) write(6,127) namn(n),(cord(n,i),I=1,3)
127    format(/,5x,'for atom ',a4,' (coordinates :',3F10.5,')',/)
      NDEG(N) = 1
      NNUC = NNUC + 1
      IONE=1
      WRITE(LUSCR) NAMN(N),IONE, (CORD(N,I),I = 1,3), CHARGE(N)
      NAME = NAMN(N)
      WRITE(FILEOUT,'(A,I3,X,A2,X,A2,X,I3,X,3F20.10)') 'o',NNUC,
     &      NAME(1:2),NAME(3:4),1,(CORD(N,I),I=1,3)
      DO 25 M=N,NMAX
      DO 26 KB=1,MAXLOP
      IF (IFBA(M).EQ.0 .and. IFBA(N).EQ.0) GOTO 26
      IF(AND(KB,OR(MULNUC(N),MULNUC(M))).NE.1) GO TO 26
      IF(M.EQ.N.AND.KB.EQ.1) GO TO 26
      HKAB = FMULT(AND(MULNUC(M),MULNUC(N)))
      DIST=0.0D0
      IF(M.EQ.N) HKAB=0.5D0*HKAB
      DO 27 I=1,3
      X(I) = CORD(M,I)*PTAND(ISYTYP(I),KB)
   27 DIST=DIST+(CORD(N,I)-CORD(M,I)*PTAND(ISYTYP(I),KB))**2
      IF(DIST .LT. 0.2 .AND. CHARGE(M)*CHARGE(N).NE.0.0) THEN
      WRITE (6,29) N,KB,M,DIST
   29 FORMAT (12H THE NUCLEUS ,I3,20H IS TOO CLOSE TO THE ,I3,
     1 32H*TH TRANSFORMATION OF NUCLEUS NO,I3, 6H DIST=,E15.5)
      CALL ERREX
      ENDIF
      RO=DSQRT(DIST)
      AVS=0.5291771*RO
      if(iprint.ge.1)write(6,128) namn(m),kb,avs
128   format(5x,a4,i4,5x,'distance is ',f6.3)
      IF (M .EQ. N) THEN
        NNUC = NNUC + 1
        WRITE(LUSCR) NAMN(N), KB, (X(I),I = 1,3), CHARGE(N)
        WRITE(FILEOUT,'(A,I3,X,A4,X,I3,X,3F20.10)') 'c',NNUC,NAMN(N),KB,
     &     (X(I),I=1,3)
      ENDIF
      POTNUC=POTNUC+CHARGE(M)*CHARGE(N)*HKAB/RO
   26 CONTINUE
   25 CONTINUE
      WRITE(FILEOUT,'(A)') '* AO list'
      NSYM=0
      IF (LPRT) WRITE (LUT,1112)
1112  format(/,'   Group multiplication table :',/)
      DO 1113 I=1,MAXLOP
1113  IF (LPRT) WRITE (LUT,1114) (EOR(I,J),J=1,MAXLOP)
1114  FORMAT (6x,8I4)
      IF (LPRT) WRITE (LUT,1121)
1121  format(/,'   Gaussian basis information :',//,
     &       3x,'atom',17x,'exponent',6x,'coefficients')
      IAB=0
      KNTR=0
      JBBAS = 0
      DO 1140 IA=1,KMAX
      NHKTA=NHKT(IA)
      NAK=(NHKTA*(NHKTA+1)*(NHKTA-1))/6
      KHKTA=KHKT(NHKTA)
      MULA=MUL(IA)
      LAM=LAMN(IA)
      JRSIA = JRS(IA)
      JSTRIA = JSTRT(IA)
      NSTRIA = NSTRT(IA)
      DO 1141 NA=1,KHKTA
      NUCA=NUCO(IA)
      NRCA = NRCO(IA)
      NAKNA=NAK+NA
      NAKO=KWO(NAKNA)
C
C....    CONSTRUCT PSEUDO-SEGMENTED BASIS FOR PROPERTIES AND ECP
C....    PROGRAMS BY DUPLICATING PRIMITIVES IN CONTRACTED FUNCTIONS.
C
      DO 81144 JAA = 1,NRCA
      DO 81142 KB=1,MAXLOP
      IF(AND(KB,MULA).NE.1) GO TO 81142
      KBBAS = 0
      DO 81154 JA = 1,NUCA
      EX = ALPHA(JSTRIA+JA)
      CON = CONTUN(JRSIA+JA+(JAA-1)*NUCA)
      IF (ABS(CON) .GT. 1.E-8) THEN
      JBBAS = JBBAS + 1
      KBBAS = KBBAS + 1
      WRITE(LUSCR) LAM, KB, NAKO, EX, CON
      ENDIF
81154 CONTINUE
      KNTR = KNTR + 1
      ICDS(KNTR) = KBBAS
81142 CONTINUE
81144 CONTINUE
      DO JAA = 1, NRCA
      DO KB=1,MAXLOP
         IF(AND(KB,MULA).EQ.1) THEN
         WRITE(FILEOUT,'(1X,A2,X,A2,X,3I3,X,A5,3I5)') LAM(1:2),
     &        LAM(3:4),KB,
     &        KHKTA,NA,NAKO,NUCA,NRCA,JAA
         ENDIF
      ENDDO
      ENDDO
      DO 1142 KB=1,MAXLOP
      IF(AND(KB,MULA).NE.1) GO TO 1142
      IF (LPRT) WRITE (LUT,1143) LAM,KB,NAKO
1143  FORMAT (/1X,A4,I3,3X,A5)
      DO 1144 JA=1,NUCA
      IAB=IAB+1
      IOVERP = ' '
      IF (JA .EQ. 1) IOVERP = '+'
      IF (LPRT) WRITE (LUT,1145) IOVERP, IAB, ALPHA(JSTRIA+JA),
     1                     (CONTUN(JRSIA+JA+(JAA-1)*NUCA),JAA = 1,NRCA)
1145  FORMAT(A1,I18,F15.6,8F12.7,/,38X,8F12.9)
1144  CONTINUE
1142  CONTINUE
1141  CONTINUE
1140  CONTINUE

      CALL PUTREC (20, 'JOBARC', 'NBASCART', 1, (norb-1))

      if(iprint.ge.1)then
       write (lut,828) norb-1,norxb
       write (lut,828) (newind(i), i=1,norb-1)
828    format (25i3)
      endif
      ji=0
      norxb  = 0
      do 486 i=1,maxlop
      nxbas(i) = 0
        do 485 j=1,nbas(i)
        ji=ji+1
        if (newind(ji) .le. 0) go to 485
        nxbas(i) = nxbas(i) + 1
        norxb = norxb+1
        LABORB(NORXB) = I*1024 + NXBAS(I)
  485  continue
  486  continue
      IF (LPRT) WRITE (LUT,1115)
1115  FORMAT (1H1,//12X,29H SYMMETRY TRANSFORMATION INFO )
      NORB=NORB-1
      IF (LPRT) WRITE (LUT,1122) NORXB
1122  FORMAT (10X,I5,34H SYMMETRY ADAPTED BASIS FUNCTIONS  )
C
C Trap the cases where the No. of Symmetry adapted basis functions
C greater than 256 on 32 bit machines. Ajith August 1997.
C

      IF (IINTFP .EQ. 2 .AND. NORXB .GT. 256) THEN
         Write(6,*)
         Write(6,*)' The number of symmetry adapted basis functions =',
     &              NORXB
         Write(6,*)' More than 256 basis functions are not allowed'
         Write(6,*)' in this architecture.'
         Call ERREX
      ENDIF
C
      I=0
      NSYM=0
      DO 1126 LA=1,MAXLOP
      NBI=NXBAS(LA)
      IF(NBI.EQ.0) GO TO 1126
      IF (LPRT) WRITE (LUT,1124) LA
1124  FORMAT (//10X,33Hirreducible representation number ,I2/)
      NSYM=NSYM+1
      DO 1116 ILA=1,NBI
      I=I+1
      J=JTRAN(I)
      IF (LPRT) WRITE (LUT,1000)I,ILA,MAMN(I),JJPRX(I),J,
     1 (ITRAN(I,K),CTRAN(I,K),K=1,J)
1116  CONTINUE
1000  FORMAT (2I5,2X,A6,2X,A5,I4,8(I4,F6.1)/(29X,8(I4,F6.1)))
1126  CONTINUE
      do 487 j=1,maxlop
      npar(j) = nxbas(j)*(nxbas(j)+1)/2
      if(j .eq. 1) goto 487
      nparsu(j) = nparsu(j-1) + npar(j-1)
      mst(j) = mst(j-1) + nxbas(j-1)
      MSTOLD(J) = MSTOLD(J-1) + NBAS(J-1)
487   continue
4870  WRITE(NTAP) TITLE,BLANKF,MAXLOP,(NXBAS(I),I=1,MAXLOP),POTNUC,
     1            IDOSPH
      REWIND LUSCR
      WRITE(NTAP) POLINP(1)
      WRITE(NTAP) POLINP(2)
      WRITE(NTAP) NNUC, IZILCH, IZILCH, IZILCH
      DO 300 I = 1, NNUC
      READ(LUSCR) LAM, II, CX, CY, CZ, CH
      WRITE(NTAP) LAM, II, CX, CY, CZ, CH
300   CONTINUE
      WRITE(NTAP) JBBAS, KNTR, IZILCH, IZILCH
      WRITE(NTAP) (ICDS(I),I = 1,KNTR), IZILCH, IZILCH, IZILCH
      DO 301 I = 1, JBBAS
      READ(LUSCR) LAM, KB, NAKO, EX, CON
      WRITE(NTAP) LAM, KB, NAKO, EX, CON
301   CONTINUE
      WRITE(NTAP) SYMINP(1)
      WRITE(NTAP) SYMINP(2)
      WRITE(NTAP) NORXB, IZILCH, IZILCH, IZILCH
      DO 302 I = 1, NORXB
c     if (newind(i)  .le.  0) goto 302
      J = JTRAN(I)
#ifdef _DEBUG_LVLM
      Write(6,"(i3)") J
      Write(6,"(I3,1x,I3,1X,F5.1,1x,i3,1X,F5.1,1x,I3,F5.2)") 
     &          i,(iTRAn(I,k), CTRAN(I,K), K = 1,J)
#endif
      WRITE(NTAP) J, (ITRAN(I,K), CTRAN(I,K), K = 1,J), IZILCH
302   CONTINUE
C
C     After the actual transformation, write the rest of the
C     matrix (deleted orbitals) necessary to make it square.
C     Earlier we placed these at the very end of the arrays.
C     
C     NOTE: NOrbDel is one too _small_ because if the way it is updated.
C     Hence, the limit NOrbDel+1.
C
      DO 3021 I = MxCBF, NOrbDel+1, -1
      J = JTRAN(I)
      WRITE(NTAP) J, (ITRAN(I,K), CTRAN(I,K), K = 1,J), IZILCH
 3021 CONTINUE

C     This is the full cartesian basis AO -> SO transformation.
C     Include in this info the number of function in each irrep too,
C     because there is nowhere else to get it.
C
C     NOTE: Write this record ONLY if it is different from the
C     computational basis.
C
      If ( LHarm .OR. NOrb .ne. Norxb) then
         Write (NTap) CSymIn(1)
         Write (NTap) CSymIn(2)
         Write (NTap) NOrb, IZILCH, IZILCH, IZILCH,
     &      (NBas(i), i=1, MAXLOP)
C
         Do 3022 I = 1, NOrb
            J = JTranC(I)
            Write (NTap) J, (ITranC(I,K), CTranC(I,K), K = 1,J), IZILCH
 3022    Continue
      EndIf

      Write (LUTR1,*) "test1"
      DO I = 1, NOrb
         J = JTranC(I)
         Write (LUTR1,'(I4,20(I5,F11.7))')
     &         J, (ITranC(I,K), CTranC(I,K), K = 1,J)
      ENDDO
C
CJFS LABELS - JUNK WRITTEN SO THAT SEARCH ROUTINE WILL WORK AND
C  NOT READ OFF END OF RECORD
C
      WRITE(NTAP)LABINP(1)
      WRITE(NTAP)LABINP(2)
      DO 3333 I=1,NORXB
       WRITE(NTAP)I,MAMN(I),JPRX(I),JUNK
3333   CONTINUE
      WRITE(NTAP) POPINP(1)
      WRITE(NTAP) POPINP(2)
      WRITE(NTAP) NMAX, NORXB, IZILCH, IZILCH
      DO 303 I = 1,NMAX
      CHARGT(I) = CHARGE(I)*NDEG(I)
303   CONTINUE
      WRITE(NTAP) NMAX, (NAMN(I), CHARGT(I), I = 1, NMAX), IZILCH
      DO 304 J = 1, NORXB
c     if (newind(j) .le. 0) goto 304
      DO 305 I = 1, NMAX
      IF (NAMN(I) .EQ. MAMN(J)) GOTO 306
305   CONTINUE
306   CONTINUE
      DO 307 K = 1, NHT*(NHT+1)*(NHT+2)/6
      IF (KWO(K) .EQ. JPRX(J)) GOTO 308
307   CONTINUE
308   CONTINUE
      WRITE(NTAP) J, I, K, IZILCH
304   CONTINUE
      REWIND LUSCR
C
c      WRITE (6,97) POTNUC
c   97 FORMAT (////10X,'NUCLEAR REPULSION ENERGY=',E13.7////)
      write(6,97) potnuc
97    format(T3,'@READIN-I, Nuclear repulsion energy : ',F15.10,' a.u.')
C...  REDUCE THE RANGE OF SYMMETRY LOOPS IF POSSIBLE
      IMX = 1.4427 * (0.47624 + ALOG(1.+MAXLOP-MAXLO))
      MAXLOT = 2**IMX
      CLOSE(UNIT=INU)
C     
C     Now we must restore NEWIND to what everything else expects:
C     deleted functions are represented by a zero instead of
C     a negative number.
C     
      Do 1999 i = 1, NOrb
         NewInd(i) = Max( NewInd(i), 0)
 1999 Continue
      WRITE(FILEOUT,'(A)') '* end'
      CLOSE(FILEOUT)
      CLOSE(FILEOUT2)
      CLOSE(LUTR1)
c
cOLD  call figmem(kmax,nhkt,nuco,nrco)
      call figmem(kmax)
#ifdef _TRACE_VMOL
      write(*,*) '@VMOL: Leaving file ',__FILE__
#endif
      RETURN
      END
