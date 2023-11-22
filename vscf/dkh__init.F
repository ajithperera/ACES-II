#include "flags.h"
         SUBROUTINE  DKH__INIT
     +
     +                   (IORDER,NBAS,ISIZE,
     +                    ILNBUF,IBUF,BUF,
     +                    KE,V,SFULL,DVD,
     +                    KIN,NAI,PVP,OVL,
     +                    SCRTRI,ZCORE,
     +                    ZMAX,
     +
     +                                 COREHAM )
     +
     +             
C------------------------------------------------------------------------
C  OPERATION   : DKH__INIT
C  MODULE      : DOUGLAS-KROLL-HESS TRANSFORMATION
C  MODULE-ID   : DKH
C  SUBROUTINES : DKH__MAIN_TRANS
C                                    
C  DESCRIPTION : Interface that connects VSCF to the DKH library
C
C                  Input:
C
C                    IORDER       =  order of the transformation (1-5)
C                    NBAS         =  number of atomic orbital basis
C                                    functions
C                    ISIZE        =  NBAS * (NBAS+1) / 2
C                    ILNBUF       =  buffer length space for integrals
C                    IBUF         =  Integer array containing integral
C                                    indices
C                    BUF          =  array used to read the integrals
C                    KE           =  symmetry packed kinetic energy integrals
C                    V            =  symmetry packed nuclear attraction
C                                    integrals
C                    DVD          =  symmetry packed pVp integrals
C                    SFULL        =  symmetry packed overlap integrals
C                    OVL          =  scratch space for expanded overlap
C                                    integrals for a given irrep
C                    PVP          =  scratch space for expanded pVp
C                                    integrals for a given irrep
C                    KIN          =  scratch space for expanded kinetic
C                                    energy integrals for a given irrep
C                    NAI          =  scratch space for expanded nuclear
C                                    attraction integrals for a given irrep
C                    ZMAX         =  maximum floating point memory
C                    ZCORE        =  double precision scratch space
C
C                  Output:
C
C                    COREHAM      =  transformed one electron
C                                    hamiltonian matrix
C
C  
C  AUTHOR      : Thomas Watson Jr.
C------------------------------------------------------------------------
C
C
C             ...Declare variables and include files!
C
C
      IMPLICIT    NONE

#include "aces.h"
#include "iuhf.com"
#include "symm2.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"

      INTEGER    IORDER,IPRINT,IERR
      INTEGER    ITMP1,ITMP2,IEND
      INTEGER    NBAS,ISIZE,IOFFSET,NDIM
      INTEGER    INT,NUT,IRR
      INTEGER    ILNBUF,ILENGTH,ICENT
      INTEGER    LUINT,LUPVP,IRWND 
      INTEGER    ZMAX, ND, NCENTERS 

      PARAMETER  ( LUINT = 10 )
      PARAMETER  ( LUPVP = 30 )

      INTEGER           IBUF (1:ILNBUF)
      DOUBLE PRECISION  BUF  (1:ILNBUF)

      DOUBLE PRECISION  COREHAM (1:ITRILN(NIRREP+1))

      DOUBLE PRECISION  SFULL (1:ITRILN(NIRREP+1))
      DOUBLE PRECISION  V     (1:ITRILN(NIRREP+1))
      DOUBLE PRECISION  KE    (1:ITRILN(NIRREP+1))
      DOUBLE PRECISION  DVD   (1:ITRILN(NIRREP+1))

      DOUBLE PRECISION  SCRTRI (1:ISIZE)

      DOUBLE PRECISION  KIN (1:MXIRR2)
      DOUBLE PRECISION  NAI (1:MXIRR2)
      DOUBLE PRECISION  PVP (1:MXIRR2)
      DOUBLE PRECISION  OVL (1:MXIRR2) 

      DOUBLE PRECISION  ZCORE (1:ZMAX)

      CHARACTER*80   FNAME

      LOGICAL BOPENED
      LOGICAL CONTRACT
C
C
C------------------------------------------------------------------------
C
C
C             ...Set the print variable
C
C
      
        CONTRACT =  (IFLAGS2(h_IFLAGS2_uncontract) .eq.0)

C The contracted basis DKH one-particle term is computed using OED
C integrals in a separate module (dkh). It also works when the basis
C is uncontracted but I will leave the orginal path for those calcs. 
C Ajith Perera, 05/2017.

        IF (CONTRACT) THEN
           Ndim=(Nbas+1) * Nbas /2 
           Ndim=Itriln(Nirrep+1)
           Write(6,"(a,a)") "  Adding contracted DKH core-",
     +                      "Hamiltonian from OED/DKH module"
           CALL GETREC(20,"JOBARC","DKH_INTS",Ndim*IINTFP,COREHAM)
           RETURN
        ENDIF 

        IPRINT = IFLAGS (1)
C
C
C             ...Grab the kinetic energy, electron-nuclear
C                attraction, and <p.Vp> integrals
C
C
         CALL GFNAME ('IIII    ', FNAME, ILENGTH)
         INQUIRE     (FILE=FNAME(1:ILENGTH), OPENED=BOPENED)
         IF (.NOT. BOPENED) THEN
            OPEN(LUINT, FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     +           ACCESS='SEQUENTIAL')
         END IF

         CALL LOCATE (LUINT, 'OVERLAP ')
         CALL ZERO   (SFULL,ITRILN(NIRREP+1))
         NUT = ILNBUF
         DO WHILE (NUT .EQ. ILNBUF)
            READ(LUINT) BUF, IBUF, NUT
            DO INT = 1, NUT 
               SFULL(IBUF(INT)) = BUF(INT)
            END DO
         END DO

         CALL LOCATE(LUINT, 'KINETINT')
         CALL ZERO(KE,ITRILN(NIRREP+1))
         CALL ZERO(V,ITRILN(NIRREP+1))
         NUT = ILNBUF
         DO WHILE (NUT .EQ. ILNBUF)
            READ(LUINT) BUF, IBUF, NUT
            DO INT = 1, NUT
               KE(IBUF(INT)) = BUF(INT)
               V(IBUF(INT)) = COREHAM(IBUF(INT)) - BUF(INT)
            END DO
         END DO

         CALL GFNAME('VPOUT   ', FNAME, ILENGTH)
         INQUIRE(FILE=FNAME(1:ILENGTH), OPENED=BOPENED)
         IF(.NOT. BOPENED) THEN
            OPEN(LUPVP, FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     +           ACCESS='SEQUENTIAL')
         END IF

         CALL  ZERO(SCRTRI, ISIZE)
         IERR  = 1 
         IRWND = 0
C
C Obtain the number of symmetry unique centers and have a loop over 
C the number of centers. The PVP integral is constructed per atom 
C basis and they need to be summed up. The DKH was correct 
C only for atoms up until this fix is installed. 09/2016, Ajith Perera. 
C 
         CALL GETREC(20, "JOBARC", "NREALATM", 1, NCENTERS) 

         DO ICENT = 1, NCENTERS 

         CALL  SEEKLB('   PVP  ', IERR, IRWND)

         NUT = ILNBUF
         DO WHILE (NUT .EQ. ILNBUF)
            READ(LUPVP) BUF, IBUF, NUT
            DO INT = 1, NUT
               SCRTRI(ibuf(int)) = SCRTRI(ibuf(int)) + buf(int)
            END DO
         END DO

         IRWND = 1

         ENDDO 
         
         IF (NIRREP .GT. 1) THEN
            ITMP1 = 1
            ITMP2 = ITMP1 + NBAS * NBAS
            IEND  = ITMP2 + NBAS * NBAS
            CALL  EXPND2 (SCRTRI, ZCORE (ITMP1), NBAS)
            CALL  DKH__SYMPACK(ZCORE (ITMP1),
     +                         ZCORE (ITMP2),
     +                         NBAS, NBFIRR)

            CALL  DCOPY(ISIZE, ZCORE (ITMP2), 1, DVD, 1)
         ELSE
            CALL  DCOPY(ISIZE, SCRTRI, 1, DVD, 1)
         END IF
C
C
C             ...Call the main DKH Driver!
C
C
         DO IRR = 1, NIRREP

            IF (NBFIRR (IRR) .GT. 0) THEN

                NDIM    =  NBFIRR(IRR)
                ISIZE   =  ITRILN(IRR)
                IOFFSET =  ITRIOF(IRR)

                CALL  EXPND2 (KE    (IOFFSET), KIN, NDIM)
                CALL  EXPND2 (V     (IOFFSET), NAI, NDIM)
                CALL  EXPND2 (DVD   (IOFFSET), PVP, NDIM)
                CALL  EXPND2 (SFULL (IOFFSET), OVL, NDIM)

#ifdef _DEBUG_LVLM
                Write(6,*) " Check sum of PVP ints"
                call checksum ("@-DKH__init", PVP, NDIM*NDIM)
                call output(PVP,1,NDIM,1,NDIM,NDIM,NDIM,1)
#endif 
                
                CALL  DKH__MAIN_TRANS
     +
     +                    ( IORDER,NDIM,
     +                      IPRINT,ZMAX,
     +                      ZCORE,
     +                      OVL,PVP,
     +
     +                                IERR,
     +                                KIN,NAI )
     +
                CALL  SQUEZ2
     +
     +              ( ZCORE (1), COREHAM (IOFFSET), NDIM )
     +
     +
            ENDIF

         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
