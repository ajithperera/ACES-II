         SUBROUTINE  DKH__INIT
     +
     +                   (IORDER,NBAS,
     +                    KIN,NAI,PVP,OVL,
     +                    ZCORE,
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

      INTEGER    IORDER,IPRINT,IERR
      INTEGER    ITMP1,ITMP2,IEND
      INTEGER    NBAS,ISIZE,IOFFSET,NDIM
      INTEGER    INT,NUT,IRR
      INTEGER    ILNBUF,ILENGTH,ICENT
      INTEGER    LUINT,LUPVP,IRWND 
      INTEGER    ZMAX, ND, NCENTERS 

      DOUBLE PRECISION  COREHAM (1:NBAS*NBAS)

      DOUBLE PRECISION  KIN (1:NBAS,1:NBAS)
      DOUBLE PRECISION  NAI (1:NBAS,1:NBAS)
      DOUBLE PRECISION  PVP (1:NBAS,1:NBAS)
      DOUBLE PRECISION  OVL (1:NBAS,1:NBAS)

      DOUBLE PRECISION  ZCORE (1:ZMAX)

      CHARACTER*80   FNAME

      LOGICAL BOPENED
C
C
C------------------------------------------------------------------------
C
C
C             ...Set the print variable
C
C
C             ...Grab the kinetic energy, electron-nuclear
C                attraction, and <p.Vp> integrals
C
C

       NDIM   = NBAS

#ifdef _DEBUG_LVLM
       IPRINT = 20
#else
       IPRINT = 0
#endif 
       CALL  DKH__MAIN_TRANS
     +
     +           ( IORDER,NDIM,
     +             IPRINT,ZMAX,
     +             ZCORE,
     +             OVL,NAI,PVP,IERR,KIN) 
     +
       CALL DCOPY 
     +
     +      ( NDIM*NDIM, ZCORE (1), 1, COREHAM, 1)
     +
     +
#ifdef _DEBUG_LVLM
      Write(6,*) "The DKH transformed one-electron Hamiltonian."
      call output(coreham,1,ndim,1,ndim,ndim,ndim,1)
#endif
C
C             ...ready!
C
C
         RETURN
         END
