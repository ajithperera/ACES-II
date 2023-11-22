C
      PROGRAM RSPDENS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POP,VRT,DIRPRD
CJDW KKB stuff
      INTEGER POP2,VRT2
CJDW END
      INTEGER POPRHF,POPDOC,VRTRHF,VRTDOC
      LOGICAL CIS,HFM3,EOM,LTRP
      LOGICAL GABCD,RELAXED,TRULY_NONHF
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      LOGICAL DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,TRIP2
C SG 11/20/97
      LOGICAL MRCC, SEWARD,SYMMETRY 
      LOGICAL NONHF_TERMS_EXIST, STANDARD_NONCAN
      logical bRedundant
C
      COMMON//ICORE(1)
      COMMON/EXCITE/CIS,EOM
      COMMON/ADD/SUM
      COMMON/INFO/NOCCO(2), NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMLOC/ISYMOFF(8,8,25)
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */
#include "flags.com"
#include "flags2.com"
      COMMON/METH/ MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON/LTRIP/LTRP
      COMMON/QRHFINF/POPRHF(8),VRTRHF(8),NOSH1(8),NOSH2(8),
c&two line mod
     &               POPDOC(8),VRTDOC(8),NAI,N1I,NA2,
     &               NUMISCF,NUMASCF,ISPINP,ISPINM,IQRHF
      COMMON/DERIV/DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,
     &              TRIP2,GABCD,RELAXED,TRULY_NONHF
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/ISTART/I0,ICRSIZ
CJDW KKB stuff
      COMMON /SHIFT/   ISHIFT
      COMMON /SYMDROP/ NDRPOP(8),NDRVRT(8)
      COMMON /DROPGEO/ NDRGEO
      COMMON /SYM/     POP(8,2),VRT(8,2),MTAA,MTBB,MD1AA,MD1BB,
     &                                             MD2AA,MD2BB
      COMMON /SYM2/    POP2(8,2),VRT2(8,2),NTAA,NTBB,ND1AA,ND1BB,
     &                                               ND2AA,ND2BB
      COMMON /SYMPOP/  IRP_DM(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMPOP2/ IRPDPD(8,22)
#include "parcc.com"
C
corig COMMON/SYM/POP(8,2),VRT(8,2),NTAA,NTBB,ND1AA,ND1BB,ND2AA,ND2BB
corig COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
C     At the beginning of a dropped core run, /SYM/, /SYMPOP/ will have
C     dropped core dimensions. /SYM2/ and /SYMPOP2/ are set in RECALAT.
C     When xdens is entered in dropped core runs, JOBARC will be JOBARC_DM.
C     Lists 1-200 or 1-290 will have dropped core dimensions. Allocation
C     of D, I, X is done in terms of NTAA, NTBB, etc, i.e. full MO space
C     dimensions.
CJDW END
C
      DATA HALF,ONE /0.5D0,1.0D0/
C
      sum=0.
C
      CALL CRAPSI(ICORE,IUHF,0)
      MAXCOR=ICRSIZ
      SYMMETRY = (iflags(h_IFLAGS_symmetry).eq.2)
C      If (SYMMETRY) THEN
C         Write(6,"(3a)") " This stupid code that print out density",
C     &              " contributions from each term does not"
C         Write(6,"(2a)") " with symmetry. Resubmit with symmetry",
C     &                   " turned off."
C         Write(6,*) 
C         Call Errex
C      Endif 
C
C FILL COMMON BLOCK /METH/
C
      CALL SETMET(IUHF)
      if (CCSD .OR. CCD .OR. TRIP1 .OR. TRIP2) Then

         call parread(iuhf)
         if (ispar) then
           If (TRIP1 .OR. TRIP2) Then
              Write(6,"(a,a)") " Analytical gradients for triple",
     &        "exciations pCC methods have not been verified"
              call aces_exit(1)
           Endif 
           write(6,*) ' Perform a parameterized CC Density calculations'
           Write(6,*)
           write(6,2010) paralpha
           write(6,2011) parbeta
           write(6,2012) pargamma
           write(6,2013) pardelta
           write(6,2014) parepsilon
 2010      format(' PCCSD   alpha parameter : ', F14.6)
 2011      format(' PCCSD    beta parameter : ', F14.6)
 2012      format(' PCCSD   gamma parameter : ', F14.6)
 2013      format(' PCCSD   delta parameter : ', F14.6)
 2014      format(' PCCSD epsilon parameter : ', F14.6)
           if (coulomb) Write(6,"(a,a)") " The Coulomb integrals are ",
     $                    "used in W(mbej) intermediate."
           write(6,*)
           Fae_scale    = (Paralpha - 1.0D0)
           Fmi_scale    = (Parbeta  - 1.0D0)
           Wmnij_scale  = Pargamma
           Wmbej_scale  = Pardelta
           Gae_scale    = Paralpha
           Gmi_scale    = Parbeta 
         else
           write(6,*) '  Perform a regular CC density calculations'
           write(6,*)
           write(6,*)
           Fae_scale    = 0.0D0
           Fmi_scale    = 0.0D0
           Wmnij_scale  = 1.0D0
           Wmbej_scale  = 1.0D0
           Gae_scale    = 1.0D0
           Gmi_scale    = 1.0D0

         endif
      endif

      bRedundant = (iflags2(h_IFLAGS2_single_store).eq.0)
C
C NEW FLAG FOR RHF AND UHF BASED MBPT(3)
C
      HFM3=MBPT3.AND.(.NOT.ROHF)
C
C SG 11/20/97 If a MRCC calculation, the density matrices are already on disk
C
      MRCC = IFLAGS2(h_IFLAGS2_program) .EQ. 3
      if (mrcc) then
         cis=.false.
      endif
      SEWARD = (IFLAGS(h_IFLAGS_integrals) .EQ.4)
C
C CALCULATE THE GAMMA INTERMEDIATES FOR THE VARIOUS METHODS
C
C NOTE THAT FOR MBPT(2) NO GAMMA INTERMEDIATES ARE REQUIRED, SO SKIP 
C THE CALL
C
C The recent noise about RELAX vs RESPONSE density convinced me that
C it is better to give the user the choice to specifiy what kind of
C density they want. The capabilities to generate the RESPONSE density
C matrix always existed in the code (vcceh: Ajith Perera; EOM-properties
C J. F. Stanton). Perhaps the seed for the complaints lies in the
C difficulties that one has to go through to get at the RESPONSE
C density matrix.
C
C The following call to gamma drive is needed only for RELAXED densities.
C If the perturbed orbitals are canonical then we need the call to gamma
C drive regardless. Also, first-order properties calculate only the
C necessary gamma elements for the relaxed density. Ajith Perera 08/2004.

      CALL GETREC(20,'JOBARC','NBASTOT ',1,NBAS)
      CALL RECALAT(IUHF,NIRREP,NDROP0,ICORE(i0),
     &             ICORE(i0+NBAS*(IUHF+1)),IINTFP,
     &             IFLAGS(h_IFLAGS_qrhf_gen))
C
      IF (IFLAGS(35).NE.0) CALL INCOR(I0,ICRSIZ,IUHF)
      MAXCOR=ICRSIZ

C ALLOCATE MEMORY FOR THE DENSITY MATRIX AT THE TOP OF ICORE
C
      MXCOR=MAXCOR

c Nevin 10/20/95 added for alignment
c      if(mod(mxcor,2).ne.0)mxcor=mxcor-1
C
C MEMORY FOR THE OCCUPIED-OCCUPIED BLOCK
C
      IDOO=i0+MXCOR-(ND1AA+IUHF*ND1BB)*IINTFP
      MXCOR=MXCOR-(ND1AA+IUHF*ND1BB)*IINTFP
C
C MEMORY FOR THE VIRTUAL-VIRTUAL BLOCK
C
      IDVV=IDOO-(ND2AA+IUHF*ND2BB)*IINTFP
      MXCOR=MXCOR-(ND2AA+IUHF*ND2BB)*IINTFP
C
C MEMORY FOR THE VIRTUAL-OCCUPIED BLOCK
C
      IDOV=IDVV-(NTAA+IUHF*NTBB)*IINTFP
      MXCOR=MXCOR-(NTAA+IUHF*NTBB)*IINTFP
C
C MEMORY FOR THE XOV INTERMEDIATE
C
C FOR NON HF WE HAVE TO ALLOCATE ADDITIONAL MEMORY FOR THE XOV INTERMEDIATE
C FOR HF CASES EQUIVALENCE MEMORY FOR XIA WITH THAT FOR DOV
C
      IF(QRHF.OR.NONHF.OR.ROHF) THEN
       IXOV=IDOV-(MAX(NTAA,NTBB)+IUHF*MAX(NTAA,NTBB))*IINTFP
       MXCOR=MXCOR-(MAX(NTAA,NTBB)+IUHF*MAX(NTAA,NTBB))*IINTFP
       LXOV=(1+IUHF)*MAX(NTAA,NTBB)
      ELSE
       IXOV=IDOV
       LXOV=0
      ENDIF
C
C MEMORY FOR THE IOV INTERMEDIATE
C
      IIOV=IXOV-(NTAA+IUHF*NTBB)*IINTFP
      MXCOR=MXCOR-(NTAA+IUHF*NTBB)*IINTFP
C
C MEMORY FOR THE IVO INTERMEDIATE
C
      IIVO=IIOV-(NTAA+IUHF*NTBB)*IINTFP
      MXCOR=MXCOR-(NTAA+IUHF*NTBB)*IINTFP
C
C MEMORY FOR THE IOO INTERMEDIATE
C
      IIOO=IIVO-(ND1AA+IUHF*ND1BB)*IINTFP
      MXCOR=MXCOR-(ND1AA+IUHF*ND1BB)*IINTFP
C
C MEMORY FOR THE IVV INTERMEDIATE
C
      IIVV=IIOO-(ND2AA+IUHF*ND2BB)*IINTFP
      MXCOR=MXCOR-(ND2AA+IUHF*ND2BB)*IINTFP
C
C SET ITOP SO THAT NON-HF WILL WORK
C
      ITOP=IIVV
C
C ZERO THE ALLOCATED SPACE FOR THE INTERMEDIATES IOO AND IVV
C
      LENGTH=LXOV+3*NTAA+2*ND1AA+2*ND2AA+
     &                IUHF*(3*NTBB+2*ND1BB+2*ND2BB)
C
C Programmers are setting the NONHF flag even when there might not
C be any non-HF terms. TRULY_NONHF is true only when the reference
C is ROHF, TWODET, or QRHF. For HF references, ANALYZE_FOCK checks
C for the existence of non-HF terms and sets NONHF_TERMS_EXIST
C accordingly.
C
      IF (.NOT.TRULY_NONHF .AND. NONHF) THEN
         CALL ANALYZE_FOCK(ICORE(IDOO), ICORE(IDVV), ICORE(IDOV),
     &                     ND1AA, ND2AA, NTAA, NONHF_TERMS_EXIST)
      ELSE
         NONHF_TERMS_EXIST = TRULY_NONHF
      END IF
C
      CALL IZERO(ICORE(IIVV),LENGTH*IINTFP)
C
C COMPUTE THE OCCUPIED-OCCUPIED BLOCK OF D
C      
       Write(6,"(a)") "D(I,J) and D(i,j) blocks"
       Write(6,"(a)") "------------------------"
       CALL DENSOO(ICORE(IDOO),ICORE(i0),MXCOR,IUHF)

C COMPUTE THE VIRTUAL-VIRTUAL BLOCK OF D
C      
       Write(6,"(a)") "D(A,B) and D(a,b) blocks"
       Write(6,"(a)") "------------------------"
       CALL DENSVV(ICORE(IDVV),ICORE(i0),MXCOR,IUHF)
C
C FOR NON HF CCSD COMPUTE DIRECT CONTRIBUTION TO THE OCCUPIED
C VIRTUAL BLOCK OF THE D, ROHF IS HERE TREATED AS NON-HF 
C
C The response density has only the direct VO contributions. The
C orbital relaxation terms for the relaxed density are calculated
C in MAKDIA1. At this point we have everything for the response
C density matrix. 09/2004, Ajith Perera.
C
       Write(6,"(2a)") "1/2(D(A,I)+D(I,A)) and /2(D(a,i)+D(i,a))",
     &                 " blocks"
       Write(6,"(2a)") "--------------------------------------------",
     &                "---"
       IF(NONHF_TERMS_EXIST.OR..NOT.RELAXED)THEN
         CALL DENSVO(ICORE(IDOV),ICORE(i0),MXCOR,IUHF,bRedundant)
       ENDIF
       Write(6,"(a)") "Upper case letters: Alpha spin"
       Write(6,"(a)") "Lower case letters: Beta  spin"
#ifdef _DEBUG_LVLM
       Write(6,*)
       Write(6,*) "The checksums of the Density blocks"
       Call checksum("DENSOO ", ICORE(IDOO),ND1AA+IUHF*ND1BB)
       Call checksum("DENSVV ", ICORE(IDVV),ND2AA+IUHF*ND2BB)
       Call checksum("DENSOV ", ICORE(IDOV),NTAA+IUHF*NTBB)
#endif 
C
      call aces_fin  
      STOP
      END 
