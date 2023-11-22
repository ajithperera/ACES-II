C
C VECTORIZED CC/MBPT PROGRAM FOR CALCULATING ONE- AND TWO-PARTICLE
C DENSITY MATRICES.
C
C THIS PROGRAM CALCULATES THE RELAXED DENSITY MATRIX FOR MBPT AND
C CC THEORY AS WELL COMPUTES THE QUANTITIES I(IA),I(IJ), AND I(AB)
C REQUIRED FOR GRADIENT CALCULATION
C IN ADDITION THE INTERMEDIATES GAMMA WERE FORMED AND STORED, THEY
C ARE USED TO CALCULATE THE TWO-ELECTRON-INTEGRAL DERIVATIVE              
C CONTRIBUTION TO THE GRADIENT
C THIS PROGRAM IS EFFICIENTLY VECTORIZED (BY WRITING ALL TERMS
C AS MATRIX PRODUCTS) AND USES EXPLICITELY SYMMETRY.
C 
C FOR NOTATION SEE :
C
C  J. GAUSS, J.F. STANTON, R.J. BARTLETT, J.C.P. 95 (1991) 2623
C
C AND ALSO :
C
C  E.A. SALTER, PH.D. THESIS, U.F. 1988
C  E.A. SALTER, G.W. TRUCKS, R.J. BARTLETT, J.C.P. 90 (1989) 1752
C
C THE IMPLEMENTATION IS MAINLY BASED ON THE THESIS OF SALTER AND 
C A RECENT FORMALISM DEVELOPED INDEPENDENTLY BY J. GAUSS IN COLOGNE
C (J. GAUSS AND D. CREMER, ADV. QUANT. CHEM. 23 (1992) 205). WITH
C RESPECT TO THE FORMALISM OF SALTER, THERE ARE SOME SUBTLE MODI-
C FICATION, MAINLY THE QUADRUPLE TERMS HAVE BEEN REFORMULATED AS
C FIRST SUGGESTED BY :
C
C  J. GAUSS AND D. CREMER ( CHEM. PHYS. LETT. 138 (1987) 131)
C
C THE DIFFERENCE TO THE ORIGINAL IMPLEMENTATION EITHER OF MBPT
C AS WELL AS OF CC GRADIENTS IS THAT THIS PROGRAM USES EXPLICITELY
C SYMMETRY AND IS VECTORIZED. ACTUALLY, THIS SHOULD BE THE MOST
C EFFICIENT IMPLEMENTATION WHICH CAN BE ACHIEVED ON MODERN SUPER
C COMPUTER.
C
C THIS PROGRAM HANDLES RHF, UHF, ROHF, AND QRHF REFERENCE FUNCTION.
C
C FOR ROHF REFERENCE FUNCTIONS SEE:
C
C  J. GAUSS, W.J. LAUDERDALE, J.F. STANTON, J.D. WATTS,
C  AND R.J. BARTLETT, C.P.L. 182 (1991) 207
C
C THE SOLUTION OF THE ROHF-CPHF Z-VECTOR EQUATIONS IN A 
C SEMICANONICAL REPRESENTATION HAS BEEN IMPLEMENTED AS AN
C ADDITIONAL OPTION:
C
C  W.J. LAUDERDALE, J.F. STANTON, J. GAUSS, J.D. WATTS,
C  R.J. BARTLETT, J.C.P. 97 (1992) 6606
C
C FOR ROHF-MBPT(3), SEE J. GAUSS, UNPUBLISHED NOTES.
C FOR QRHF REFERENCE FUNCTIONS SEE:
C
C  J. GAUSS, J.F. STANTON, AND R.J. BARTLETT, J.C.P. 95 (1991) 2639
C
C SOME FURTHER EXTENSIONS ARE DUE TO J.F. STANTON (UNPUBLISHED).
C
C FOR AN EFFICIENT TREATMENT OF TRIPLE SUBSTITUTIONS WE
C INCLUDED CANONICAL PERTURBED ORBITALS FOLLOWING A 
C SUGGESTION BY T.J. LEE AND A.P. RENDELL (J.CHEM.PHYS. 
C 94, 6229(1991).
C HOWEVER, REGARDING ALL OTHER ASPECTS TRIPLE EXCITATIONS 
C ARE TREATED AS FIRST SUGGESTED AND IMPLEMENTED BY GAUSS
C AND CREMER (C.P.L. 153 (1988) 303) AND INDEPENDENTLY BY 
C BARTLETT'S GROUP (C.P.L. ... (1988) ...). 
C APPROPRIATE REFERENCES FOR TRIPLES IN ACES II ARE:
C
C  J.D. WATTS, J. GAUSS, AND R.J. BARTLETT C.P.L. 200 (1992) 1
C  J.D. WATTS, J. GAUSS, AND R.J. BARTLETT J.C.P., 98 (1993) 8718
C
C GRADIENTS ARE ALSO AVAILABLE FOR UNITARY CC METHODS (P.SZALAY,
C UNPUBLISHED).
C
C DENSITIES FOR EXCITED STATE APPROACHES ARE OFFERED FOR CIS
C (TDA) IN FORM OF CIS-GRADIENTS AS WELL FOR EOM-CC APPROACHES
C 
C  J.F. STANTON, J.C.P. 99 (1993) XXXX
C  J.F. STANTON AND J. GAUSS J.C.P. 100 (1994) XXXX
C  J.F. STANTON AND J. GAUSS, TO BE PUBLISHED.
C 
C VARIOUS APPROACHES HAVE BEEN INCLUDED TO REDUCE DISK SPACE REQUIREMENTS.
C TOWARDS THE IMPLEMENTATION OF `DIRECT CC GRADIENTS' AS A FIRST STEP
C AO-BASED ALGORITHMS HAVE BEEN IMPLEMENTED (DEC. 93 -- FEBR. 94,
C J. GAUSS AND J.F. STANTON, UNPUBLISHED)
C
C TD-CC EXTENSION ADDED BY P.G. SZALAY (ROUTINES WITH 'OS')
C
C THIS PROGRAM HAS BEEN DEVELOPED AT THE QUANTUM THEORY PROJECT
C IN 1990-1992 BY J.F. STANTON AND J. GAUSS. EXTENSIONS HAVE BEEN
C MADE IN KARLSRUHE AS WELL IN AUSTIN.
C
      PROGRAM VDENS
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
      LOGICAL MRCC, SEWARD
      LOGICAL NONHF_TERMS_EXIST, STANDARD_NONCAN
      logical bRedundant
      LOGICAL pCCD
      LOGICAL LCC,OO_CC
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
           write(6,*)' Perform a parameterized CC Density calculations.'
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
           write(6,*) ' Perform a regular CC density calculations.'
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
      CALL Getrec(0,"JOBARC","PCCD_RUN",Len,Junk)
      pCCD = .False. 
      If (Len .Gt. 0) pCCD = .True.
C
C NEW FLAG FOR RHF AND UHF BASED MBPT(3)
C
      HFM3=MBPT3.AND.(.NOT.ROHF)
C
C SG 11/20/97 If a MRCC calculation, the density matrices are already on disk
C
      MRCC  = IFLAGS2(h_IFLAGS2_program) .EQ. 3
      LCC   = IFLAGS(h_IFLAGS_calc) .EQ. 5 .OR. 
     &        IFLAGS(h_IFLAGS_calc) .EQ. 6
      OO_CC = (IFLAGS(h_IFLAGS_calclevel).EQ.52 .OR.
     &         IFLAGS(h_IFLAGS_calclevel).EQ.53 .OR. 
     &         IFLAGS(h_IFLAGS_calclevel).EQ.54) 
      if (mrcc) then
         cis=.false.
      endif
      IF (PCCD .OR. OO_CC) CANON = .FALSE.
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
C
      IF (IFLAGS(35).NE.0) CALL INCOR(I0,ICRSIZ,IUHF)
      MAXCOR=ICRSIZ

C The following change is made to S^2 printed for response density
C matrices.
CSSS      IF (.NOT.MBPT2 .AND. .NOT.MRCC .AND. (RELAXED.OR.CANON)) THEN

      IF (.NOT.MBPT2 .AND. .NOT.MRCC) THEN
       CALL GAMDRV(ICORE(i0),MAXCOR,IUHF,LCC)
       IF (IFLAGS(35).NE.0) CALL INCOR(I0,ICRSIZ,IUHF)
       MAXCOR=ICRSIZ
      ENDIF
c
#ifdef _DEBUG_LVL0
CSSS      Call Checkgamma(Icore(I0),Icrziz,Iuhf)
CSSSS      CALL CHECKGAM1(ICORE(i0),16,116,TWO,IUHF,2,VRT)
CSSSS      CALL CHECKGAM1(ICORE(i0),133,233,TWO,IUHF,2,VRT)
CSSSS      call checkgam1(icore(i0),13,113,ONE,IUHF,1,POP)
CSSSS      CALL CHECKGAM(ICORE(i0),23,123,TWO)
CSSSS      CALL CHECKGAM(ICORE(i0),25,125,TWO)
CSSSS      CALL CHECKGAM(ICORE(i0),18,118,TWO)
CSSSS      call checkgam1(icore(i0),13,113,ONE,IUHF,1,POP)
CSSSS      call checkgam1(icore(i0),30,130,two,iuhf,2,vrt)
#endif 

CJDW KKB stuff
      CALL GETREC(20,'JOBARC','NBASTOT ',1,NBAS)
      CALL RECALAT(IUHF,NIRREP,NDROP0,ICORE(i0),
     &             ICORE(i0+NBAS*(IUHF+1)),IINTFP,
     &             IFLAGS(h_IFLAGS_qrhf_gen))
CJDW END
C
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
      IF (.NOT.CIS .AND. .NOT.EOM .AND. .NOT.MRCC)THEN
C
C COMPUTE THE OCCUPIED-OCCUPIED BLOCK OF D
C      
       CALL DENSOO(ICORE(IDOO),ICORE(i0),MXCOR,IUHF)

       IF(IFLAGS(1).GE.15.AND.NDRGEO.EQ.0)THEN
        write(6,*) ' D intermediates after densoo '
        CALL DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IIVO),ICORE(IIOV),
     &              IUHF,"AL")
       ENDIF
C
C COMPUTE THE VIRTUAL-VIRTUAL BLOCK OF D
C      
       CALL DENSVV(ICORE(IDVV),ICORE(i0),MXCOR,IUHF)
       IF(IFLAGS(1).GE.15.AND.NDRGEO.EQ.0)THEN
        write(6,*) ' D intermediates after densvv '
        CALL DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IIVO),ICORE(IIOV),
     &              IUHF,"AL")
       ENDIF
C
C FOR NON HF CCSD COMPUTE DIRECT CONTRIBUTION TO THE OCCUPIED
C VIRTUAL BLOCK OF THE D, ROHF IS HERE TREATED AS NON-HF 
C
C The response density has only the direct VO contributions. The
C orbital relaxation terms for the relaxed density are calculated
C in MAKDIA1. At this point we have everything for the response
C density matrix. 09/2004, Ajith Perera.
C
       IF(NONHF_TERMS_EXIST.OR..NOT.RELAXED)THEN
        CALL DENSVO(ICORE(IDOV),ICORE(i0),MXCOR,IUHF,bRedundant)
        IF(IFLAGS(1).GE.15.AND.NDRGEO.EQ.0)THEN
         write(6,*) ' D intermediates after densvo '
         CALL DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IDOV),ICORE(IIOV),
     &               IUHF,"AL")
        ENDIF
       ENDIF
#ifdef _DEBUG_LVL0
       Write(6,*)
       Write(6,*) "The checksums of the Density blocks"
       Call checksum("DENSOO ", ICORE(IDOO),ND1AA+IUHF*ND1BB)
       Call checksum("DENSVV ", ICORE(IDVV),ND2AA+IUHF*ND2BB)
       Call checksum("DENSOV ", ICORE(IDOV),NTAA+IUHF*NTBB)
       Write(6,*)
       write(6,*) ' Density matrices  '
       CALL PCCD_DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IDOV),
     &                 ICORE(IIOV),IUHF,"AL")
       Write(6,*)
#endif 
C
      ELSE
C
C FOR CIS (TDA), GET DENSITY MATRIX FROM GAMLAM FILE
C THE CONTRIBUTIONS ARE CALCULATED IN VEE
C
       IF(CIS)CALL GETDCIS(ICORE(IDOO),ICORE(IDVV),IUHF)
       IF (EOM .OR. MRCC) THEN
        CALL GETD(ICORE(IDOO),ICORE(IDOO+IINTFP*MD1AA),
     &           ICORE(IDVV),ICORE(IDVV+IINTFP*MD2AA),
     &           ICORE(IDOV),ICORE(IDOV+IINTFP*MTAA),
     &           IUHF,ROHF.OR.QRHF.OR.NONHF)
       ENDIF
C
       IF(IFLAGS(1).GE.15.AND.NDRGEO.EQ.0)THEN
        write(6,*) ' D intermediates after getd '
        CALL DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IDOV),ICORE(IIOV),
     &              IUHF,"AL")
       ENDIF
C
      ENDIF
C
C GENERATE FOR ALL METHODS EXCEPT MBPT3 AND MBPT2 HERE THE GAMMA1
C AMPLITUDES. ROHF-MBPT(3) REQUIRES CALL TO GAMMA1.
C
      IF (.NOT.CIS .AND. .NOT.MBPT2 .AND. .NOT.HFM3

C Then following change is made to get S^2 printed for response
C densities 
CSSS     &   .AND. .NOT.MRCC  .AND. (RELAXED.OR.CANON)) THEN
     &   .AND. .NOT.MRCC) THEN
#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(a)") "Checksums of G(ij,ab)"
#endif
       CALL GAMMA1(ICORE(i0),MXCOR,IUHF)
      ENDIF
#ifdef _DEBUG_LVL0
CSSS      Call Checkgamma(Icore(I0),Icrsiz,Iuhf)
#endif 
C 
CJDW KKB stuff
c-------------------------------------------------------------------
      ISHIFT = 0 
      ISHIFT2= 0 
      if (ndrgeo.ne.0) then 
       ISHIFT = 300  
       if (iuhf.ne.0) then 
        call Icopy(MD1BB, ICORE(IDOO+IINTFP*MD1AA),1,
     &                    ICORE(IDOO+IINTFP*ND1AA),1)
        call expand(ICORE(IDOO+IINTFP*ND1AA),ICORE(i0),NIRREP,4,2,ND1BB)
        call Icopy(MD2BB,ICORE(IDVV+IINTFP*MD2AA),1,
     &                   ICORE(IDVV+IINTFP*ND2AA),1)
        call expand(ICORE(IDVV+IINTFP*ND2AA),ICORE(i0),NIRREP,1,2,ND2BB)
        call Icopy (MTBB, ICORE(IDOV+IINTFP*MTAA),1,
     &                    ICORE(IDOV+IINTFP*NTAA),1)
        call expand(ICORE(IDOV+IINTFP*NTAA),ICORE(i0),NIRREP,3,2,NTBB) 
       endif
       call expand(ICORE(IDOO),ICORE(i0),NIRREP,4,1,ND1AA) 
       call expand(ICORE(IDVV),ICORE(i0),NIRREP,1,1,ND2AA) 
       call expand(ICORE(IDOV),ICORE(i0),NIRREP,3,1,NTAA) 
       IF(IFLAGS(1).GE.15.AND.NDRGEO.EQ.0)THEN
        write(6,*) ' D intermediates after dmove and expand'
        CALL DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IDOV),ICORE(IIOV),
     &              IUHF,"AL")
       ENDIF
      endif 
C
      IF(NDROP0.NE.0)THEN
        CALL PCHANGE(IUHF,NIRREP,1,IFLAGS(h_IFLAGS_qrhf_gen))
      ENDIF
c-------------------------------------------------------------------
CJDW END
C
C At this point we have everything we need for first-order property
C calculations. Let's write the RESPONSE density matrix to the disk
C and quit. Note that "RESPONSE" density is possible only for first
C order property calculations. Note that the geometry optimizations
C automaticaly turns on the "RELAXED" density option. I hate this
C GO TO statements. This code can be cleaned up bit more. Ajith Perera
C
C OLD
C     IF(NONHF) GO TO 1500
C NEW
C I don't see the merits of the JDW's jump to 1500 for NONHF methods.
C This JDW jump correspond to strange logic in setmet (see line 232)
C where NONHF is set to .TRUE. when REFERENCE  flag is set to ROHF
C and NONHF flag is on  (very strange logic!!). What was acomplished
C by doing all this is to generate response density matrix for NONHF
C density matrix calculations. It can be acomplished now without having
C to go through this strange logic. 09/04, Ajith Perera.
      IF (.NOT.RELAXED) GO TO 1500
C END
C
C FOR GRADIENT CALCULATION CALCULATE ALSO I(IJ) AND I(AB)
C FOR NONHF METHODS THESE TERMS ARE NEEDED ANYWAY, FOR 
C GRADIENT AND PROPERTY CALCULATIONS
C
C OLD
C     IF(GRAD.OR.QRHF)THEN
C NEW
c GRAD was always true
      IF (.TRUE.) THEN
C END
C
C LOAD I INTERMEDIATES CALCULATED IN EOM CODE FOR "DIRECT"
C G(ABCD) ALGORITHM
C
CSSS      IF(EOM.AND..NOT.GABCD)THEN
C GABCD=.TRUE. only when GAMMA_ABCD is set to DIRECT. What that implies
C is not to store G(AB,CD) (known as G2) on disk. Prior twisted logic
C that set GABCD was confusing and now it is readable and makes sense.
C
      IF (EOM.AND.GABCD) THEN
       DO 50 ISPIN=1,1+IUHF
        CALL GETLST(ICORE(IIVV+IINTFP*(ISPIN-1)*ND2AA),1,1,1,ISPIN,92)
        CALL GETLST(ICORE(IIOV+IINTFP*(ISPIN-1)*NTAA),1,1,1,ISPIN,93)
50     CONTINUE
      ENDIF
C
C NOTE THE INTERMEDIATES ARE CALCULATED HERE IN ORDER TO
C BE ABLE TO RESORT SOME OF THE GAMMA INTERMEDIATES TO
C A NEW ORDERING AS REQUIRED FOR THE CALCULATION OF IOV
C
C FIRST TERM IS REQUIRED FOR ALL METHODS
C
      IF(IFLAGS(1).GE.15)THEN
       write(6,*)
       write(6,*) ' i intermediates before INTO0 '
       CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"OO")
      ENDIF
C
C JCP-95-2623-91, Eqn. 43, 3-6 terms when (p,q = i,j; INTO1;(r,s,t = m,e,f
C INTO2; (r,s,t = m,n,o), INTO3; (r,s,t = m,e,f), INTO4 (r,s,t = e,f,g)
C INTO5; (r,s,t = m,n,e), INTO6; (r,s,t = m,n,e). Note that when r,s, t
C is m,e,f or m,n,e the i and j can permute to give two distinct terms.
C That is the reason behind INTO1,INTO3 and INTO5, INTO6 involve same
C r,s,t indices. 09/2004, Ajith Perera
C
      IF(.NOT.CIS)CALL INTO1(ICORE(IIOO),ICORE(i0),MXCOR,IUHF)
      IF(IFLAGS(1).GE.15)THEN
       write(6,*) ' i intermediates after  INTO1 '
       CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"OO")
      ENDIF
C
C SECOND AND THIRD TERM ARE REQUIRED FOR MBPT3 AND CCD
C WHILE ALL OTHER METHODS REQUIRE ALL TERMS
C
      IF(.NOT.MBPT2)THEN
       IF(.NOT.CIS)CALL INTO2(ICORE(IIOO),ICORE(i0),MXCOR,IUHF)
      IF(IFLAGS(1).GE.15)THEN
       write(6,*) ' i intermediates after  INTO2 '
       CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"OO")
      ENDIf
       CALL INTO3(ICORE(IIOO),ICORE(i0),MXCOR,IUHF,bRedundant)
      IF(IFLAGS(1).GE.15)THEN
       write(6,*) ' i intermediates after  INTO3 '
       CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"OO")
      ENDIF
       IF(.NOT.CIS.AND..NOT.HFM3.AND..NOT.CCD) THEN
C Recognize that ABCDTYPE and GAMMA_ABCD are two different flags. In
C principal one can use them independently, but that makes little sense
C (if you can store ABCD, then it does not matter whether it is GAMMA or
C bare ABCD). So, we need to be sensible here and when the user invokes
C one, then the other should get set by default. The comment below
C computing both occupied-occupied and virtual-occupied contributions
C together in one place for direct algorithms is correct (see below).
C Ajith Perera, 05/05.
        IF(IFLAGS(h_IFLAGS_abcdtype).NE.2) THEN
C
C IF, THE <AB||CD> CONTRIBUTION TO XIA IS CALCULATED
C WITH AO INTEGRALS, SKIP INTO4. THIS CONTRIBUTION
C IS THEN A BY-PRODUCT IN XINT6AO
C
         CALL INTO4(ICORE(IIOO),ICORE(i0),MXCOR,IUHF)
      IF(IFLAGS(1).GE.15)THEN
       write(6,*) ' i intermediates after  INTO4-MO '
       CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"OO")
      ENDIF
        ENDIF
        CALL INTO5(ICORE(IIOO),ICORE(i0),MXCOR,IUHF)
      IF(IFLAGS(1).GE.15)THEN
       write(6,*) ' i intermediates after  INTO5 '
       CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"OO")
      ENDIF
        CALL INTO6(ICORE(IIOO),ICORE(i0),MXCOR,IUHF)
      IF(IFLAGS(1).GE.15)THEN
       write(6,*) ' i intermediates after  INTO6 '
       CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"OO")
      ENDIF
       ENDIF
      ENDIF
C
C IN THE CASE OF CANONICAL PERTURBED ORBITALS, SKIP
C THE CALL TO DXAINI (IT IS CALLED AFTER X(I,J) HAS
C BEEN FORMED AND THE FULL BLOCKS OF D(I,J) AND D(A,B)
C ARE AVAILABLE.
C
C This is the 7th term in the Eqn. 43 of  JCP-95-2623-91.
C For perturbative triples perturb orbitals are set to be canonical
C and this term looks different (JCP-98-8718-93 Eqn. 70) and they
C are calculated below in the routines FORMXIJ, FORMXAB, ADDIA1, ADDIA2,
C XAINT1, XINT2, and DXAINI. 09/2004, Ajith Perera.
C
      IF(.NOT.CANON) THEN
       CALL DXAINI(ICORE(IDOO),ICORE(IDOV),ICORE(IDVV),
     &             ICORE(IIOO),ICORE(i0),MXCOR,IUHF)
      ENDIF
C
C FIRST TERM REQUIRED FOR ALL METHODS
C
C JCP-95-2623-91, Eqn. 43, 3-6 terms when (p,q = a,b; INTV1;(r,s,t = m,n,e
C INTV2; (r,s,t = e,f,g), INTV3; (r,s,t = m,n,e,), INTV4 (r,s,t = m,e,f)
C INTV5; (r,s,t = m,e,f), INTV6; (r,s,t = m,n,o). Note that when r,s,t
C is m,e,f or m,n,e the a and b can permute to give two distinct terms.
C That is the reason behind INTV1,INTO3 and INTV4, INTV5 involve same
C r,s,t indices.
C
      IF(.NOT.CIS)CALL INTV1(ICORE(IIVV),ICORE(i0),MXCOR,IUHF)
      IF (IFLAGS(1).GE.15) THEN
         write(6,*)
         write(6,*) ' I intermediates after INTV1 '
         CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &               ICORE(IIOV),IUHF,"VV")
       ENDIF
C
C SECOND AND THIRD TERM ARE REQUIRED FOR MBPT3 AND CCD
C WHILE ALL OTHER METHODS REQUIRE ALL TERMS
C
      IF(.NOT.MBPT2)THEN
       IF(.NOT.CIS) THEN
C
C CHOOSE BETWEEN THE VARIOUS ALGORITHMS FOR THE INTV2 TERM
C Choose the direct route when G(ab,cd) is not on disk: GABCD=.TRUE.
C
        IF(.NOT.GABCD) THEN
         CALL INTV2(ICORE(IIVV),ICORE(i0),MXCOR,IUHF)
         IF (IFLAGS(1).GE.15) THEN
            write(6,*) ' I intermediates after INTV2MO '
            CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"VV")
         ENDIF
        ELSE
         CALL DRVAOVV(ICORE(IIVV),ICORE(I0),MXCOR,IUHF,
     &                0,1,143,63,243,260)
         CALL REWORK63(ICORE(I0),MXCOR,IUHF)
         IF (IFLAGS(1).GE.15) THEN
            write(6,*) ' I intermediates after INTV2-DRVAOVV '
            CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"VV")
         ENDIF
        ENDIF
       ENDIF
C
       CALL INTV3(ICORE(IIVV),ICORE(i0),MXCOR,IUHF,bRedundant)
       IF (IFLAGS(1).GE.15) THEN
           write(6,*) ' I intermediates after INTV3 '
           CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"VV")
       ENDIF
       IF(.NOT.CIS.AND..NOT.HFM3.AND..NOT.CCD) THEN
        CALL INTV4(ICORE(IIVV),ICORE(i0),MXCOR,IUHF)
         IF (IFLAGS(1).GE.15) THEN
            write(6,*) ' I intermediates after INTV4'
            CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"VV")
         ENDIF
        CALL INTV5(ICORE(IIVV),ICORE(i0),MXCOR,IUHF)
         IF (IFLAGS(1).GE.15) THEN
            write(6,*) ' I intermediates after INTV5 '
            CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"VV")
         ENDIF
        CALL INTV6(ICORE(IIVV),ICORE(i0),MXCOR,IUHF)
         IF (IFLAGS(1).GE.15) THEN
            write(6,*) ' I intermediates after INTV6 '
            CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"VV")
         ENDIF
       ENDIF
      ENDIF
      ENDIF
C
C GO ON WITH DENSITY STUFF
C
C COMPUTE THE INTERMEDIATE I(IA)
C
C JCP-95-2623-91, Eqn. 43, 3-6 terms when (p,q = i,a; INTOV1;(r,s,t = m,e,f
C INTOV2; (r,s,t = e,f,g), INTOV3; (r,s,t = m,n,e,), INTOV4 (r,s,t = m,e,f)
C INTOV5; (r,s,t = m,e,f), INTOV6; (r,s,t = m,n,o). Note that when r,s,t
C is m,e,f or m,n,e the a and b can permute to give two distinct terms.
C That is the reason behind INTV1,INTO3 and INTV4,INTV5 involve same
C r,s,t indices. 09/2004, Ajith Perera
C
      IF(.NOT.CIS)CALL INTOV1(ICORE(IIOV),ICORE(i0),MXCOR,IUHF)
      IF (IFLAGS(1).GE.15) THEN
          write(6,*)
          Write(6,*)' I intermediates after INTOV1 '
          CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                      ICORE(IIOV),IUHF,"OV")
      ENDIF 
C
C FOR MBPT(3) AND CCD IN ADDITION THE SECOND AND THIRD
C TERM ARE REQUIRED, FOR ALL OTHER METHODS ALL TERMS ARE REQUIRED
C
      IF(.NOT.MBPT2)THEN
C
C RESORT FIRST THE G(IA,JB) AMPLITUDES
C
       CALL SORTGAM(ICORE(i0),MXCOR,IUHF)
CJDW KKB stuff
C
C     Note that this is always called. The value of 2 means that /SYM/ is
C     copied to /SYM2/.
C
      CALL PCHANGE(IUHF,NIRREP,2,IFLAGS(h_IFLAGS_qrhf_gen))
CJDW END
       IF(.NOT.CIS) THEN
C Choose the direct route to when G(AB,CD) is not on disk:GABCD=.TRUE.
        IF(.NOT.GABCD) THEN
         CALL INTOV2(ICORE(IIOV),ICORE(i0),MXCOR,IUHF)
         IF (IFLAGS(1).GE.15) THEN
            Write(6,*)' I intermediates after INTOV2MO '
            CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"OV")
         ENDIF
        ELSE
         CALL DRVAOOV(ICORE(IIOV),ICORE(I0),MXCOR,IUHF)
         IF (IFLAGS(1).GE.15) THEN
            Write(6,*)' I intermediates after INTOV2-AO'
            CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"OV")
         ENDIF
        ENDIF
       ENDIF
       CALL INTOV3(ICORE(IIOV),ICORE(i0),MXCOR,IUHF)
       IF (IFLAGS(1).GE.15) THEN
           Write(6,*)' I intermediates after INTOV3'
           CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                 ICORE(IIOV),IUHF,"OV")
       ENDIF
       IF(.NOT.CIS.AND..NOT.HFM3.AND..NOT.CCD) THEN
        CALL INTOV4(ICORE(IIOV),ICORE(i0),MXCOR,IUHF)
       IF (IFLAGS(1).GE.15) THEN
           Write(6,*)' I intermediates after INTOV4'
           CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                 ICORE(IIOV),IUHF,"OV")
       ENDIF
        CALL INTOV5(ICORE(IIOV),ICORE(i0),MXCOR,IUHF,bRedundant)
       IF (IFLAGS(1).GE.15) THEN
           Write(6,*)' I intermediates after INTOV5'
           CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                 ICORE(IIOV),IUHF,"OV")
       ENDIF
        CALL INTOV6(ICORE(IIOV),ICORE(i0),MXCOR,IUHF)
       IF (IFLAGS(1).GE.15) THEN
           Write(6,*)' I intermediates after INTOV6'
           CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                 ICORE(IIOV),IUHF,"OV")
       ENDIF
       ENDIF
      ENDIF
C
C COMPUTE THE INTERMEDIATE X(IA)
C
C NOTE THE FIRST AND SECOND CONTRIBUTION ARE FOR ALL METHODS IDENTICAL
C
C JCP-95-2623-91, The 43 the 7th term. When the perturb orbitals
C are set to canonical (ROHF non iterative triples), the Eqn, 70
C of JCP-98-8718-93 applies(skip at this point to do latter).
C A very important to note that what we are calculating is I(a,i)
C contirbutions. They will eventually turned into the X(a,i) intermediate.
C The XINT1; (r,s = m,n), XINT2; (r,s = e,f). 09/2004, Ajith Perera.
C
      IF(.NOT.CANON) THEN
       CALL XINT1(ICORE(IIVO),ICORE(IDOO),ICORE(i0),MXCOR,IUHF)
       IF (IFLAGS(1).GE.15) THEN
          Write(6,*)
          Write(6,*)' I intermediates after XINT1'
          CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                ICORE(IIOV),IUHF,"VO")
       ENDIF
       CALL XINT2(ICORE(IIVO),ICORE(IDVV),ICORE(i0),MXCOR,IUHF)
       IF (IFLAGS(1).GE.15) THEN
          Write(6,*)' I intermediates after XINT2'
          CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                ICORE(IIOV),IUHF,"VO")
       ENDIF
      ENDIF
C
C JCP-95-2623-91, The 43 the 7th term. Only when D(O,V) is involoved.
C XINT9 (r,s = e,m). The direct contribution of D(O,V) is applicable
C only when F(a,i) is non zero.
C
      IF(NONHF_TERMS_EXIST) THEN
       CALL XINT9(ICORE(IIVO),ICORE(IDOV),ICORE(i0),MXCOR,IUHF,
     &            bRedundant)
       IF (IFLAGS(1).GE.15) THEN
          Write(6,*)' I intermediates after XINT9'
          CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                ICORE(IIOV),IUHF,"VO")
       ENDIF
      ENDIF
C
C  FOR MBPT(3) AND CCD IN ADDITION THE FOURTH AND FIFTH
C  TERM ARE REQUIRED, FOR ALL OTHER METHODS ALL TERMS ARE REQUIRED
C
C The formation of X(a,i) intermediates. See Eqns. 49 in  JCP-95-2623-91
C followed by Eqn. 43. The X(a,i) is needed for Z-vector method
C for orbital relaxations.
C
C JCP-95-2623-91, Eqn. 43, 3-6 terms when (p,q = a,i; XINT4;(r,s,t = m,n,o
C XINT5; (r,s,t = n,e,f), XINT6; (r,s,t = e,f,g), XINT7 (r,s,t = m,n,e)
C XINT8; (r,s,t = m,n,e), XINT3; (r,s,t = m,e,f). Notice the symmetry
C in XINT5, XINT3 and XINT4 and XINT7. 09/2004, Ajith Perera.
C
      IF(.NOT.MBPT2)THEN
       IF(.NOT.CIS)CALL XINT4(ICORE(IIVO),ICORE(i0),MXCOR,IUHF)
       IF (IFLAGS(1).GE.15) THEN
          Write(6,*)' I intermediates after XINT4'
          CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                ICORE(IIOV),IUHF,"VO")
       ENDIF
       CALL XINT5(ICORE(IIVO),ICORE(i0),MXCOR,IUHF)
       IF (IFLAGS(1).GE.15) THEN
          Write(6,*)' I intermediates after XINT5'
          CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                ICORE(IIOV),IUHF,"VO")
       ENDIF
       IF(.NOT.CIS.AND..NOT.HFM3.AND..NOT.CCD) THEN
C
C CHOOSE BETWEEN THE TWO POSSIBLE ALGORITHM FOR THE XINT6 CONTRIBUTION
C (I.E., AO OR MO BASED)
C
C Compute both occupied-occupied and virtual-occupied contributions to I
C from G(ABCI) (G6) when <AB||CD> is treated directly. See the comments
C before the call to INTO4. Ajith Perera, 05/05.
C
        IF (IFLAGS(h_IFLAGS_abcdtype).EQ.2) THEN
         CALL XINT6AO(ICORE(IIVO),ICORE(IIOO),ICORE(i0),MAXCOR,1,IUHF)
         IF (IFLAGS(1).GE.15) THEN
            Write(6,*) ' XI intermediates after XINT6AO'
            CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"VO")
         ENDIF
        ELSE 
         CALL XINT6(ICORE(IIVO),ICORE(i0),MXCOR,IUHF)
         IF (IFLAGS(1).GE.15) THEN
            Write(6,*)  ' XI intermediates after  XINT6MO'
            CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"VO")
         ENDIF
        ENDIF
C
        CALL XINT7(ICORE(IIVO),ICORE(i0),MXCOR,IUHF,bRedundant)
       IF (IFLAGS(1).GE.15) THEN
          Write(6,*)' I intermediates after XINT7'
          CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                ICORE(IIOV),IUHF,"VO")
       ENDIF
        CALL XINT8(ICORE(IIVO),ICORE(i0),MXCOR,IUHF)
        IF (IFLAGS(1).GE.15) THEN
          Write(6,*)' I intermediates after XINT8'
          CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                ICORE(IIOV),IUHF,"VO")
        ENDIF
       ENDIF
      ENDIF
C
C THE THIRD CONTRIBUTION IS REQUIRED FOR ALL METHODS (EXCEPT CIS)
C
      IF(.NOT.CIS) then 
        CALL XINT3(ICORE(IIVO),ICORE(i0),MXCOR,IUHF)
        IF (IFLAGS(1).GE.15) THEN
          Write(6,*)' I intermediates after XINT3'
          CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                ICORE(IIOV),IUHF,"VO")
        ENDIF
      endif
C
C ADD FINALLY THE F(P,R)*D(Q,R) CONTRIBUTION TO THE
C INTERMEDIATES
C
C This term is JCP-95-2623-91, Eqn. 43's first two
C terms for (p,q = a,b; i,j; a,i) and what is calculated here
C is contributions to the I(a,b), I(i,j) and I(i,a). 09/2004, Ajith Perera.
C
      CALL DFINI(ICORE(IDOO),ICORE(IDVV),ICORE(IDOV),ICORE(IIOO),
     &           ICORE(IIVV),ICORE(IIOV),ICORE(IIVO),ICORE(i0),
     &           MXCOR,IUHF,NONHF_TERMS_EXIST)

      IF(IFLAGS(1).GE.15)THEN
       write(6,*)
       write(6,*) ' i intermediates after DFINI '
       CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"AL")
      ENDIF
      IF(IFLAGS(1).GE.10)THEN
       write(6,*) ' D matrices '
       CALL DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"AL")
      ENDIF
C
C IN THE CASE OF CANONICAL PERTURBED ORBITALS, FORM HERE
C X(IJ) AND X(AB) AS WELL AS THE FINAL I(IJ) AND I(AB).
C SOLVE FOR D(IJ) AND D(AB) AND AUGMENT X(AI)
C
C Once again the perturb canonical orbitals are choosen for
C ROHF non-iterative triples along with the semicanonical option.
C What that does is to slightly change the formulas for I and as
C a result X intermediates.  09/2004, Ajith Perera.
C
      IF(CANON) THEN
C
       IEVAL=i0
       IXOO=IEVAL+IINTFP*(NOCCO(1)+NVRTO(1))
C
C
      IF(IFLAGS(1).GE.15)THEN
       Write(6,*) "I, D, X before formxij"
       CALL DUMPIT(ICORE(IIOO),ICORE(IDVV),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"AL")
       CALL DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"AL")
       CALL DUMPIT(ICORE(IXOO),ICORE(IEVAL),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"AL")
      ENDIF
C
C The Eqn. 66 of JCP-98-8718-93, first term. Form the D(i,j) contribution.
C
       CALL FORMXIJ(ICORE(IIOO),ICORE(IDOO),
     &              ICORE(IXOO),ICORE(IEVAL),IUHF,QRHF)
      IF(IFLAGS(1).GE.15)THEN
       Write(6,*) "I, D, X after formxij"
       CALL DUMPIT(ICORE(IIOO),ICORE(IDVV),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"AL")
       CALL DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"AL")
       CALL DUMPIT(ICORE(IXOO),ICORE(IEVAL),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"AL")
      ENDIF
C
C The Eqn. 70 of JCP-98-8718-93, part of the second contribution.
C A contribution to the I(a,i) intermediate. Once again notice that
C only so called non-HF methods contribute.
C
       IF(NONHF_TERMS_EXIST) THEN
        IFAI=IXOO+IINTFP*(ND1AA+ND1BB)
        CALL ADDIA1(ICORE(IXOO),ICORE(IIVO),ICORE(IFAI))
       ENDIF
C
       IXVV=IXOO
C
C The Eqn. 66 of JCP-98-8718-93, second  term. Form the D(a,b) contribution
C
       CALL FORMXAB(ICORE(IIVV),ICORE(IDVV),
     &              ICORE(IXVV),ICORE(IEVAL),IUHF,QRHF)
C
C The Eqn. 70 of JCP-98-8718-93, part of the second contribution.
C A contribution to the I(a,i) intermediate. Once again notice that
C only so called non-HF methods contribute.
C
       IF(NONHF_TERMS_EXIST) THEN
        IFAI=IXVV+IINTFP*(ND2AA+ND2BB)
        CALL ADDIA2(ICORE(IXVV),ICORE(IIOV),ICORE(IFAI))
       ENDIF
C
C These are the contributions to the I(a,i) intermediate that we
C did not included previously (see above). The reason is that the
C form of D(i,j) and D(a,b) is different (see Eqn. 66 in JCP-98-8718-93).
C Since now we have the correct form of D(i,j) and D(a,b) (after
C calls to FORMXIJ and FORMXAB, we can do the now. What is calculated
C is contribution to X(a,i) intermediate via as in Eqn. 70, the
C second term in the first Equation. The XINT1; (r,s = m,n),
C XINT2; (r,s = e,f).
C
       IF(IFLAGS(1).GE.15)THEN
       Write(6,*) "I, and D before XINT1 and XINT2"
       CALL DUMPIT(ICORE(IIOO),ICORE(IDVV),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"AL")
       CALL DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"AL")
       ENDIF

       CALL XINT1(ICORE(IIVO),ICORE(IDOO),ICORE(i0),MXCOR,IUHF)
       CALL XINT2(ICORE(IIVO),ICORE(IDVV),ICORE(i0),MXCOR,IUHF)

       IF(IFLAGS(1).GE.15)THEN
       Write(6,*) "I, and D after XINT1 and XINT2"
       CALL DUMPIT(ICORE(IIOO),ICORE(IDVV),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"AL")
       CALL DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"AL")
       ENDIF
C
C  ADD THE CONTRIBUTION DUE TO DOO AND DVV TO IOO
C
C This is the 7th term in the Eqn. 43 of JCP-95-2623-91.
C As discussed above, now we have the correct form of the
C D(i,j) and D(a,b) for perturb orbitals canonical (FORMXIJ,
C FORMXAB), so we can do them here now.
C
       CALL DXAINI(ICORE(IDOO),ICORE(IDOV),ICORE(IDVV), 
     &             ICORE(IIOO),ICORE(i0),MXCOR,IUHF,NONHF_TERMS_EXIST)
C
      ENDIF
C
      IF(IFLAGS(1).GE.10)THEN
       Write(6,*) "Intermediates before solving Z vec eqns"
       CALL DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),ICORE(IIOV),
     &             IUHF,"AL")
      ENDIF

      IF (OO_CC .OR. PCCD) THEN
         CALL PUTREC(20,"JOBARC","ROTGRDOO",ND1AA+IUHF*ND1BB,
     &               ICORE(IIOO))
         CALL PUTREC(20,"JOBARC","ROTGRDVV",ND2AA+IUHF*ND2BB,
     &               ICORE(IIVV))
         CALL PUTREC(20,"JOBARC","ROTGRDOV",NTAA+IUHF*NTBB,
     &               ICORE(IIOV))
         CALL PUTREC(20,"JOBARC","ROTGRDVO",NTAA+IUHF*NTBB,
     &               ICORE(IIVO))
      ENDIF 
C
C FORM X(A,I) REQUIRED FOR SOLVING THE Z-VECTOR EQUATIONS
C
C The relevent paper is JCP-95-2639-91 and the equations are 15 and
C 16 for the QRHFX. The equations 24 applies for QRHFZ. Beware that
C QRHFX and QRHFZ are doing specilized tasks. Unless you know
C exactly what you doing do not turn on the QRHF flag. All the
C other methods, the equation 48 of JCP-95-2623-91 applies. Note
C the memory location IIVO has the I(a,i) contributions and the
C and the  IIOV has the I(i,a) contributions. Now, form the
C I(i,a) - I(a,i) and copy it to the IXOV memory location. A comment
C about the test with NONHF_TERMS_EXIST. If true what does this tells
C is that the user has only specified NONHF flag, but analysis of the
C Fock matrix revels that it indeed a NONHF type fock matrix. Then
C follow the usual path. If that is not the case, user is doing
C something strange but not invalid. In that case, memory for the D(a,i)
C (direct contributions is not required). Since memory is allocated
C based on NONHF=.TRUE. before the analysis of the Fock matrix, proper
C location for IXOV is IDOV (see the memory allocation section for D(a,i) at the
C begining). Obviously, if this happen memory allocated for D(a,i) is simply a
C waste, but it is such small amount and can be safely ignored.
C 09/2004, Ajith Perera.
C
      IF(QRHF) THEN
C
C FOR NON-HF WAVE FUNCTION THIS IS A LITTLE BIT COMPLICATED
C NONHFX EXTRACTS THE RELEVANT TERMS FOR QRHF-CCSD CALCULATIONS
C
C ALLOCATE FIRST ADDITIONAL MEMORY FOR THE X1A AND X1B ARRAYS
C
C
       IX1=IIVV-2*IINTFP*N1I
       MXCOR=MXCOR-IINTFP*2*N1I
       IX2=IX1-2*IINTFP*NA2
       MXCOR=MXCOR-IINTFP*2*NA2
       ITOP=IX2
       CALL IZERO(ICORE(IX1),2*IINTFP*N1I)
       CALL IZERO(ICORE(IX2),2*IINTFP*NA2)
C
       CALL QRHFX(ICORE(IIOO),ICORE(IIVV),ICORE(IIOV),ICORE(IIVO),
     &            ICORE(IIOO+IINTFP*ND1AA),ICORE(IIVV+IINTFP*ND2AA),
     &            ICORE(IIOV+IINTFP*NTAA),ICORE(IIVO+IINTFP*NTAA),
     &            ICORE(IXOV),ICORE(IXOV+IINTFP*NAI),ICORE(IX1),
     &            ICORE(IX1+IINTFP*N1I),ICORE(IX2),
     &            ICORE(IX2+IINTFP*NA2),ICORE(i0))
C
C  FORM HERE NOW THE Z1I (OR ZA1) CONTRIBUTION TO THE DENSITY
C  AND CALCULATE XAI(TWIDDLE)
C
       IEV=ITOP-IINTFP*(NOCCO(1)+NVRTO(1))
       CALL QRHFZ(ICORE(IXOV),ICORE(IXOV+IINTFP*NAI),ICORE(IX1),
     &            ICORE(IX1+IINTFP*N1I),ICORE(IX2),
     &            ICORE(IX2+IINTFP*NA2),ICORE(IEV),ICORE(i0),
     &            MXCOR,IUHF)
C
      ELSE
       IF (.NOT.NONHF_TERMS_EXIST) IXOV=IDOV
       CALL DCOPY((NTAA+IUHF*NTBB),ICORE(IIVO),1,ICORE(IXOV),1)
       CALL SAXPY((NTAA+IUHF*NTBB),ONE,ICORE(IIOV),1,ICORE(IXOV),1)
       ITOP=IIVV
      ENDIF
C
C  SOLVE THE Z-VECTOR EQUATION WITHIN CPHF THEORY
C
C Now we have the correct X(a,i) at the correct memory location.
C So, we are ready to compute the orbital relaxation contribution
C to the virtual-occupied block of the density matrix. For the QRHF
C we need to make sure that elements relevent to the special orbital
C (s?) is taken care of for the appropriate blocks of the density
C matrices. That is done in QRHFD. When the QRHFD is finished
C everything is done for the density calculations and the FRMFUL
C with 'REL" as last argument write out the NHFDENSA and NHFDENSB records.
C The I intermediates for the gradients (GRAD only) is built according to
C the Eqn. 29 and 30 in JCP-95-2623-91 paper. If it is not a QRHF (any
C other reference) the Eqns. 53,54 and 55 in JCP-95-2623-91 in applies.
C They are computed in DAIINI. Aftre DAIINI the I intermediate is
C ready for use. For ROHF cases we also need to add the the orbital
C relaxation terms to the occ-vrt block of the density matrix.
C 09/2004, Ajith Perera
C
C Our density code was not able to handle non-canonical (better say
C f(ab),f(ij) off-diagnoal are non zero). The new MAKDIA1_4STD_NONCAN
C routine correctly handle Z-vector solving for those situations. 
C Note that there are changes to how X(a,i) is built. They are
C sperad out in various call (DFINI in particular). With these 
C changes both UHF and RHF standard non-canonical situations are 
C treated correctly. Since the NONHF_TERMS_EXIST is determined 
C by analyzing the Fock matrix no additional key-words are
C not necessary but NONHF=ON is necessary for VCC and Lambda
C codes. 04/2014. Ajith perera. 
C
      IF (NONHF_TERMS_EXIST .AND. .NOT. TRULY_NONHF) THEN 
 
          Write(6,"(2x,a,a)")"Warning!! Enter into exclusively reserved"
     &                     ," block that is only for non-digonal"
          Write(6,"(12x,a)") "elements in FOCK but not ROHF or QRHF." 

          CALL MAKDIA1_4STD_NONCAN(ICORE(IXOV),ICORE(i0),MXCOR,IUHF,
     &                             bRedundant)
      ELSE
          CALL MAKDIA1(ICORE(IXOV),ICORE(i0),MXCOR,IUHF,bRedundant)
      ENDIF
C
C  FOR HF CASES MAKDIA1 RETURNS THE OCCUPIED-VIRTUAL BLOCK OF THE
C  DENSITY MATRIX AT ADDRESS  IDOV (=IXOV), FOR NON-HF CASES 
C  IDOV NE IXOV
C
C  FORM NOW FINAL INTERMEDIATES AND SAVE QUANTITIES
C
      IF(QRHF) THEN
C
C  FOR NON-HF CASES THE INDIRECT CONTRIBUTION (ICORE(IXOV)) 
C  HAS TO BE TREATED SEPARATELY FROM THE DIRECT CONTRIBUTION
C  SINCE FOR GRADIENT RUNS THEY HAVE TO CONTRACTED WITH
C  DIFFERENT FOCK MATRIX DERIVATIVES
C
C  FORM FULL UNPACKED NON-DIRECT PART OF THE MO DENSITY MATRIX
C  AND WRITE IT TO DISK
C
       MXCOR1=MXCOR
       NBAS=NOCCO(1)+NVRTO(1)
       IDFULL=ITOP-2*NBAS*NBAS*IINTFP
       MXCOR=MXCOR-2*NBAS*NBAS*IINTFP
       IDROO=IDFULL-(ND1AA+ND1BB)*IINTFP
       IDRVV=IDROO-(ND2AA+ND2BB)*IINTFP
       IDROV=IDRVV-(NTAA+NTBB)*IINTFP
       CALL IZERO(ICORE(IDROV),IINTFP*(ND1AA+ND1BB+ND2AA
     &            +ND2BB+NTAA+NTBB))
       MXCOR=MXCOR-(ND1AA+ND1BB+ND2AA+ND2BB+NTAA+NTBB)*IINTFP
C
C  DISTRIBUTE NOW THE XOV AND X1 ELEMENTS CALCULATED AMONG THE 
C  VARIOUS BLOCKS OF THE DENSITY MATRIX
C
C  THIS IS HOWEVER A BOTHER !
C
       CALL QRHFD(ICORE(IXOV),ICORE(IXOV+IINTFP*NAI),ICORE(IX1),
     &            ICORE(IX1+IINTFP*N1I),ICORE(IX2),
     &            ICORE(IX2+IINTFP*NA2),ICORE(IDROO),
     &            ICORE(IDROO+IINTFP*ND1AA),ICORE(IDRVV),
     &            ICORE(IDRVV+IINTFP*ND2AA),ICORE(IDROV),
     &            ICORE(IDROV+IINTFP*NTAA),ICORE(i0))
       CALL FRMFUL(ICORE(IDROO),ICORE(IDRVV),ICORE(IDROV),ICORE(IDFULL),
     &             ICORE(i0),MXCOR,NBAS,IUHF,'REL')
C
C  TO CALCULATE THE REMAINING CONTRIBUTION TO I(P,Q), IT
C  IS HERE ADVANTAGEOUS TO CONSTRUCT THE SPIN-SYMMETRIZED DENSITY
C  MATRIX.  THE ADDRESSES BELOW REFER TO THIS QUANTITY, WHICH IS 
C  COMPUTED WITHIN QRHFI.  
C  
       IDSOO=IDROV-IINTFP*(ND1AA+ND1BB)
       IDSVV=IDSOO-IINTFP*(ND2AA+ND2BB)
       IDSOV=IDSVV-IINTFP*(NTAA+NTBB)
       MXCOR=MXCOR-IINTFP*(ND1AA+ND1BB+ND2AA+ND2BB+NTAA+NTBB)
C
       CALL QRHFI (ICORE(IIOO),ICORE(IIVV),ICORE(IIOV),
     &             ICORE(IDROO),ICORE(IDRVV),ICORE(IDROV),
     &             ICORE(IIOO+IINTFP*ND1AA),ICORE(IIVV+IINTFP*ND2AA),
     &             ICORE(IIOV+IINTFP*NTAA),ICORE(IDROO+IINTFP*ND1AA),
     &             ICORE(IDRVV+IINTFP*ND2AA),ICORE(IDROV+IINTFP*NTAA),
     &             ICORE(IDSOO),ICORE(IDSVV),ICORE(IDSOV),      
     &             ICORE(IDSOO+IINTFP*ND1AA),ICORE(IDSVV+IINTFP*ND2AA),
     &             ICORE(IDSOV+IINTFP*NTAA),
     &             ICORE(i0),ICORE(i0),MXCOR,
c&new line
     &             ICORE(IIVO),ICORE(IIVO+IINTFP*NTAA))
C
C RESET CORE MEMORY
C
       MXCOR=MXCOR1
C
      ELSE
C
C  NOTE, THAT DAINI TAKES CARE OF ALL THE ROHF STUFF
C
       CALL DAIINI(ICORE(IXOV),ICORE(IIOO),ICORE(IIOV),ICORE(i0),
     &             MXCOR,IUHF,NONHF_TERMS_EXIST)
C
C FOR ROHF, ADD DIRECT AND IN-DIRECT CONTRIBUTIONS
C
      IF(ROHF .OR. NONHF_TERMS_EXIST) THEN
        IF (IUHF .NE. 0) THEN
            CALL SAXPY(NTAA+NTBB,ONE,ICORE(IXOV),1,ICORE(IDOV),1)
        ELSE
            CALL SAXPY(NTAA,ONE,ICORE(IXOV),1,ICORE(IDOV),1)
        ENDIF 
      ENDIF
C
      ENDIF
C
C
C  FORM FULL UNPACKED MO-DENSITY AND WRITE IT TO DISK
C
C PRESERVE MXCOR SINCE THIS MXCOR FOR ALL LATER STEPS
C
1500  CONTINUE
C
C The following two FRMFUL calls built and write the complete
C density matrix and the I intermediates, depending on the last
C argument. The argument 'DEN' will write out RELDENSA and RELDENSB records.
C The argument INT will write out the I intermediate. The I intermediate
C is required only for the gradient calculations. 09/2004, Ajith Perera.
C
#ifdef _DEBUG_LVL0
       Write(6,*) 
       Write(6,*) "The checksums of the Density blocks"
       Call checksum("DENSOO ", ICORE(IDOO),ND1AA+IUHF*ND1BB)
       Call checksum("DENSVV ", ICORE(IDVV),ND2AA+IUHF*ND2BB)
       Call checksum("DENSOV ", ICORE(IDOV),NTAA+IUHF*NTBB)
#endif
      MXCOR1=MXCOR
      NBAS=NOCCO(1)+NVRTO(1)
      IDFULL=ITOP-(1+IUHF)*NBAS*NBAS*IINTFP
      MXCOR=MXCOR-(1+IUHF)*NBAS*NBAS*IINTFP
      CALL FRMFUL(ICORE(IDOO),ICORE(IDVV),ICORE(IDOV),
     &            ICORE(IDFULL),ICORE(i0),MXCOR,NBAS,IUHF,'DEN')
cmn
c Add density correction if STEOM calculation
c
      CALL GETREC(20,'JOBARC','STEOM_SZ',1,NSTEOM)
      if (nsteom .ne. 0) then
         IF (NSTEOM .NE. NBAS*NBAS) THEN
            WRITE(6,*) ' SOMETHING WRONG STEOM CORRECTION IN DENS'  
            CALL ERREX
         ENDIF
      ISTEOM=IDFULL-(1+IUHF)*NBAS*NBAS*IINTFP
      MXCOR=MXCOR-(1+IUHF)*NBAS*NBAS*IINTFP
      CALL GETREC(20,'JOBARC','STEOM_D ',NBAS*NBAS*IINTFP,
     $     ICORE(ISTEOM))
cmn      call output(icore(idfull), 1, nbas, 1, nbas, nbas, nbas, 1)
      CALL SAXPY(NBAS*NBAS, ONE, ICORE(ISTEOM), 1, ICORE(IDFULL), 1)   
cmn      call output(icore(idfull), 1, nbas, 1, nbas, nbas, nbas, 1)
      endif
C
C  FORM FULL UNPACKED INTERMDIATES AND WRITE IT TO DISK
C
C PRESERVE MXCOR SINCE THIS MXCOR FOR ALL LATER STEPS
C
      IF (GRAD) THEN
      MXCOR=MXCOR1
      NBAS=NOCCO(1)+NVRTO(1)
      IIFULL=ITOP-(1+IUHF)*NBAS*NBAS*IINTFP
      MXCOR=MXCOR-(1+IUHF)*NBAS*NBAS*IINTFP
      CALL FRMFUL(ICORE(IIOO),ICORE(IIVV),ICORE(IIOV),
     &            ICORE(IIFULL),ICORE(i0),MXCOR,NBAS,IUHF,'INT')
      ENDIF

#ifdef _DEBUG_LVLM
      Call check_1body_dens(Icore(Idoo),Icore(Idvv),Icore(Idov),
     &                      Icore(Iioo),Icore(Iivv),Icore(Iiov),
     &                      Icore(I0),Mxcor,nbas,Iuhf)
#endif 
C
C DIAGONALIZE DENSITY AND FORM NATURAL ORBITALS IF REQUESTED.
C Marcel wants naturla occupation numbers always printed. 
C 
cSSS      IF(IFLAGS(h_IFLAGS_rdo).NE.0)THEN
       I000=i0
       I010=I000+NBAS*NBAS*IINTFP
       MEMLEFT=I0+MAXCOR-I010
       CALL NATORB(ICORE(I000),ICORE(I010),NBAS,MEMLEFT,IUHF,NDROP0)
cSSS      ENDIF
C
C OLD CODE, IF G(AB,CD) IS NOT CALCULATED, IT IS CALCULATED
C ON THE FLY IN ANTI.
C
C IF G(AB,CD) IS NOT CALCULATED, COMPUTE HERE 
C AND STORE IT ON LIST 231 -- 233
C
c      IF(IFLAGS(h_IFLAGS_props).NE.1.AND.(.NOT.GABCD)) THEN
C
C CREATE LISTS 231 -- 233 IF THEY ARE NOT THERE     
C
c       IF(IFLAGS(h_IFLAGS_abcdtype).EQ.2) THEN
c        IMODE=0
c        CALL INIPCK(1,13,13,233,IMODE,0,1)
c        IF(IUHF.NE.0) THEN
c         CALL INIPCK(1,1,1,231,IMODE,0,1)
c         CALL INIPCK(1,2,2,232,IMODE,0,1)
c        ENDIF
c       ENDIF
C 
c       IF(EOM)CALL EOMDABCD(ICORE(i0),MAXCOR,IUHF)
c       CALL GAMMA2(ICORE(i0),MAXCOR,IUHF,100)
C
c      ENDIF
C
C FOR OPEN-SHELL CCSD CALCULATIONS, COMPUTE SPIN MULTIPLICITY
C
      IF(NDRGEO.EQ.0)THEN
C
C Extending the S2EXPECT to work when GABCD=DIRECT is rather
C straightforward. At this point I am not going to spend time on it
C since the application of S2EXPECT results are rather limited. AP 5/5
C The conditional RELAXED is need because two-particle density matrix
c is not created for response first-order properties. 
C Ajith Perera, 04/2013. 
C
      IF (.NOT. EOM .AND..NOT.GABCD .AND. CCSD .AND. (.NOT. LTRP) .AND.
C The following change 
CSSS     &    IFLAGS(h_IFLAGS_ref).NE.0 .AND. RELAXED) THEN
     &    IFLAGS(h_IFLAGS_ref).NE.0) THEN
       CALL S2EXPECT(ICORE(i0),MAXCOR,IUHF,'R',.FALSE.,0)
      ENDIF
C
      ENDIF
C
C FOR CC AND MBPT(4) SECOND DERIVATIVE CALCULATIONS, SCALE H4 BY A
C FACTOR OF 1/2. FOR CC, TRANSPOSE THE H4 ARRAY.
C
      IF(IFLAGS(h_IFLAGS_deriv_lev).EQ.2) THEN
C
       CALL FIXH4(ICORE(i0),MAXCOR,IUHF)
C
      ENDIF
C
      WRITE(6,3001)
3001  FORMAT('  Density calculation successfully completed.')
c-------------------------------------------
      if (ndrgeo.ne.0) then 
       call resetscf(iuhf,nirrep,ndrop0,ICORE(i0),IINTFP)
      endif 
c-------------------------------------------
C For Molcas gradient calculations let's built the AO
C basis density matrix and the I(p,q) intermediates here. For
C vdint based analytical derivative calculations these steps are
C done internally. The total SCF density is already in the JOBARC file (see 
C vscf.f, blt_efctvfcok.f in vscf module). Ajith Perera 07/2000.
C Lets built these records anyway to get density plot codes to 
C work. These records should not interfere with others. 
C Ajith Perera, 05/2015. 
C
CSSS      IF (SEWARD) THEN
         IF (NDRGEO .EQ.0 .AND. .NOT. PCCD) CALL BLTAO_QUNTITES
     +                                           (ICORE(i0),
     +                                            MAXCOR, IUHF, 
     +                                            NONHF)
CSSS      ENDIF
C
      IF (OO_CC .OR. PCCD) THEN
          CALL FORM_HTAU(ICORE(I0),MAXCOR,IUHF)
      ENDIF 

      IF (IFLAGS(35).NE.0) THEN
         CALL ACES_AUXCACHE_FLUSH
         CALL ACES_AUXCACHE_RESET
      END IF
      call aces_fin  
      STOP
      END 
