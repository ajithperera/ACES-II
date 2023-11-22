C
C VECTORIZED PROGRAM TO SOLVE THE LAMBDA EQUATIONS 
C WITHIN CC/MBPT THEORY
C
C TD-CC EXTENSION ADDED BY P.G. SZALAY (ROUTINES WITH 'OS')
C       (SEE SZALAY AND BARTLETT, JCP 101, XXX (1994))
C
C EXTENDED IN FEBRUARY 1991 TO CALCULATE IN ADDITION
C HBAR FOR EOM-CC OR FOCK-SPACE CC METHODS.
C
C IN JUNE, 1991, THE HBAR-CALCULATING ROUTINES ARE PUT INTO THIS PRODUCTION
C VERSION.
C
C WRITTEN BY J. GAUSS AND J.F. STANTON
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */

      PROGRAM VLAMCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER RLECYC,POP,VRT,DIRPRD,DISTSZ
c& new line
      LOGICAL OSFLAG
      LOGICAL CC,SING,IHBAR,IHBAR2,IHBAR3,ROHF
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,CC2
      LOGICAL TERM1, bredundant,TERM2,TERM3,TERM4,TERM5,TERM6
      LOGICAL TRIPNI1,TRIPNI2,TRIPIT,T3STOR,PRESNT,INIT
      LOGICAL HBAR_4LCCSD, DO_TAU
      LOGICAL ACT_SPC_CC,SING_4TCC
      LOGICAL pCCD,pCCDS,pCCDTS,pCCDTSD
      LOGICAL LCC
      CHARACTER*4 ACT
      DIMENSION ECORR(3)
      COMMON // ICORE(1)
      COMMON /ISTART/I0,ICRSIZ
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /FLAGS/IFLAGS(100)
      COMMON /FLAGS2/IFLAGS2(500)
      COMMON /HBAR/IHBAR,IHBAR2,IHBAR3
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &             CC2
      COMMON /REF/ROHF
cAP - 500 is the maximum number of CC iterations (not basis functions)
      COMMON /ENERGY/ ENERGY(500,2),IXTRLE(500)
c&one line modif
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /SYMLOC/ ISYMOFF(8,8,25)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
C
c&new lines
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
C
      INTEGER FSDPDAN,FSDPDNA,FSDPDAA,FSDPDIN,FSDPDNI,FSDPDII,FSDPDAI,
     $   FSDPDIA
      COMMON /FSSYMPOP/ FSDPDAN(8,22),FSDPDNA(8,22),FSDPDAA(8,22),
     &                  FSDPDIN(8,22),FSDPDNI(8,22),FSDPDII(8,22),
     &                  FSDPDAI(8,22),FSDPDIA(8,22)
C
      COMMON /ACTIRR/ IRREPOA,IRREPOB,IRREPVA,IRREPVB
c&end
C
      COMMON /TRIPLES/ TRIPNI1,TRIPNI2,TRIPIT,T3STOR
      COMMON /T3IOOF/ IJKPOS(8,8,8,2),IJKLEN(36,8,4),IJKOFF(36,8,4),
     1                NCOMB(4)
      COMMON /AUXIO / DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     1                NRECS,LUAUX
C
#include "parcc.com"
C
      EQUIVALENCE (IFLAGS(2),METHOD)
C
      DATA ONE,ONEM /1.0D0,-1.0D0/
C
      WRITE(6,1000)
C
C  INITIALIZE CRAPS
C 
      CALL CRAPSI(ICORE,IUHF,0)
      MAXCOR=ICRSIZ

C Set pCCD and LCC (LCCD and LCCSD) flags 

      pCCD    = IFLAGS(h_IFLAGS_calc) .EQ. 52
      LCC     = IFLAGS(h_IFLAGS_calc) .EQ. 5 .OR. 
     &         IFLAGS(h_IFLAGS_calc) .EQ. 6
      pCCDS   = IFLAGS(h_IFLAGS_calc) .EQ. 55
      pCCDTS  = IFLAGS(h_IFLAGS_calc) .EQ. 56
      pCCDTSD = IFLAGS(h_IFLAGS_calc) .EQ. 58
      If (pCCD) THEN
         METHOD=8
      ELSEIF (pCCDS) THEN
         METHOD=10
      ELSEIF (pCCDTS) THEN
         METHOD=10
      ELSEIF (pCCDTSD) THEN
         METHOD=10
      ENDIF
      bredundant = iflags2(155).eq.0
C
C  Pseudoenergy printing flag (for CMPENG2).
C
      IF (IFLAGS(1).GE.10) THEN
         IPRTENG = 1
      ELSE
         IPRTENG = 0
      END IF
C
C  FILL /METH/ WHICH DETERMINES THE METHOD
C
      CALL SETMET(pCCD,pCCDS,pCCDST)

      if (method .eq. 10 .or. method .eq. 8
     $     .or. method .eq. 22 .or. method  .eq. 13
     $     .or. method .eq. 14 .or. method  .eq. 16) then
         coulomb = .false. 
         call parread(iuhf)
         if (ispar) then
           write(6,*) ' Parameterized lambda calculation is performed '
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
           write(6,*) ' Regular lambda calculation is performed'
           write(6,*)
           Fae_scale    = 0.0D0
           Fmi_scale    = 0.0D0
           Wmnij_scale  = 1.0D0
           Wmbej_scale  = 1.0D0
           Gae_scale    = 1.0D0
           Gmi_scale    = 1.0D0
         endif
      endif
#ifdef _DEBUG_LVL0
      Write(6,"(2a,1x,3(1x,L))") " @-Vlamcc, Hbar control variables ",
     +                           " IHBAR,IHBAR2,IHBAR2:",
     +                             IHBAR,IHBAR2,IHBAR3
#endif 
c
      if (ispar .and. (.not. bredundant) .and.
     $     abs(parbeta-1.0d0) .gt. 1.0d-6) then
         write(6,*) ' PCCSD not supported with NO_REDUNDANT option on '
         write(6,*) ' Need to fix subroutine WMBEJNR for this to work'
         write(6,*) ' PCCSD Calculation will still work if parbeta=1.0'
         call aces_exit(1)
      endif
C
C  INITIALIZE SOME LAMBDA LISTS
C
      CALL INILAM(IUHF)
C
C ADDITIONAL STUFF FOR THE LOGIC DEALING WITH THE METHOD
C
      CC=CCD.OR.CCSD.OR.QCISD.OR.CC2
      SING=M4SDQ.OR.M4SDTQ.OR.QCISD.OR.CCSD.OR.ROHF.OR.CC2
c& new lines
      CALL GETREC(-1,'JOBARC','OSCALC  ',1,IOS)
      IF(IOS.EQ.0) THEN
         OSFLAG=.FALSE.
      ELSE
         OSFLAG=.TRUE.
      ENDIF
c& new line end
C
C-----------------------------------------------------------------------
C      Open triples file if necessary.
C-----------------------------------------------------------------------
      T3STOR = METHOD.EQ.17 .OR. METHOD.EQ.18                       .OR.
     &         ( METHOD.EQ.22 .AND. IFLAGS(87).EQ.3 .AND.
     &          (IFLAGS2(124).EQ. 4 .OR. IFLAGS2(124).EQ. 6 .OR.
     &           IFLAGS2(124).EQ. 8 .OR. IFLAGS2(124).EQ.10 .OR.
     &                                   IFLAGS2(124).EQ.12) )      .OR.
     &   ( (METHOD.EQ.14 .OR. METHOD.EQ.16 .OR. METHOD.EQ.33
     &                                     .OR. METHOD.EQ.34) .AND.
     &      IFLAGS(87).EQ.3                                       ) .OR.
     &   (  METHOD.EQ.13 .AND. IFLAGS(87).EQ.3 .AND. IFLAGS(91).GT.0 )
C
      IF(T3STOR)THEN
        IRREPX = 1
        IF( (METHOD.GE.14.AND.METHOD.LE.16) .OR.METHOD.EQ.33 .OR.
     &                                          METHOD.EQ.34    )THEN
          PRESNT = .TRUE.
          INIT   = .FALSE.
          NLIST  = 12
        ELSEIF(METHOD.EQ.13)THEN
          PRESNT = .FALSE.
          INIT   = .TRUE.
          NLIST  = 12
        ELSEIF(METHOD.EQ.22)THEN
          PRESNT = .TRUE.
          INIT   = .FALSE.
          NLIST  =  4
        ELSEIF(METHOD.EQ.17 .OR. METHOD.EQ.18)THEN
          PRESNT = .TRUE.
          INIT   = .TRUE.
          NLIST  =  8
        ENDIF
        CALL AUXIOI(PRESNT,INIT,IRREPX,NLIST,IUHF)
      ENDIF
C-----------------------------------------------------------------------
      IF(.NOT.IHBAR)THEN
c& new lines
C
C     INITIALIZATION FOR OPENSHELL SINGLETT CALCULATION
C
       IF(OSFLAG) THEN
          CALL OSLINIC(ICORE,ICRSIZ,IUHF,i0)
          MAXCOR=ICRSIZ
       ENDIF
C
c& end
C
C INITIALIZE STUFF FOR RLE CONVERGENCE ACCELERATION
C
       IF(CC) THEN
        CALL ZERO(ENERGY,1000)
        CALL IZERO(IXTRLE,500)
       ENDIF
C
C SET UP SOME INFORMATION BEFORE CALCULATION.
C
#ifdef _TAILORED
      act_spc_cc = .true.
      act = "READ"
      Write(6,"(a)") "First update of L1/L2 lists 100+90(1,2),44,45,46"
      SING_4TCC = .TRUE.
      SING_4TCC = .NOT.(IFLAGS(h_IFLAGS_calc) .EQ. 8 .OR.
     &                  IFLAGS(h_IFLAGS_calc) .EQ. 5 .OR.
     &                  pCCD .OR. pCCDTS .OR. pCCDTSD)
      if (act_spc_cc .and. (act .eq."READ" .or. act .eq."ZERO"))
     &   call modf_amps(ICORE(I0),MAXCOR,IUHF,0,SING_4TCC,'L',
     &                  act,"AMPS")
#endif
       IEVERY=5
       IF(CC)THEN
        NCYCLE=IFLAGS(7)
        NKEEP=IFLAGS(12)
        ICONTL=IFLAGS(4)
       ENDIF
       IF(IFLAGS(1).GE.10)THEN
        WRITE(6,801)
        CALL AMPSUM(ICORE(i0),MAXCOR,IUHF,100,SING,'L')
       ENDIF
801    FORMAT('  Initial lambda amplitudes: ')
       CALL CMPENG2(ICORE(i0),MAXCOR,143,0,ECORR,ENERGY(1,1),
     &              ENERGY(1,2),IUHF,IPRTENG)
       CALL CMPENG2(ICORE(i0),MAXCOR, 60,0,ECORR,ENERGY(1,1),
     &              ENERGY(1,2),IUHF,IPRTENG)

      ENDIF
C
C The above ENDIF is for IHBAR.
C

      IF (IFLAGS(35).NE.0) CALL INCOR(I0,ICRSIZ,IUHF)
      MAXCOR=ICRSIZ
C
C BLOCK OF CODE FOR THIRD-ORDER ROHF AND FOURTH-ORDER CALCULATIONS.
C
      IF(.NOT.CC .AND. .NOT. LCC)THEN
C
C THIRD-ORDER LOGIC: COMPUTE LAMBDA(2) AMPLITUDES
C (ONLY REQUIRED FOR ROHF REFERENCE FUNCTIONS)
C
       IF(ROHF) THEN
C
C FOR FOURTH-ORDER, RESORT FIRST T2[1] AMPLITUDES ON 33--39
C
         IF(.NOT.MBPT3) CALL RNABIJ(ICORE(i0),MAXCOR,IUHF,'T')
C
C CALCULATE MISSING CONTRIBUTIONS TO SECOND-ORDER LAMBDAS
C
          CALL ROM3L(ICORE(i0),MAXCOR,IUHF)
C 
C FOR ROHF-MBPT(4), CALCULATE THIRD-ORDER LAMBDA
C
        IF(.NOT.MBPT3) THEN
c         CALL ROM4L(ICORE(i0),MAXCOR,IUHF)
        ENDIF
C
C ALL IS DONE
C
        WRITE(6,1010)
        IF (IFLAGS(35).NE.0) THEN
           CALL ACES_AUXCACHE_FLUSH
           CALL ACES_AUXCACHE_RESET
        END IF
        call aces_fin
        STOP
       ENDIF

C The above ENDIF ends the RPHF only logic.
C
C FOURTH-ORDER LOGIC : COMPUTE DELTA T3(SDQ) AMPLITUDES
C
C  THE DELTA T3(Q) CONTRIBUTION IS ALREADY ON THE L2 ( = T2) INCREMENT LIST
C  CALCULATE ONLY DELTA T3(D) AND DELTA T3(S) HERE
C
C   FIRST DENOMINATOR "DEWEIGHT" THE QUAD CONTRIBUTION
C   (SHOULD BE FIXED IN A LATER VERSION)
C
       CALL DEWQ(ICORE(i0),MAXCOR,IUHF)
c       IF(IUHF.EQ.1) THEN
c       CALL CHECKGAM(ICORE,44,61)
c       CALL CHECKGAM(ICORE,45,62)
c       ENDIF
c       CALL CHECKGAM(ICORE,46,63)
C	
C  RESORT FIRST LAMBDA AMPLITUDES
C
       if(bredundant) then 
          CALL RNABIJ(ICORE(i0),MAXCOR,IUHF,'L')
       endif 
C
C CALCULATE THE DELTA T3(D) CONTRIBUTION
C
       CALL L2INL2(ICORE(i0),MAXCOR,IUHF,0)
C
C CALCULATE THE DELTA T3(S) CONTRIBUTION
C
c       IF(IUHF.EQ.1) THEN
c       CALL CHECKGAM(ICORE,44,61)
c       CALL CHECKGAM(ICORE,45,62)
c       ENDIF
c       CALL CHECKGAM(ICORE,46,63)
       IF(M4SDQ.OR.M4SDTQ) THEN
        CALL L1INL2A(ICORE(i0),MAXCOR,IUHF,.TRUE.)
        CALL L1INL2B(ICORE(i0),MAXCOR,IUHF,.TRUE.)
       ENDIF
C
C CALCULATE THE NEW DELTA T3(SDQ) CONTRIBUTION
C
       CALL NEWT2(ICORE(i0),MAXCOR,IUHF)
C
C CALCULATE THE ENERGY, IT SHOULD EQUAL E4(S) + E4(D) + E4(Q)
C
       CALL CMPENG2(ICORE(i0),MAXCOR,60,0,ECORR,ENERGY(2,1),
     &              ENERGY(2,2),IUHF,IPRTENG)
C
       CALL FORMT2(ICORE(i0),MAXCOR,IUHF,43,143)
C
C GENERATE W LISTS FOR INTERMEDIATES
C
       CALL STLLST(ICORE(i0),MAXCOR,IUHF)
C
C GENERATE V AND G INTERMEDIATES, REQUIRED FOR THE X(IJ,AB) CONTRIBUTION
C
        CALL LGENINT(ICORE(i0),MAXCOR,IUHF)

C Ends the perturbation theory logic
C
      ELSE
C 
C Lets compute the Hbar for linear CCSD in a seperate routine. 
C We have not yet call genint in CC, so the lists are simply
C in the state of LCCSD. 
   
      HBAR_4LCCSD =  (IFLAGS(2) .EQ. 6 .OR. IFLAGS2(117) .EQ. 7)

CSSSS      IF (HBAR_4LCCSD) CALL MAKE_LCCSD_HBAR(ICORE(i0),MAXCOR,IUHF)

C
C CC CALCULATION LOGIC BEGINS HERE.
C
C RESORT FIRST LAMBDA AMPLITUDES 
C
       if(bredundant) then
          CALL RNABIJ(ICORE(i0),MAXCOR,IUHF,'L')
       endif
c&new lines
C
C RESORT ALSO T2 LAMBDA AMPLITUDES 
C
       IF(OSFLAG .and. bRedundant)CALL RNABIJ(ICORE(i0),MAXCOR,IUHF,'T')
c&new lines end
C
C GENERATE W LISTS FOR INTERMEDIATES
C
       CALL STLLST(ICORE(i0),MAXCOR,IUHF)
C
C CALCULATE FOR CCSD THE W(AB,EF) INTERMEDIATES
C
C
C CALCULATE MODIFIED W(MBEJ) INTERMEDIATE
C
CSSS#ifdef _DCC_FLAG
       If (Ispar) Then
          Write(6,"(a)") "            -     =               "
          Write(6,"(a)") " ***The DCC W and W is constructed" 
          Call Dcc_fixes_4lambda(Icore(i0),Maxcor,Iuhf)
CSSS#else
       Else
          CALL FORMWL(ICORE(i0),MAXCOR,IUHF) 
       Endif 
CSSS#endif 
C
C CALCULATE W(IFMN) AND W(EFAM) INTERMEDIATES FOR CCSD (see the
C comments in setmet.F for RPA and DRPA (CCSD is set for those
C as well. Ajith Perera, 07/2013. 
C CC2 option is added in 04/2015. See also the changes in WTWTW.
C Note that IFLAGS(38).EQ.0 means HF methods (f(a,i) = 0.0)
C
       IF(CCSD .OR. CC2 .OR. CCD)THEN

CSSS#ifndef _DCC_FLAG
       If (.NOT. Ispar) CALL WTWTW(ICORE(i0),MAXCOR,IUHF)
CSSS#endif 

#ifdef _DEBUG_LVLM
        call check_wps(Icore(i0),Maxcor,iuhf)
#endif 
        TERM1 = .TRUE.
        TERM2 = .TRUE.
        TERM3 = .TRUE.
        IF (CC2 .AND. IFLAGS(38).EQ.0) TERM3=.FALSE.
        TERM4 = .TRUE.
        TERM5 = .TRUE.
        IF (CC2) TERM5 = .FALSE. 
        TERM6 = .TRUE.

        DO_TAU = .TRUE.
        IF (HBAR_4LCCSD) DO_TAU=.FALSE.

        CALL FORMW4(ICORE(i0),MAXCOR,IUHF,TERM1,TERM2,TERM3,
     &              TERM4,TERM5,TERM6,DO_TAU,.TRUE.)

CJDW 2/6/95. Problems on Sun : TERM2 in FORMW5 can be changed, i.e. there
C            is an attempt to modify the constant ".TRUE.".
CJDW 3/29/96. Add CCSD=.TRUE. argument to FORMW4, FORMW5.
C
        TERM1=.TRUE.
        TERM2=.TRUE.
        TERM3=.TRUE.
        IF (CC2 .AND. IFLAGS(38).EQ.0) TERM3=.FALSE.
        TERM4=.TRUE.
        TERM5=.TRUE.
        TERM6=.TRUE.
        CALL FORMW5(ICORE(i0),MAXCOR,IUHF,TERM1,TERM2,TERM3,
     &              TERM4,TERM5,TERM6,DO_TAU,.TRUE.)

c************************************************************
c from now on lists 33-39 (newly created taus) are on 233-239
c*************************************************************

       ENDIF
C
C       IF(CCSD .OR. CC2) CALL FORMW1(ICORE(i0),MAXCOR,IUHF,.TRUE.)
C
       IF(IHBAR)THEN
C
C CALCULATE REMAINING HBAR CONTRIBUTIONS
C IFLAGS2(123)=2=HBARABCI=ON, IFLAGS2(122)=2=HBARABCD=ON
C IFLAGS(87)=3=EXCITE=EOMEE,METHOD is equivalence to 
C CALC (Iflags(2)=14,18,33,34 CCSDT-1b, CCSDT, CC3,CCSDT-T1T2
C (I doubt that anyone has checked any of those obscure 
C  methods). 
C
         IF (IFLAGS2(123) .EQ. 2) THEN
           CALL W5T1ABCD(ICORE(i0),MAXCOR,IUHF)
         ENDIF
C
C This block is for iterative (CCSDT-1b, CCSDT, CC3,CCSDT-T1T2??)
C methods, Not sure how many of them actually works.
C
         IF((METHOD.GE.14.AND.METHOD.LE.18).OR.METHOD.EQ.33.OR.
     &                                         METHOD.EQ.34    )THEN
           CALL W45T3DRV(ICORE(i0),MAXCOR,IUHF,106,126,13,4,.FALSE.)
         ENDIF

         IF (IFLAGS2(122) .EQ. 2) THEN
           CALL   FORMW1(ICORE(i0),MAXCOR,IUHF,.TRUE.)
           CALL HBARABCD(ICORE(i0),MAXCOR,IUHF,.TRUE.)
         ENDIF
C
        CALL HBRIJKA0(ICORE(i0),MAXCOR,IUHF)
        IF (.NOT. HBAR_4LCCSD) CALL FIXFBAR(ICORE(i0),MAXCOR,IUHF)
        CALL MODAIBC(ICORE(i0),MAXCOR,IUHF,ONE)

#ifdef _DEBUG_LVL0
CSSS        Call checkhbar(Icore(i0),Maxcor,iuhf)
#endif 

        write(*,5000)
 5000   FORMAT('  Calculation of Hbar completed.')
C
C-----------------------------------------------------------------------
C     Close triples file.
C-----------------------------------------------------------------------
        IF(T3STOR)THEN
          CALL AUXIOO(.FALSE.)
        ENDIF
C-----------------------------------------------------------------------
        WRITE(6,1010)
        IF (IFLAGS(35).NE.0) THEN
           CALL ACES_AUXCACHE_FLUSH
           CALL ACES_AUXCACHE_RESET
        END IF
C 10/2021
        IF (pCCD) THEN
           IFLAGS(h_IFLAGS_calc) = 8
           CALL PUTREC(20,"JOBARC","PCCD_RUN",1,1)
           CALL PUTREC(20,"JOBARC","IFLAGS  ",100,IFLAGS)
        ENDIF
        call aces_fin
        STOP
       ENDIF
C
C The above ENDIF is for IHBAR. If only IHBAR is set the calculation
C ends here (no lambda iterations).
C
C**********************************************************************
#ifdef _DEBUG_LVL0
        write(6,*) "Hbar before CC loop"
CSSS        Call checkhbar(Icore(i0),Maxcor,iuhf)
#endif 
      IF (LCC) THEN
         Write(6,"(a)") " Lambda is constructed for linear CC methods."
         DO ISPIN = 3, 3-2*IUHF, -1
            CALL MVSYM(ICORE(I0),MAXCOR,IUHF,1,43+ISPIN,143+ISPIN)
         ENDDO
         DO ISPIN=1,IUHF+1
             CALL UPDMOI(1,NT(ISPIN),ISPIN,190,0,0)
         END DO

         CALL GETLST(ICORE(I0), 1, 1, 1, 1, 90)
         CALL PUTLST(ICORE(I0), 1, 1, 1, 1, 190)

         IF (IUHF .NE. 0) THEN
             CALL GETLST(ICORE(I0), 1, 1, 1, 2, 90)
             CALL PUTLST(ICORE(I0), 1, 1, 1, 2, 190)
         ENDIF
         CALL LGENINT(ICORE(i0),MAXCOR,IUHF)
         CALL ACES_FIN
         STOP
      ENDIF 
C
C START OF CC LOOP
C
C  INITIALIZE THE RLE CONVERGENCE ACCELERATION
C
       RLECYC=0
       IF(MOD(IFLAGS(21),2).EQ.0)THEN
        CALL DRRLE(ICORE(i0),MAXCOR,IUHF,RLECYC,.TRUE.)
       ENDIF
       IF(IFLAGS(21).EQ.1)THEN
        CALL DIISLST(1,IUHF,METHOD.GE.6.AND.METHOD.NE.8)
       ENDIF

#ifdef _TTRANISL
C
C Caution: The following block is there to T^(t) = L games
C that Rod once in a while like to play.
C
      Write(6,*) 
      Write(*,"(a)") " Caution: The T^(t) = L approximation is invoked"
      IF (IFLAGS(19) .EQ. 0) THEN
         Write(*,"(a)") " Only respone density is constructed for L=T"
         Write(*,"(a)") " Use the density=response option"
         Call aces_exit(1)
      ENDIF 

      DO ISPIN = 3, 3-2*IUHF, -1
         CALL MVSYM(ICORE(I0),MAXCOR,IUHF,1,43+ISPIN,143+ISPIN)
      ENDDO
      DO ISPIN=1,IUHF+1
          CALL UPDMOI(1,NT(ISPIN),ISPIN,190,0,0)
      END DO

      CALL GETLST(ICORE(I0), 1, 1, 1, 1, 90)
      CALL PUTLST(ICORE(I0), 1, 1, 1, 1, 190)

      IF (IUHF .NE. 0) THEN
          CALL GETLST(ICORE(I0), 1, 1, 1, 2, 90)
          CALL PUTLST(ICORE(I0), 1, 1, 1, 2, 190)
      ENDIF
      CALL LGENINT(ICORE(i0),MAXCOR,IUHF)
      Write(*,"(a)") " Caution: The true lambada is not constructed"
      CALL ACES_FIN
      STOP
#endif
       DO 55 ICYCLE=1,NCYCLE
        RLECYC=RLECYC+1
C
C FILL LAMBDA(2) INCREMENTS WITH <IJ||AB> INTEGRALS
C
C AND ZERO THE LAMBDA(1) INCREMENTS (ONLY CCSD AND QCISD)
C
        CALL INITIN(ICORE(i0),MAXCOR,IUHF)
        IF(CCSD.OR.QCISD .OR.CC2) CALL INITSN(ICORE(i0),MAXCOR,IUHF)
c& new lines
C
C     CALCULATE OS CONTRIBUTION WITH THE H(EFF) DERIVATIVE
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals before multiplication"
      Call check_leom(Icore(i0),Maxcor,Iuhf)
#endif 
      IF(OSFLAG) THEN
         CALL OSL1(ICORE(i0),MAXCOR,IUHF)
         CALL OSL2(ICORE(i0),MAXCOR,IUHF)
      ENDIF
C& new line ends
C
C GENERATE V AND G INTERMEDIATES.
C
        IF (.NOT. CC2) THEN
           CALL LGENINT(ICORE(i0),MAXCOR,IUHF)
        ENDIF 
C
C COMPUTE F AND G INTERMEDIATE CONTRIBUTION TO LAMBDA(2) INCREMENT.
C
        CALL F1INL2(ICORE(i0),MAXCOR,IUHF)
        CALL F2INL2(ICORE(i0),MAXCOR,IUHF)

CSSS#ifdef _DCC_FLAG  
C Compute the contributions from F(ae) and F(me) stored at 91, 92 (9,10)
C location and add that to L2. That will make the overall scalle of 
C F(ae) and F(me) to be the correct scalling based on paralpha and parbeta.
C
        If (Ispar) Then
           Call Pdcc_f1inl2(Icore(i0),Maxcor,Iuhf,1.0D0)
           Call Pdcc_f2inl2(Icore(i0),Maxcor,Iuhf,1.0D0)
        Endif 
cSSS#endif 

#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals before G"
      Call check_leom(Icore(i0),Maxcor,Iuhf)
#endif
        CALL G1INL2(ICORE(i0),MAXCOR,IUHF)
        CALL G2INL2(ICORE(i0),MAXCOR,IUHF)

#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after G"
      Call check_leom(Icore(i0),Maxcor,Iuhf)
#endif 
C
        IF(CCSD.OR.QCISD)THEN
C
C L1 CONTRIBUTION TO L2
C
         CALL L1INL2A(ICORE(i0),MAXCOR,IUHF,.FALSE.)
         CALL L1INL2B(ICORE(i0),MAXCOR,IUHF,.FALSE.)

#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after L1inL2"
      Call check_leom(Icore(i0),Maxcor,Iuhf)
#endif
         If (Ispar) CALL RESTORE_CC_WMBEJ(ICORE(I0),MAXCOR,IUHF)
C
         CALL L1INL1(ICORE(i0),MAXCOR,IUHF,1)
         IF(IUHF.NE.0)CALL L1INL1(ICORE(i0),MAXCOR,IUHF,2)
C
         If (Ispar) CALL RESTORE_PDCC_WMBEJ(ICORE(I0),MAXCOR,IUHF)

#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after L1inL1"
      Call check_leom(Icore(i0),Maxcor,Iuhf)
#endif
C
C  FOR CCSD CALCULATE IN ADDITION
C
         IF(CCSD .AND. .NOT. CC2) THEN

CSSS#ifdef _DCC_FLAG
      If (Ispar) Then
C For PCC we need to reconstruct the standard G intermediates before
C computing the contributions to L1.
         
         IRREPX      = 1
         IHHA_LENGTH = IRPDPD(IRREPX,21)
         IHHB_LENGTH = IRPDPD(IRREPX,22)
         IPPA_LENGTH = IRPDPD(IRREPX,19)
         IPPB_LENGTH = IRPDPD(IRREPX,20)

         Call zerlst(ICORE(I0),IHHA_LENGTH,1,1,1,191)
         Call zerlst(ICORE(I0),IPPA_LENGTH,1,1,1,192)
         If (iuhf .ne. 0) then
            call zerlst(ICORE(I0),IHHB_LENGTH,1,1,2,191)
            call zerlst(ICORE(I0),IPPB_LENGTH,1,1,2,192)
         Endif

         Call Formg1(Icore(I0),Maxcor,Iuhf,1.0D0)
         Call Formg2(Icore(I0),Maxcor,Iuhf,1.0D0)
      Endif 
CSSS#endif
C
C 1)  G*T1 CONTRIBUTIONS TO L1
C
          CALL GT1INL1(ICORE(i0),MAXCOR,IUHF)
C
C 2)  G CONTRIBUTIONS TO L1
C
          CALL GWINL1(ICORE(i0),MAXCOR,IUHF)

         ENDIF
        ENDIF

#ifdef _DEBUG_LVLM
      Write(6,*) "The lambda residuals after GinL1"
      Call check_leom(Icore(i0),Maxcor,Iuhf)
      Call checkintms(Icore(i0),Maxcor,Iuhf)
#endif
C
C DO W INTERMEDIATE CONTRIBUTION TO T2 EQUATION (SAME AS THIRD ORDER CODE).
C
         CALL L2INL2(ICORE(i0),MAXCOR,IUHF,0)
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after l2lad"
      Call check_leom(Icore(i0),Maxcor,Iuhf)
#endif
         IF(CCSD .AND. .NOT. CC2 )
     &      CALL L1W1(ICORE(i0),MAXCOR,IUHF,1,143,90,26,60)
C
C COMPUTE THE CONTRIBUTION OF DOUBLES TO T1 AND
C DENOMINATOR WEIGHT THE NEW INCREMENTS 
C (ONLY CCSD AND QCISD METHODS)
C
        IF(CCSD.OR.QCISD.OR.CC2) THEN
          CALL L2INL1(ICORE(i0),MAXCOR,IUHF)
        ENDIF

        IF (pCCD .or. PCCDS .or. pCCDTS) THEN
C Ajith Perera, 10/2021.
        Write(6,"(20x,a)") " ----------Warning------------"
        Write(6,"(3a)")  "  A pCCD or pCCD like calculation is being",
     +                   " perforemed and the off-diagonal"
        write(6,"(a)")   " blocks of T2 is set to zero"
         Write(6,"(20x,a)") " -----------------------------"
         call pccd_reset_vcc(Icore(i0),Maxcor,Iuhf,63)
        ENDIF 
#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after l2inl1"
      Call check_leom(Icore(i0),Maxcor,Iuhf)
#endif
C
C DENOMINATOR WEIGHT T2 INCREMENTS TO FORM NEW T2.
C
        CALL NEWT2(ICORE(i0),MAXCOR,IUHF)

#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after newT2"
      Call check_leom(Icore(I0),Maxcor,Iuhf)
#endif

c&new lines
C
C  REMOVE L(Nm,Mn)
C
        IF(OSFLAG) THEN
           LISTL2=146
           IRREP=DIRPRD(IRREPOA,IRREPOB)
           NUMSYL=FSDPDAA(IRREP,ISYTYP(2,LISTL2))
           CALL FSPUT(ZILCH,1,NUMSYL,1,IRREP,LISTL2,'AAAA')
        ENDIF
c&new end
C
C WRITE OUT RMS AND MAX DIFFERENCES FOR NEW T2 AND PUT NEW ERROR VECTOR
C
cpr
#ifdef _TAILORED
      act_spc_cc = .true.
      act = "READ"
      SING_4TCC = .TRUE.
      SING_4TCC = .NOT.(IFLAGS(h_IFLAGS_calc) .EQ. 8 .OR.
     &                  IFLAGS(h_IFLAGS_calc) .EQ. 5 .OR.
     &                  pCCD .OR. pCCDTS .OR. pCCDTSD)
      Write(6,"(2a)")"First update of L1/L2 Res. lists 100+90(3,4)",
     &               " 61,62,63"
      if (act_spc_cc .and. (act .eq."READ" .or. act .eq."ZERO"))
     &   call modf_amps(ICORE(I0),MAXCOR,IUHF,0,SING_4TCC,'L',
     &                  act,"RESD")
#endif
        CALL DRTSTS(ICORE(i0),MAXCOR,ICYCLE,IUHF,ICONVG,ICONTL,
     &              SING,100,'L')
        CALL CMPENG2(ICORE(i0),MAXCOR,60,2,ECORR,ENERGY(ICYCLE+1,1),
     &               ENERGY(ICYCLE+1,2),IUHF,IPRTENG)
C
C When no convergence accelerations are used, we should be 
C able to proceed with the simple standard "JACOBI" procedure.
C This option was not included in the orginal implementation.
C In order to have this option, in addition to the following
C two calls, a small change in DIJAB.F in intprc.f was required.
C The change in DIJAB.F is to get the reciprocal of (D(i) - D(a))
C written into the list (9, 64) and the following calls update the 
C 190 and 143 lists. The list (9, 64) is used in l2inl1.f. 
C Ajith Perera 12/2003. 
C
        IF (ICONVG .NE. 0 .AND. IFLAGS(21).EQ. 3) THEN
           CALL DRMOVE(ICORE(i0),MAXCOR,IUHF,100,SING)
           CALL RNABIJ(ICORE(i0),MAXCOR,IUHF,'L')
        ENDIF 

        IF(IFLAGS(21).EQ.1)THEN
         CALL DODIIS0(ICORE(i0),MAXCOR/IINTFP,IUHF,1,ICYCLE,
     &                ICONVG,ICONTL,SING,144,61,190,0,90,2,70,
     &                '     ')
         CALL DRMOVE(ICORE(i0),MAXCOR,IUHF,100,SING)
          if(bRedundant)
     &      CALL RNABIJ(ICORE(i0),MAXCOR,IUHF,'L')
         IF(ICONVG.EQ.0)THEN
          CALL CMPENG2(ICORE(i0),MAXCOR,60,2,ECORR,
     &                 ENERGY(ICYCLE+1,1),
     &                 ENERGY(ICYCLE+1,2),IUHF,IPRTENG)

C#ifdef _DEBUG_LVLM
C Write the Lambada amplitudes to the disk.  This is primarily
C for Jesse Lutz work on TD-CCSD, Ajith Perera, 10/2015
C
         IF (IFLAGS2(171) .EQ. 1) CALL DDMPLGSS(ICORE(I0), MAXCOR, 
     &                                IUHF, 100, "LGUESS  ")
C#endif
C
C This call was moved out of IFLAGS(21) logic. Only needs to
C know whether the Lambda is converged, how does not matter!
C 04/2018.
CSSS        If (Ispar) Call Post_dccl_mods(Icore(I0),Maxcor,Iuhf)

         ENDIF
        ENDIF

        If (Ispar .and. iconvg .eq. 0) Call 
     &     Post_dccl_mods(Icore(I0),Maxcor,Iuhf)
#ifdef _DEBUG_LVLM
CSSS      Write(6,*)
CSSS      Write(6,"(a)") " Antisymmetry check of T2 and L2 Amps."
CSSS      Call Modf_T2L2_inL(ICore(I0),Maxcor,Iuhf)
#endif
C
        IF(ICONVG.EQ.0) THEN

C IHBAR2 logic. There is converged L1 and L2 are in the list 143.
C IFLAGS2(123)=2=HBARABCI=ON, IFLAGS2(122)=2=HBARABCD=ON
C IFLAGS(87)=3=EXCITE=EOMEE,IFLAGS(2)=22 is CCSD(T)
C IFLAGS2(124)=NT3EOMEE (controls EOM triples logic).

         IF(IHBAR2)THEN
           IF (IFLAGS2(123) .EQ. 2) THEN
             CALL W5T1ABCD(ICORE(i0),MAXCOR,IUHF)
           ENDIF
           IF( IFLAGS(2).EQ.22.AND.IFLAGS(87).EQ.3 .AND.
     &         (IFLAGS2(124).EQ.4.OR.IFLAGS2(124).EQ.6 .OR.
     &          IFLAGS2(124).EQ.8.OR.IFLAGS2(124).EQ.10    ) )THEN
             CALL W45T3DRV(ICORE(i0),MAXCOR,IUHF,384,388,13,0,.TRUE.)
           ENDIF
           IF (IFLAGS2(122) .EQ. 2) THEN
             CALL   FORMW1(ICORE(i0),MAXCOR,IUHF,.TRUE.)
             CALL HBARABCD(ICORE(i0),MAXCOR,IUHF, .TRUE.)
           ENDIF
           CALL HBRIJKA0(ICORE(i0),MAXCOR,IUHF) 
           CALL FIXFBAR(ICORE(i0),MAXCOR,IUHF)
           CALL MODAIBC(ICORE(i0),MAXCOR,IUHF,ONE)

         ELSE IF(IHBAR3) THEN

          CALL MODF(ICORE(i0),MAXCOR,IUHF,1,1,90,0,93,0,91,0,92,0)
C
c&new lines
         ELSE IF(OSFLAG) THEN
C
C     MODIFY L2(Nm,Mn)

          CALL OSL1(ICORE(i0),MAXCOR,IUHF)
c&new end
         ENDIF
#ifdef _TAILORED
       call getrec(0,'JOBARC','NDROPGEO',Ilength,Ijunk)
       If (Ilength .Gt. 0) act_spc_cc = .true.
       act = "SAVE"
       if (act_spc_cc .and. act .eq. "SAVE")
     &    call modf_amps(ICORE(I0),MAXCOR,IUHF,0,SING,'L',ACT,"AMPS")
#endif
C-----------------------------------------------------------------------
C     Close triples file.
C-----------------------------------------------------------------------
        IF(T3STOR)THEN
          CALL AUXIOO(.FALSE.)
        ENDIF
C-----------------------------------------------------------------------
         IF (IFLAGS(35).NE.0) THEN
            CALL ACES_AUXCACHE_FLUSH
            CALL ACES_AUXCACHE_RESET
         END IF
         CALL FINISH(ICYCLE+1,pCCD)
#ifdef _TAILORED
       call getrec(0,'JOBARC','NDROPGEO',Ilength,Ijunk)
       If (Ilength .Gt. 0) act_spc_cc = .true.
       act = "SAVE"
       if (act_spc_cc .and. act .eq. "SAVE")
     &    call modf_amps(ICORE(I0),MAXCOR,IUHF,0,SING,'L',ACT,"AMPS")
#endif
        ENDIF
C
C DO THE RLE EXTRAPOLATION.
C
        IF(MOD(IFLAGS(21),2).EQ.0)THEN
         CALL DRRLE(ICORE(i0),MAXCOR,IUHF,RLECYC,.TRUE.)
         CALL DRMOVE(ICORE(i0),MAXCOR,IUHF,100,SING)
         if(bRedundant)
     &     CALL RNABIJ(ICORE(i0),MAXCOR,IUHF,'L')
        ENDIF
C
C PRINT AMPLITUDES EVERY IEVERY CYCLES.
C
       IF(IFLAGS(1).GE.10)THEN
        IF(MOD(ICYCLE,IEVERY).EQ.0) THEN
         CALL AMPSUM(ICORE(i0),MAXCOR,IUHF,100,SING,'L')
        ENDIF
       ENDIF
55     CONTINUE
       WRITE(6,*)'  CC-Lambda did not converge !!! '
C
       CALL AMPSUM(ICORE(i0),MAXCOR,IUHF,100,SING,'L')
       call aces_exit(1)
      ENDIF

      IF (IFLAGS(35).NE.0) THEN
         CALL ACES_AUXCACHE_FLUSH
         CALL ACES_AUXCACHE_RESET
      END IF

      call aces_fin
      WRITE(6,1010)
      STOP
 1000 FORMAT(/,77('-'),/,T30,' Entering xlambda',/,77('-'),/)
 1010 FORMAT(/,77('-'),/,T30,' Exiting xlambda',/,77('-'))
      END
