      SUBROUTINE DWMBEJ(ICORE,MAXCOR,SPCASE,IUHF,C1,C2,C3,C4)
C
C Subroutine has been modified to accomodate PCCSD modifications (ring diagrams, beta).
c
C THIS ROUTINE AND DEPENDENTS COMPUTE THE RING-TYPE W(MBEJ)
C  INTERMEDIATE FOR ALL SIX POSSIBLE SPIN CASES.  THE ALGORITHM
C  ASSUMES IN-CORE STORAGE OF SYMMETRY PACKED TARGET, T2 AND W
C  VECTORS AND USES THE DPD SYMMETRY APPROACH TO EVALUATE THE
C  CONTRACTIONS.
C
C SPIN ORBITAL EQUATION FOR THIS INTERMEDIATE (EVENTUAL MODIFICATION
C                     TO INCLUDE T1 IS NECESSARY FOR CCSD):
C
C         W(mbej)=<mb||ej>-(1/2)*SUM T(jn,fb) <mn||ef>
C                                n,f
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION C1,C2,C3,C4,X,TWO,ONE,ONEM,HALF,ZILCH,TWOM
      CHARACTER*4 SPCASE,TAUTYP1,TAUTYP2
      DIMENSION ICORE(MAXCOR),IOFFT1(8,2)
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC,HBAR_4LCCSD,ADC2
      LOGICAL IJAB_DEBUG
      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC
      COMMON /INFO/ NOCA,NOCB,NVRTA,NVRTB
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT1AA,NT1BB,NF1AA,NF1BB,NF2AA,
     &             NF2BB
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
#include "parcc.com"
c
      DATA ZILCH /0.0/
      DATA ONE /1.0/
      DATA HALF/0.5/
      DATA TWO /2.0/
      DATA TWOM /-2.0/
      DATA ONEM/-1.0/
      NNP1O2(I)=I*(I+1)/2
C
      X=1.0
      ADC2 = .FALSE.
      ADC2 =  (IFLAGS(2) .EQ. 1 .AND. IFLAGS(87) .EQ. 3 .AND.
     &         IFLAGS2(117) .EQ. 10)
      HBAR_4LCCSD =  .FALSE.
      IF(C3.NE.0.0)CALL GETT1(ICORE,MAXCOR,MXCOR,IUHF,IOFFT1)
      IF (IFLAGS(2).EQ.6 .OR. IFLAGS2(117).EQ.7) HBAR_4LCCSD = .TRUE.
      IF (ADC2) THEN
         HBAR_4LCCSD = .TRUE. 
         C4 = 0.0D0
      ENDIF 

#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,"(a)") " The DCC parameters read from parcc file"
      Write(6,"(a,4a)") " The spin-case :", SPCASE
      Write(*,"(a,2L4)") " The ISPAR and COULOMB", ISPAR, COULOMB
      Write(*,"(a,F6.3)") " The Fae scale : ", Fae_scale 
      Write(*,"(a,F6.3)") " The Fmi scale: " , Fmi_scale 
      Write(*,"(a,F6.3)") " The Wmnij_scale: ", Wmnij_scale
      Write(*,"(a,F6.3)") " The Wmbej_scale: ", Wmbej_scale
      Write(6,"(a)") " The parcc pCC path is used" 
      Write(6,*)
      call check_ints(icore,maxcor,iuhf)
#endif
C
C********************************************************************
C
C DO SPIN CASE AAAA AND BBBB:
C
C       W(MBEJ)=<MB||EJ>-(1/2)*SUM [T(JN,FB)<MN||EF>-T(Jn,fB)<Mn|Ef>]
C
C       W(mbej)=<mb||ej>-(1/2)*SUM [T(jn,fb)<mn||ef>-T(jN,Fb)<mN|eF>]
C
C TERMS A AND B BOTH HANDLED BELOW.
C
      IF(SPCASE.EQ.'AAAA'.OR.SPCASE.EQ.'BBBB')THEN
       IF(SPCASE.EQ.'AAAA')THEN
        ISPIN=1
        TAUTYP1='AAAA'
       ELSEIF(SPCASE.EQ.'BBBB')THEN
        ISPIN=2
        TAUTYP1='BBBB'
       ENDIF
       MAXSIZ=0
       LSTINT=22+ISPIN
       LISTQ =53+ISPIN
       LISTWA=18+ISPIN
       IF (ISPAR .AND. COULOMB)  LISTWA = 118 + ISPIN
       LISTTA=33+ISPIN
       LISTWB=19-ISPIN  
       LISTTB=35+ISPIN
       IF(IUHF.EQ.0)LISTTB=37
       DO 100 IRREP=1,NIRREP 
        DSSYWA=IRPDPD(IRREP,ISYTYP(1,LISTWA))
        DSSYTA=IRPDPD(IRREP,ISYTYP(1,LISTTA))
        NMSYWA=IRPDPD(IRREP,ISYTYP(2,LISTWA))
        NMSYTA=IRPDPD(IRREP,ISYTYP(2,LISTTA))
        DSSYWB=IRPDPD(IRREP,ISYTYP(1,LISTWB))
        DSSYTB=IRPDPD(IRREP,ISYTYP(1,LISTTB))
        NMSYWB=IRPDPD(IRREP,ISYTYP(2,LISTWB))
        NMSYTB=IRPDPD(IRREP,ISYTYP(2,LISTTB))
        DISSYQ=IRPDPD(IRREP,ISYTYP(1,LISTQ))
        NUMSYQ=IRPDPD(IRREP,ISYTYP(2,LISTQ))
        I000=1
        I010=I000+IINTFP*DISSYQ*NUMSYQ
        I020=I010+IINTFP*NMSYWA*DSSYWA
        I030=I020+IINTFP*NMSYTA*DSSYTA
        I030=MAX(I030,I010+DISSYQ*NUMSYQ*IINTFP)
cjp    
        if(i030+DISSYQ*NUMSYQ*IINTFP .ge. MAXCOR) then
        write(6,*)'@DWMBEJ-E insufficient memory', 
     &      i030+DISSYQ*NUMSYQ*IINTFP, MAXCOR
        stop
        endif
        CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
        NSIZE=DISSYQ*NUMSYQ
        IF(C1.NE.0.0)THEN
          CALL GETLST(ICORE(I000),1,NUMSYQ,1,IRREP,LSTINT) 
         IF(UCC) THEN
          CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)
          CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
         ENDIF  
        ENDIF
        IF(C4.NE.0.0)THEN
         CALL GETLST(ICORE(I020),1,NMSYTA,1,IRREP,LISTTA)

         if(ispar .and. coulomb) then 
CSSS#ifdef _DCC_FLAG
           call dscal(nmsyta*dssyta,Wmbej_scale,icore(i020),1)
CSSS#else
CSSS           call dscal(nmsyta*dssyta,2.0D0*parbeta,icore(i020),1)
CSSS#endif 
         else if (ispar) then 
           call dscal(nmsyta*dssyta,Wmbej_scale,icore(i020),1)
         endif 

         IF(C3.NE.0.0 .AND. .NOT. HBAR_4LCCSD)THEN
          CALL GETARGS(LISTTA,I1,I2)
          CALL F2TAU(ICORE(I020),ICORE(IOFFT1(1,I1)),
     &                ICORE(IOFFT1(1,I2)),DSSYTA,NMSYTA,
     &                POP(1,I1),POP(1,I2),VRT(1,I1),VRT(1,I2),
     &                IRREP,TWO,TAUTYP1)
         ENDIF

#ifdef _DEBUG_LVLM
      if (ispar .and. coulomb) then
      Write(*,"(a)") "The Coulomb only <MN|EF> Tau(IJ,EF) is computed"
      Write(*,"(a,I4,1x,a)") "The Coulomb only lists",listwa,"is used"
      write(*,"(a,1x,a,1x,a)") "The spin-case", "ABAB", "BABA"
      endif 
#endif 
         CALL CNTRCT(ICORE(I000),ICORE(I010),ICORE(I020),
     &               DISSYQ,NUMSYQ,DSSYTA,NMSYTA,DSSYWA,NMSYWA,
     &               IRREP,LISTWA,LISTWA,LISTTA,'WxT',X,-C4,0)
         IF(UCC) THEN
          CALL DT(TWO,ICORE(I000),DISSYQ,NUMSYQ,'ADD',icore(i030))
          CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)
          CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
         ENDIF  
         I020=I010+IINTFP*NMSYWB*DSSYWB
         I030=I020+IINTFP*NMSYTB*DSSYTB
         CALL GETLST(ICORE(I020),1,NMSYTB,1,IRREP,LISTTB)

         if(ispar .and. coulomb) then

CSSS#ifdef _DCC_FLAG
           call dscal(nmsyta*dssyta,Wmbej_scale,icore(i020),1)
CSSS#else
CSSS           call dscal(nmsyta*dssyta,2.0D0*parbeta,icore(i020),1)
CSsS#endif 
         else if (ispar) then 
           call dscal(nmsyta*dssyta,Wmbej_scale,icore(i020),1)
         endif

         CALL CNTRCT(ICORE(I000),ICORE(I010),ICORE(I020),
     &               DISSYQ,NUMSYQ,DSSYTB,NMSYTB,DSSYWB,NMSYWB,
     &               IRREP,LISTWB,LISTWB,LISTTB,'WxT',X,C4,0)
         IF(UCC) THEN
          CALL DT(TWO,ICORE(I000),DISSYQ,NUMSYQ,'ADD',icore(i030))
         ENDIF  
        ENDIF
        CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)
#ifdef _DEBUG_LVL0
        Write(6,"(a,i4)") "AAAA and BBBB, The target list",listq
        call checksum(SPCASE, ICORE(I000), NSIZE)
#endif
100    CONTINUE
C
C DO SPIN CASES ABAB AND BABA:
C
C       W(MbEj)=<Mb|Ej>-(1/2)*SUM [T(jn,fb)<Mn|Ef>-T(jN,bF)<MN||EF>]
C
C       W(mBeJ)=<mB|eJ>-(1/2)*SUM [T(JN,FB)<mN|eF>-T(Jn,Bf)<mn||ef>]
C    

      ELSEIF(SPCASE.EQ.'ABAB'.OR.SPCASE.EQ.'BABA')THEN
#ifdef _DEBUG_LVLM
        IJAB_DEBUG=.FALSE.
        Write(6,"(a,a)") " @-DWBEJ: List 121 and 122 instead of",
     +                   " 17 and 18 is read for <MB||EJ> Ints"
#endif
       IF(IUHF.NE.0)THEN
       IF(SPCASE.EQ.'ABAB')THEN
        LSTINT=18
CSSS        IF (IJAB_DEBUG) LSTINT=122
        LISTQ=56
        LISTQ1=57
        IF(IUHF.EQ.0)THEN
         LISTTA=34
        ELSE
         LISTTA=35
        ENDIF
        LISTTB=37
        LISTWA=18
        LISTWB=19
        IF (ISPAR .AND. COULOMB)  LISTWB = 119
        TAUTYP1='BBBB'
       ELSEIF(SPCASE.EQ.'BABA')THEN 
        LSTINT=17
CSSS        IF (IJAB_DEBUG) LSTINT=121
        LISTQ=57
        LISTQ1=56
        LISTTA=34
        LISTWA=17
        LISTTB=36
        LISTWB=20
        IF (ISPAR .AND. COULOMB)  LISTWB = 120
        TAUTYP1='AAAA'
       ENDIF
       DO 200 IRREP=1,NIRREP 
        DSSYWA=IRPDPD(IRREP,ISYTYP(1,LISTWA))
        DSSYTA=IRPDPD(IRREP,ISYTYP(1,LISTTA))
        NMSYWA=IRPDPD(IRREP,ISYTYP(2,LISTWA))
        NMSYTA=IRPDPD(IRREP,ISYTYP(2,LISTTA))
        DSSYWB=IRPDPD(IRREP,ISYTYP(1,LISTWB))
        DSSYTB=IRPDPD(IRREP,ISYTYP(1,LISTTB))
        NMSYWB=IRPDPD(IRREP,ISYTYP(2,LISTWB))
        NMSYTB=IRPDPD(IRREP,ISYTYP(2,LISTTB))
        DISSYQ=IRPDPD(IRREP,ISYTYP(1,LISTQ))
        NUMSYQ=IRPDPD(IRREP,ISYTYP(2,LISTQ))
        I000=1
        I010=I000+IINTFP*DISSYQ*NUMSYQ
        I020=I010+IINTFP*NMSYWA*DSSYWA
        I030=I020+IINTFP*NMSYTA*DSSYTA
cjp    
        if(i030+DISSYQ*NUMSYQ*IINTFP .ge. MAXCOR) then
        write(6,*)'@DWMBEJ-E insufficient memory', 
     &      i030+DISSYQ*NUMSYQ*IINTFP, MAXCOR
        stop
        endif
        CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
        NSIZE=DISSYQ*NUMSYQ
        IF(C1.NE.0.0)THEN
          CALL GETLST(ICORE(I000),1,NUMSYQ,1,IRREP,LSTINT)  
CSSS          write(6,*) "Intentionally <mb|ej> is zero"
CSSS          CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
CSSS          call checksum("MBEJ",ICORE(I000),NUMSYQ*DISSYQ)
         IF(UCC) THEN
          CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)
          CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
         ENDIF  
        ENDIF
        IF(C4.NE.0.0)THEN
         CALL GETLST(ICORE(I020),1,NMSYTA,1,IRREP,LISTTA)

         if(ispar .and. coulomb) then
CSSS#ifdef _DCC_FLAG
           call dscal(nmsyta*dssyta,Wmbej_scale,icore(i020),1)
cSSS#else
CSSS           call dscal(nmsyta*dssyta,2.0D0*parbeta,icore(i020),1)
CSSS#endif 
         else if (ispar) then
           call dscal(nmsyta*dssyta,Wmbej_scale,icore(i020),1)
         endif

         IF(C3.NE.0.0 .AND. .NOT. HBAR_4LCCSD)THEN
          CALL GETARGS(LISTTA,I1,I2)
          CALL F2TAU(ICORE(I020),ICORE(IOFFT1(1,I1)),
     &                ICORE(IOFFT1(1,I2)),DSSYTA,NMSYTA,
     &                POP(1,I1),POP(1,I2),VRT(1,I1),VRT(1,I2),
     &                IRREP,TWO,TAUTYP1)
         ENDIF

         CALL CNTRCT(ICORE(I000),ICORE(I010),ICORE(I020),
     &               DISSYQ,NUMSYQ,DSSYTA,NMSYTA,DSSYWA,NMSYWA,
     &               IRREP,LISTWA,LISTWA,LISTTA,'WxT',X,C4,0) 
         IF(UCC) THEN
          CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)
          CALL DT(TWO,ICORE(I000),DISSYQ,NUMSYQ,'NSY',ICORE(I010))
          CALL SUMSYM2(ICORE(I010),ICORE(I000),NSIZE,1,IRREP,LISTQ1)
          CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
         ENDIF  
         I020=I010+IINTFP*NMSYWB*DSSYWB
         I030=I020+IINTFP*NMSYTB*DSSYTB
         I030=MAX(I030,I010+DISSYQ*NUMSYQ*IINTFP)
cjp    
        if(i030+DISSYQ*NUMSYQ*IINTFP .ge. MAXCOR) then
        write(6,*)'@DWMBEJ-E insufficient memory', 
     &      i030+DISSYQ*NUMSYQ*IINTFP, MAXCOR
        stop
        endif
         CALL GETLST(ICORE(I020),1,NMSYTB,1,IRREP,LISTTB)

         if(ispar .and. coulomb) then
CSSS#ifdef _DCC_FLAG
           call dscal(nmsyta*dssyta,Wmbej_scale,icore(i020),1)
CSSS#else
CSSSS           call dscal(nmsyta*dssyta,2.0D0*parbeta,icore(i020),1)
CSSS#endif 
         else if (ispar) then
           call dscal(nmsyta*dssyta,Wmbej_scale,icore(i020),1)
         endif

#ifdef _DEBUG_LVLM
      if (ispar .and. coulomb) then
      Write(*,"(a)") "The Coulomb only <MN|EF> Tau(IJ,EF) is computed"
      Write(*,"(a,I4,1x,a)") "The Coulomb only lists",listwb,"is used"
      write(*,"(a,1x,a,1x,a)") "The spin-case", "ABAB", "BABA"
      endif
#endif
         CALL CNTRCT(ICORE(I000),ICORE(I010),ICORE(I020),
     &               DISSYQ,NUMSYQ,DSSYTB,NMSYTB,DSSYWB,NMSYWB,
     &               IRREP,LISTWB,LISTWB,LISTTB,'WxT',X,-C4,0)
         IF(UCC) THEN
          CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)
          CALL DT(TWO,ICORE(I000),DISSYQ,NUMSYQ,'NSY',ICORE(I010))
          CALL SUMSYM2(ICORE(I010),ICORE(I000),NSIZE,1,IRREP,LISTQ1)
         ENDIF  
        ENDIF
        IF(.NOT.UCC)
     &  CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)
#ifdef _DEBUG_LVL0
        Write(6,"(a,i4)") "The ABAB and BABA, target list",listq
        call checksum(SPCASE, ICORE(I000), NSIZE)
#endif
200    CONTINUE
C
C SPIN ADAPTED CODE FOR RHF "ABAB" SPIN CASE
C
       ELSE

#ifdef _DEBUG_LVL0
       IJAB_DEBUG = .FALSE.
       write(*,"(a)") "Entering RHF only block"
       write(*,*)
#endif 
       LSTINT=18
CSSS       IF (IJAB_DEBUG)  LSTINT=122
       LSTINTB=25
CSSS       IF (IJAB_DEBUG)  LSTINTB=147
       LISTQ=56
       LISTTA=37
       LISTTB=39
       LISTWA=18
       LISTWB=21
CSSS       IF (IJAB_DEBUG) LISTWB=123
       TAUTYP1='ABAB'
       DO 201 IRREP=1,NIRREP 
        DSSYWA=IRPDPD(IRREP,ISYTYP(1,LISTWA))
        DSSYTA=IRPDPD(IRREP,ISYTYP(1,LISTTA))
        NMSYWA=IRPDPD(IRREP,ISYTYP(2,LISTWA))
        NMSYTA=IRPDPD(IRREP,ISYTYP(2,LISTTA))
        DSSYWB=IRPDPD(IRREP,ISYTYP(1,LISTWB))
        DSSYTB=IRPDPD(IRREP,ISYTYP(1,LISTTB))
        NMSYWB=IRPDPD(IRREP,ISYTYP(2,LISTWB))
        NMSYTB=IRPDPD(IRREP,ISYTYP(2,LISTTB))
        DISSYQ=IRPDPD(IRREP,ISYTYP(1,LISTQ))
        NUMSYQ=IRPDPD(IRREP,ISYTYP(2,LISTQ))
        I000=1
        I010=I000+IINTFP*DISSYQ*NUMSYQ
        I020=I010+IINTFP*NMSYWA*DSSYWA
        I030=I020+IINTFP*NMSYTA*DSSYTA
cjp    
        if(i030+DISSYQ*NUMSYQ*IINTFP .ge. MAXCOR) then
        write(6,*)'@DWMBEJ-E insufficient memory', 
     &      i030+DISSYQ*NUMSYQ*IINTFP, MAXCOR
        stop
        endif
        IF(C4.NE.0.0)THEN
#ifdef _DEBUG_LVLM
      if (ispar .and. coulomb) then
      Write(*,*) "The spin case ", "ABAB"
      Write(*,"(a)") "spin-adapated T and <MN|EF> is used"
      endif
#endif
C
C SPIN ADAPT THE T AMPLITUDES
C
         CALL GETLST(ICORE(I020),1,NMSYTA,1,IRREP,LISTTA)

         if(ispar .and. coulomb) then
CSSS#ifdef _DCC_FLAG
           call dscal(nmsyta*dssyta,Wmbej_scale,icore(i020),1)
CSSS#else
cSSS           call dscal(nmsyta*dssyta,2.0D0*parbeta,icore(i020),1)
CSSS#endif 
         else if (ispar) then 
           call dscal(nmsyta*dssyta,Wmbej_scale,icore(i020),1)
         endif

         CALL GETLST(ICORE(I010),1,NMSYTA,1,IRREP,LISTTB)

C         if(ispar .and. coulomb) then
C            call dscal(nmsyta*dssyta,2.0D0*parbeta,icore(i020),1)
C         else if (ispar) then 
C           call dscal(nmsyta*dssyta,parbeta,icore(i020),1)
C         endif

         if(ispar .and. coulomb) then

CSSS#ifdef _DCC_FLAG
           call dscal(nmsyta*dssyta,Wmbej_scale,icore(i010),1)
CSSS#else
CSSS            call dscal(nmsyta*dssyta,2.0D0*parbeta,icore(i010),1)
CSSS#endif 
         else if (ispar) then 
           call dscal(nmsyta*dssyta,Wmbej_scale,icore(i010),1)
         endif

         IF(C3.NE.0.0 .AND. .NOT. HBAR_4LCCSD)THEN
          CALL GETARGS(LISTTA,I1,I2)
          CALL F2TAU(ICORE(I010),ICORE(IOFFT1(1,I1)),
     &                ICORE(IOFFT1(1,I2)),DSSYTA,NMSYTA,
     &                POP(1,I1),POP(1,I2),VRT(1,I1),VRT(1,I2),
     &                IRREP,TWO,TAUTYP1)
         ENDIF

         CALL SSCAL (NMSYTA*NMSYTA,TWO,ICORE(I020),1)
         CALL SAXPY (NMSYTA*NMSYTA,ONEM,ICORE(I010),1,ICORE(I020),1)

        IF (ispar .and. coulomb) Then

         CALL CNTRCT_CONLY(ICORE(I000),ICORE(I010),ICORE(I020),
     &                     DISSYQ,NUMSYQ,DSSYTA,NMSYTA,DSSYWA,
     &                     NMSYWA,IRREP,LISTWA,LISTWB,LISTTA,
     &                     'WxT',ZILCH,-C4,1)

         ELSE

         CALL CNTRCT(ICORE(I000),ICORE(I010),ICORE(I020),
     &                DISSYQ,NUMSYQ,DSSYTA,NMSYTA,DSSYWA,NMSYWA,
     &                IRREP,LISTWA,LISTWB,LISTTA,'WxT',ZILCH,-C4,1)

         ENDIF 
C
         IF(UCC) THEN
          CALL DT(TWO,ICORE(I000),DISSYQ,NUMSYQ,'ADD',icore(i030))
         ENDIF  
        ENDIF
        NSIZE=DISSYQ*NUMSYQ
C
C NOW INCLUDE THE SPIN-ADAPTED INTEGRALS
C
        IF(C1.NE.0.0)THEN
C
C I don't think we need to do anything here for ISPAR. There is no 
C antisymmetry in these integrals. 

         CALL GETLST(ICORE(I010),1,NUMSYQ,1,IRREP,LSTINT)
         CALL GETLST(ICORE(I020),1,NUMSYQ,1,IRREP,LSTINTB)
         CALL SSCAL (DISSYQ*DISSYQ,TWO,ICORE(I010),1)
         CALL SAXPY (DISSYQ*DISSYQ,ONEM,ICORE(I020),1,ICORE(I010),1) 

CSSS         Write(6,*) "Intentionally zeroing for debuggin"
CSSS         call izero(ICORE(I010),DISSYQ*DISSYQ)

         IF(C4.NE.0.0)THEN
          CALL SAXPY (DISSYQ*DISSYQ,ONE,ICORE(I010),1,ICORE(I000),1)
         ELSE
c YAU : old
c         CALL ICOPY(DISSYQ*DISSYQ*IINTFP,ICORE(I010),1,ICORE(I000),1)
c YAU : new
          CALL DCOPY(DISSYQ*DISSYQ,ICORE(I010),1,ICORE(I000),1)
c YAU : end
         ENDIF
        ENDIF
C
C WE ALSO NEED SPIN-ADAPTED T1RING STUFF FOR CCSD
C
c       IF(C3.NE.0.0)THEN
CJDW 3/20/98. The T1RING stuff must also be picked up when we have
C             non-HF SDQ-MBPT(4) or MBPT(4). Also, require C4.NE.0
C             so that we know we are calculating second-order MBEJ,
C             which implies that T1RING term has been calculated
C             and must be included.

        IF(C3.NE.0.0 .OR. HBAR_4LCCSD .OR. (IFLAGS(38).NE.0 .AND. 
     &    (IFLAGS(2).EQ.3.OR.IFLAGS(2).EQ.4) .AND. C4.NE.0))THEN

#ifdef _DEBUG_LVLM
       write(*,"(a)") "Entering RHF only T1ring"
       write(*,*)
#endif       
         CALL GETLST(ICORE(I010),1,NUMSYQ,1,IRREP,56)
         CALL GETLST(ICORE(I020),1,NUMSYQ,1,IRREP,58)
         CALL SSCAL(NSIZE,TWO,ICORE(I010),1)
         CALL SAXPY(NSIZE,ONEM,ICORE(I020),1,ICORE(I010),1)
         CALL SAXPY(NSIZE,ONE ,ICORE(I010),1,ICORE(I000),1)
         CALL PUTLST(ICORE(I000),1,NUMSYQ,1,IRREP,56)
#ifdef _DEBUG_LVL0
       write(*,"(a,I4)") "Writing RHF only ABAB to disk: ",LISTQ
       call checksum("ABAB",ICORE(I000),NSIZE)
       write(*,*)
#endif       
        ELSE
#ifdef _DEBUG_LVL0
       write(*,"(a,I4)") "Writing RHF only ABAB to disk: ",LISTQ
       call checksum("ABAB",ICORE(I000),NSIZE)
       write(*,*)
#endif       
         CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)
        ENDIF
201    CONTINUE
       ENDIF
C
C
      ELSEIF(SPCASE.EQ.'ABBA'.OR.SPCASE.EQ.'BAAB')THEN
C
C DO SPIN CASES ABBA AND BAAB:
C
C       W(MbeJ) = -<Mb|Je>+(1/2)*SUM [T(Jn,Fb)<Mn|Fe>]
C
C       W(mBEj) = -<mB|jE>+(1/2)*SUM [T(jN,fB)<mN|fE>]
C
       IF(SPCASE.EQ.'ABBA')THEN
        LISTQ=58
        LSTINT=25
        LISTTA=38
        IF(IUHF.EQ.0)LISTTA=39
        LISTWA=21
        TAUTYP1='ABAB'
       ELSEIF(SPCASE.EQ.'BAAB')THEN
        LISTQ=59
        LSTINT=26
        LISTTA=39
        LISTWA=22
        TAUTYP1='BABA'
       ENDIF
       DO 300 IRREP=1,NIRREP 
        DSSYWA=IRPDPD(IRREP,ISYTYP(1,LISTWA))
        DSSYTA=IRPDPD(IRREP,ISYTYP(1,LISTTA))
        NMSYWA=IRPDPD(IRREP,ISYTYP(2,LISTWA))
        NMSYTA=IRPDPD(IRREP,ISYTYP(2,LISTTA))
        DISSYQ=IRPDPD(IRREP,ISYTYP(1,LISTQ))
        NUMSYQ=IRPDPD(IRREP,ISYTYP(2,LISTQ))
        I000=1
        I010=I000+IINTFP*DISSYQ*NUMSYQ
        I020=I010+IINTFP*NMSYWA*DSSYWA
        I030=I020+IINTFP*NMSYTA*DSSYTA
        I030=MAX(I030,I010+DISSYQ*NUMSYQ*IINTFP)
cjp    
        if(i030+DISSYQ*NUMSYQ*IINTFP .ge. MAXCOR) then
        write(6,*)'@DWMBEJ-E insufficient memory', 
     &      i030+DISSYQ*NUMSYQ*IINTFP, MAXCOR
        stop
        endif
        CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
        NSIZE=DISSYQ*NUMSYQ
C
C LIST ORDERING OF TARGET IS SAME AS INTEGRAL INVOLVED IN
C  CONTRACTION, SO JUST USE THIS SET OF POINTERS TO FETCH 
C  PACKED INTEGRAL LISTS.
C
        IF(C1.NE.0.0)THEN
          CALL GETLST(ICORE(I000),1,NUMSYQ,2,IRREP,LSTINT)
CSSS         Write(6,*) "Intentionally zeroing for debuggin"
CSSS         call izero(ICORE(I000),DISSYQ*DISSYQ)
         IF(UCC) THEN
          CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)
          CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
         ENDIF  
        ENDIF
        IF(C4.NE.0.0)THEN
         CALL GETLST(ICORE(I020),1,NMSYTA,1,IRREP,LISTTA)

         if(ispar .and. coulomb) then
CSSS#ifdef _DCC_FLAG
           call dscal(nmsyta*dssyta,Wmbej_scale,icore(i020),1)
CSSS#else
CSSS           call dscal(nmsyta*dssyta,2.0D0*parbeta,icore(i020),1)
CSSS#endif 
         else if (ispar) then
           call dscal(nmsyta*dssyta,Wmbej_scale,icore(i020),1)
         endif

         IF(C3.NE.0.0 .AND. .NOT. HBAR_4LCCSD)THEN
          CALL GETARGS(LISTTA,I1,I2)
          CALL F2TAU(ICORE(I020),ICORE(IOFFT1(1,I1)),
     &                ICORE(IOFFT1(1,I2)),DSSYTA,NMSYTA,
     &                POP(1,I1),POP(1,I2),VRT(1,I1),VRT(1,I2),
     &                IRREP,TWO,TAUTYP1)
         ENDIF

        IF (.NOT. (ispar .and. coulomb)) THEN

         CALL CNTRCT(ICORE(I000),ICORE(I010),ICORE(I020),
     &               DISSYQ,NUMSYQ,DSSYTA,NMSYTA,DSSYWA,NMSYWA,
     &               IRREP,LISTWA,LISTWA,LISTTA,'WxT',X,C4,0)

         ENDIF 

         IF(UCC) THEN
c&          CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)
c&          CALL DT(ONE,ICORE(I000),DISSYQ,NUMSYQ,'EXT',icore(i030))
          CALL DT(TWO,ICORE(I000),DISSYQ,NUMSYQ,'ADD',icore(i030))
         ENDIF  
        ENDIF

#ifdef _DEBUG_LVL0
       write(*,"(a,I4)") "The ABBA and BAAB, the target list:",LISTQ
       call checksum(SPCASE,ICORE(I000),NSIZE)
       write(*,*)
#endif
        CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)
300    CONTINUE
      ENDIF
      RETURN
      END
