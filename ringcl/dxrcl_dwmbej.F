      SUBROUTINE DXRCL_DWMBEJ(ICORE,MAXCOR,SPCASE,IUHF,C4)
C
C THIS ROUTINE AND DEPENDENTS COMPUTE THE RING-TYPE W(MBEJ)
C INTERMEDIATE FOR ALL SIX POSSIBLE SPIN CASES FOR RING-CC CALCULATIONS.
C
C         W(mbej)=<mb||ej>-(1/2)*SUM T(jn,fb) <mn||ef>
C                                n,f
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION C4,X,TWO,ONE,ONEM,HALF,ZILCH,TWOM
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      CHARACTER*4 SPCASE,TAUTYP1
      DIMENSION ICORE(MAXCOR),IOFFT1(8,2)

      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
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
c
      DATA ZILCH /0.0/
      DATA ONE /1.0/
      DATA HALF/0.5/
      DATA TWO /2.0/
      DATA TWOM /-2.0/
      DATA ONEM/-1.0/
C
      X=1.0D0

C DO SPIN CASE AAAA AND BBBB (Antisymmetric).
C
C W(MBEJ)=(-)<BM||EJ>+(1/2)*SUM [(-)(-)T(JN,FB)<MN||EF>+(-)T(Jn,fB)<Mn|Ef>]
C
C W(mbej)=(-)<bm||ej>+(1/2)*SUM [(-)(-)T(jn,fb)<mn||ef>+(-)T(jN,Fb)<mN|eF>]
C
C As above expressions indicates for antisymmetric cases the fisrt term is 
C with -C4 (C4 is negative) to account form (-) storage of T2 and the 
C second terms is +C4 to account for the negative sign. The Overall the 
C terms is computed with negtive sign. 

C For case where T2 is symmetric, but the integrals are antisymmetric
C the lead integral is stored as -<BM||EJ>. Therefore scale the integral 
C by negtive sign. The overall contribution is computed as positive 
c sign.
C
      IF(SPCASE.EQ.'AAAA'.OR.SPCASE.EQ.'BBBB')THEN
       IF(SPCASE.EQ.'AAAA')THEN
        ISPIN=1
       ELSEIF(SPCASE.EQ.'BBBB')THEN
        ISPIN=2
       ENDIF
       MAXSIZ=0
       LSTINT=22+ISPIN
       LISTQ =53+ISPIN
       LISTWA=18+ISPIN
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

        IF (i030+DISSYQ*NUMSYQ*IINTFP .ge. MAXCOR) then
        write(6,*)'@DWMBEJ-E insufficient memory', 
     &      i030+DISSYQ*NUMSYQ*IINTFP, MAXCOR
        CALL ERREX
        ENDIF 

        CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
        NSIZE=DISSYQ*NUMSYQ
        CALL GETLST(ICORE(I000),1,NUMSYQ,1,IRREP,LSTINT) 
        CALL DSCAL(NSIZE,ONEM,ICORE(I000),1)
        CALL GETLST(ICORE(I020),1,NMSYTA,1,IRREP,LISTTA)

#ifdef _DEBUG_LVL0
        Write(6,"(a,i4)") "AAAA and BBBB-T2 list",listta 
        call checksum(SPCASE, ICORE(I020), NMSYTA*DSSYTA)
        Write(6,"(a,i4)") "AAAA and BBBB-W2 list",lstint
        call checksum(SPCASE, ICORE(I000), NUMSYQ*DISSYQ)
#endif
        CALL CNTRCT(ICORE(I000),ICORE(I010),ICORE(I020),
     &              DISSYQ,NUMSYQ,DSSYTA,NMSYTA,DSSYWA,NMSYWA,
     &              IRREP,LISTWA,LISTWA,LISTTA,'WxT',X,-C4,0)
        I020=I010+IINTFP*NMSYWB*DSSYWB
        I030=I020+IINTFP*NMSYTB*DSSYTB
        CALL GETLST(ICORE(I020),1,NMSYTB,1,IRREP,LISTTB)
        CALL CNTRCT(ICORE(I000),ICORE(I010),ICORE(I020),
     &              DISSYQ,NUMSYQ,DSSYTB,NMSYTB,DSSYWB,NMSYWB,
     &              IRREP,LISTWB,LISTWB,LISTTB,'WxT',X,-C4,0)
#ifdef _DEBUG_LVL0
        Write(6,"(a,i4)") "AAAA and BBBB-the target list",listq
        call checksum(SPCASE, ICORE(I000), NSIZE)
#endif
        CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)

100    CONTINUE
C
C DO SPIN CASES ABAB AND BABA (original)
C
C  W(MbEj)=<Mb|Ej>+(1/2)*SUM [(-)T(jn,fb)<Mn|Ef>+T(jN,bF)<MN||EF>]
C
C  W(mBeJ)=<mB|eJ>+(1/2)*SUM [(-)T(JN,FB)<mN|eF>+T(Jn,Bf)<mn||ef>]

C As above expressions indicates for antisymmetric cases the second term is
C with C4 (C4 is negative) to account form (-) storage of T2 and the
C third term is with (-)C4 (positive). The overall term is positive. 
C For symmetric cases, all the terms are positve and -C4 is used.
C No changes to above comments for cases where integrals are antisymmetric
C but T2 is symmetric.
C
      ELSEIF(SPCASE.EQ.'ABAB'.OR.SPCASE.EQ.'BABA')THEN

       IF(IUHF.NE.0)THEN
       IF(SPCASE.EQ.'ABAB')THEN
        LSTINT=18
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
       ELSEIF(SPCASE.EQ.'BABA')THEN 
        LSTINT=17
        LISTQ=57
        LISTQ1=56
        LISTTA=34
        LISTWA=17
        LISTTB=36
        LISTWB=20
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

        IF (i030+DISSYQ*NUMSYQ*IINTFP .ge. MAXCOR) then
        write(6,*)'@DWMBEJ-E insufficient memory', 
     &      i030+DISSYQ*NUMSYQ*IINTFP, MAXCOR
        CALL ERREX
        ENDIF 

        CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
        NSIZE=DISSYQ*NUMSYQ

        CALL GETLST(ICORE(I000),1,NUMSYQ,1,IRREP,LSTINT)  
        CALL GETLST(ICORE(I020),1,NMSYTA,1,IRREP,LISTTA)
#ifdef _DEBUG_LVL0
        Write(6,"(a,i4)") "AAAA and BBBB-T2 list",listta 
        call checksum(SPCASE, ICORE(I020), NMSYTA*DSSYTA)
        Write(6,"(a,i4)") "ABAB-W2 list",lstint
        call checksum(SPCASE, ICORE(I000), DISSYQ*NUMSYQ)
#endif
        CALL CNTRCT(ICORE(I000),ICORE(I010),ICORE(I020),
     &              DISSYQ,NUMSYQ,DSSYTA,NMSYTA,DSSYWA,NMSYWA,
     &              IRREP,LISTWA,LISTWA,LISTTA,'WxT',X,-C4,0) 
        I020=I010+IINTFP*NMSYWB*DSSYWB
        I030=I020+IINTFP*NMSYTB*DSSYTB
        I030=MAX(I030,I010+DISSYQ*NUMSYQ*IINTFP)

        IF (i030+DISSYQ*NUMSYQ*IINTFP .ge. MAXCOR) then
        write(6,*)'@DWMBEJ-E insufficient memory', 
     &      i030+DISSYQ*NUMSYQ*IINTFP, MAXCOR
        CALL ERREX
        ENDIF 

        CALL GETLST(ICORE(I020),1,NMSYTB,1,IRREP,LISTTB)
        CALL CNTRCT(ICORE(I000),ICORE(I010),ICORE(I020),
     &              DISSYQ,NUMSYQ,DSSYTB,NMSYTB,DSSYWB,NMSYWB,
     &              IRREP,LISTWB,LISTWB,LISTTB,'WxT',X,-C4,0)

#ifdef _DEBUG_LVL0
        Write(6,"(a,i4)") "The ABAB and BABA-target list",listq
        call checksum(SPCASE, ICORE(I000), NSIZE)
#endif
        CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)

200    CONTINUE
C
C SPIN ADAPTED CODE FOR RHF "ABAB" SPIN CASE
C
       ELSE

#ifdef _DEBUG_LVL0
       write(*,"(a)") "Entering RHF only block"
       write(*,*)
#endif 
C For RHF rCC there are two blocks. One of them is the singlet spin-adapted 
C and the other is the jtriplet spin-adapted. In the following loop
C IMULT=1 singlet SA is built (stored in list 56) and INPUT=2
C triplet SA is built (stored in list 54). 

       DO IMULT = 1, 2
       LSTINT=18
       LSTINTB=25
       IF (IMULT .EQ. 1) THEN
          LISTQ=56
       ELSE IF (IMULT .EQ. 2) THEN
          LISTQ=54
       ENDIF 
       LISTTA=37
       LISTTB=34
       LISTWA=18
       LISTWB=21
       TAUTYP1='ABAB'
       DO IRREP=1,NIRREP 
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

        IF(i030+DISSYQ*NUMSYQ*IINTFP .ge. MAXCOR) then
        write(6,*)'@DWMBEJ-E insufficient memory', 
     &      i030+DISSYQ*NUMSYQ*IINTFP, MAXCOR
        CALL ERREX 
        ENDIF 

C There are no antisymmetrized T2 amplitudes.

        CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
        IF (IMULT .EQ. 1) THEN
           CALL GETLST(ICORE(I020),1,NMSYTA,1,IRREP,LISTTA)
CSSS           CALL GETLST(ICORE(I010),1,NMSYTA,1,IRREP,LISTTB)
CSSS           CALL SSCAL (NMSYTA*NMSYTA,TWO,ICORE(I020),1)
CSSS           CALL SAXPY (NMSYTA*NMSYTA,ONEM,ICORE(I010),1,ICORE(I020),1)
        ELSE IF (IMULT .EQ.2) THEN
           CALL GETLST(ICORE(I020),1,NMSYTA,1,IRREP,LISTTB)
        ENDIF 
        
C        IF (IMULT .EQ. 2) CALL SSCAL (NMSYTA*NMSYTA,TWO,
C     &                                ICORE(I020),1)
C
C The integrals are antisymmetric. So they get singlet spin-adapated 
C in cntract for the singlet spin-adpated path. 

        IF (IMULT .EQ. 1) THEN
            CALL CNTRCT(ICORE(I000),ICORE(I010),ICORE(I020),
     &                  DISSYQ,NUMSYQ,DSSYTA,NMSYTA,DSSYWA,NMSYWA,
     &                  IRREP,LISTWA,LISTWB,LISTTA,'WxT',ZILCH,-C4,1)
        ELSEIF (IMULT .EQ.2) THEN
C The triplet spin-adapted path. Also notice that minus sign in front
C (by using C4 instead of -C4).

            CALL RCC_TRP_CNTRCT(ICORE(I000),ICORE(I010),ICORE(I020),
     &                          DISSYQ,NUMSYQ,DSSYTA,NMSYTA,DSSYWA,
     &                          NMSYWA,IRREP,LISTWA,LISTWB,LISTTA,
     &                         'WxT',ZILCH,-C4,1)
        ENDIF 
        NSIZE=DISSYQ*NUMSYQ
#ifdef _DEBUG_LVLM
        Write(6,"(a,i4)") "The ABAB and BABA-target list",listq
        call checksum(SPCASE, ICORE(I000), NSIZE)
#endif 
C
C For rCCD we have antisymmetric integrals and we need to singlet and
c triplet spin-adapt them.
C
        CALL GETLST(ICORE(I010),1,NUMSYQ,1,IRREP,LSTINT)
        CALL GETLST(ICORE(I020),1,NUMSYQ,1,IRREP,LSTINTB)
        IF (IMULT .EQ. 1) THEN
           CALL SSCAL (DISSYQ*DISSYQ,TWO,ICORE(I010),1)
           CALL SAXPY (DISSYQ*DISSYQ,ONEM,ICORE(I020),1,ICORE(I010),1) 
           CALL SAXPY (DISSYQ*DISSYQ,ONE,ICORE(I010),1,ICORE(I000),1) 
#ifdef _DEBUG_LVL0
           call checksum("2J-E :",ICORE(I010),NUMSYQ*DISSYQ)
#endif 
        ELSEIF (IMULT .EQ. 2) THEN 
           CALL SSCAL (DISSYQ*DISSYQ,ONEM,ICORE(I020),1)
           CALL SAXPY (DISSYQ*DISSYQ,ONE,ICORE(I020),1,ICORE(I000),1) 
#ifdef _DEBUG_LVL0
        call checksum("-E   :",ICORE(I020),NUMSYQ*DISSYQ)
#endif 
        ENDIF 

        CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)

#ifdef _DEBUG_LVL0
        Write(6,"(a,i4)") "The ABAB and BABA-target list",listq
        call checksum(SPCASE//" :", ICORE(I000), NSIZE)
#endif

       ENDDO
       ENDDO

       ENDIF
    
      ENDIF 
C
      RETURN
      END
