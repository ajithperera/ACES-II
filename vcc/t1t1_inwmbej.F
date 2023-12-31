      SUBROUTINE T1T1_INWMBEJ(ICORE,MAXCOR,SPCASE,IUHF,FACT)
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION TWO,ONE,ONEM,HALF,ZILCH,TWOM,FACT

      CHARACTER*4 SPCASE,TAUTYP1
      DIMENSION ICORE(MAXCOR),IOFFT1(8,2)

      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC,HBAR_4LCCSD,ADC2
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

      DATA ZILCH /0.0/
      DATA ONE /1.0D0/
      DATA HALF/0.5D0/
      DATA TWO /2.0D0/
      DATA TWOM /-2.0D0/
      DATA ONEM/-1.0D0/

      NNP1O2(I)=I*(I+1)/2

      CALL GETT1(ICORE,MAXCOR,MXCOR,IUHF,IOFFT1)
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

        if(i030+DISSYQ*NUMSYQ*IINTFP .ge. MAXCOR) then
        write(6,*)'@DWMBEJ-E insufficient memory', 
     &      i030+DISSYQ*NUMSYQ*IINTFP, MAXCOR
        stop
        endif

        NSIZE=DISSYQ*NUMSYQ 
        CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
        CALL IZERO(ICORE(I020),NMSYTA*DSSYTA*IINTFP)

        CALL GETARGS(LISTTA,I1,I2)
        CALL F2TAU(ICORE(I020),ICORE(IOFFT1(1,I1)),
     &             ICORE(IOFFT1(1,I2)),DSSYTA,NMSYTA,
     &             POP(1,I1),POP(1,I2),VRT(1,I1),VRT(1,I2),
     &             IRREP,FACT,TAUTYP1)

        CALL CNTRCT(ICORE(I000),ICORE(I010),ICORE(I020),
     &              DISSYQ,NUMSYQ,DSSYTA,NMSYTA,DSSYWA,NMSYWA,
     &              IRREP,LISTWA,LISTWA,LISTTA,'WxT',ONE,ONEM,0)

        CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)

100    CONTINUE
C
C DO SPIN CASES ABAB AND BABA:
C
C       W(MbEj)=<Mb|Ej>-(1/2)*SUM [T(jn,fb)<Mn|Ef>-T(jN,bF)<MN||EF>]
C
C       W(mBeJ)=<mB|eJ>-(1/2)*SUM [T(JN,FB)<mN|eF>-T(Jn,Bf)<mn||ef>]
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
        TAUTYP1='BBBB'
       ELSEIF(SPCASE.EQ.'BABA')THEN 
        LSTINT=17
        LISTQ=57
        LISTQ1=56
        LISTTA=34
        LISTWA=17
        LISTTB=36
        LISTWB=20
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

        NSIZE=DISSYQ*NUMSYQ 
        CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
        CALL IZERO(ICORE(I020),NMSYTA*DSSYTA*IINTFP)

        CALL GETARGS(LISTTA,I1,I2)
        CALL F2TAU(ICORE(I020),ICORE(IOFFT1(1,I1)),
     &             ICORE(IOFFT1(1,I2)),DSSYTA,NMSYTA,
     &             POP(1,I1),POP(1,I2),VRT(1,I1),VRT(1,I2),
     &             IRREP,FACT,TAUTYP1)

        CALL CNTRCT(ICORE(I000),ICORE(I010),ICORE(I020),
     &              DISSYQ,NUMSYQ,DSSYTA,NMSYTA,DSSYWA,NMSYWA,
     &              IRREP,LISTWA,LISTWA,LISTTA,'WxT',ONE,ONE,0) 

        CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)

200    CONTINUE
C
C SPIN ADAPTED CODE FOR RHF "ABAB" SPIN CASE
C
       ELSE

#ifdef _DEBUG_LVLM
       write(*,"(a)") "Entering RHF only block"
       write(*,*)
#endif 
       LSTINT=18
       LSTINTB=25
       LISTQ=56
       LISTTA=37
       LISTTB=39
       LISTWA=18
       LISTWB=21
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

        NSIZE=DISSYQ*NUMSYQ
        CALL IZERO(ICORE(I000),DISSYQ*NUMSYQ*IINTFP)
        CALL IZERO(ICORE(I010),NMSYTA*DSSYTA*IINTFP)
        CALL IZERO(ICORE(I020),NMSYTA*DSSYTA*IINTFP)

        CALL GETARGS(LISTTA,I1,I2)
        CALL F2TAU(ICORE(I010),ICORE(IOFFT1(1,I1)),
     &             ICORE(IOFFT1(1,I2)),DSSYTA,NMSYTA,
     &             POP(1,I1),POP(1,I2),VRT(1,I1),VRT(1,I2),
     &             IRREP,FACT,TAUTYP1)

        CALL SAXPY (NMSYTA*NMSYTA,ONEM,ICORE(I010),1,ICORE(I020),1)

CSSS        If (Ispar .and. Coulomb) Then
CSSS        CALL CNTRCT_CONLY(ICORE(I000),ICORE(I010),ICORE(I020),
CSSS     &                     DISSYQ,NUMSYQ,DSSYTA,NMSYTA,DSSYWA,
CSSS     &                     NMSYWA,IRREP,LISTWA,LISTWB,LISTTA,
CSSS     &                     'WxT',ZILCH,ONEM,1)
CSSS        Else
        CALL CNTRCT(ICORE(I000),ICORE(I010),ICORE(I020),
     &              DISSYQ,NUMSYQ,DSSYTA,NMSYTA,DSSYWA,NMSYWA,
     &              IRREP,LISTWA,LISTWB,LISTTA,'WxT',ZILCH,ONEM,1)
CSSS        Endif 

        CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)
201    CONTINUE
C      
       ENDIF 
       
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

        NSIZE=DISSYQ*NUMSYQ
        CALL IZERO(ICORE(I000),IINTFP*DISSYQ*NUMSYQ)
        CALL IZERO(ICORE(I020),IINTFP*NMSYTA*DSSYTA)

        CALL GETARGS(LISTTA,I1,I2)
        CALL F2TAU(ICORE(I020),ICORE(IOFFT1(1,I1)),
     &             ICORE(IOFFT1(1,I2)),DSSYTA,NMSYTA,
     &             POP(1,I1),POP(1,I2),VRT(1,I1),VRT(1,I2),
     &             IRREP,FACT,TAUTYP1)

        CALL CNTRCT(ICORE(I000),ICORE(I010),ICORE(I020),
     &              DISSYQ,NUMSYQ,DSSYTA,NMSYTA,DSSYWA,NMSYWA,
     &              IRREP,LISTWA,LISTWA,LISTTA,'WxT',ONE,ONE,0)

#ifdef _DEBUG_LVLM
       write(*,"(a,I4)") "Writing RHF only ABBA to disk:",LISTQ
       call checksum("ABAB",ICORE(I000),NSIZE)
       write(*,*)
#endif
        CALL SUMSYM2(ICORE(I000),ICORE(I010),NSIZE,1,IRREP,LISTQ)
300    CONTINUE
      ENDIF

      RETURN
      END   
