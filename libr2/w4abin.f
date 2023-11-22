      SUBROUTINE W4ABIN(ICORE,MAXCOR,IUHF,SPCASE,TERM1,TERM2,
     &                  IOFFLIST)
C
C
C THIS ROUTINE CALCULATES THE CONTRIBUTIONS
C
C           Q(ifmn) = SUM T(ef,mo) * <io||en>
C                     o,e
C
C TO THE W(ifmn) INTERMEDIATE FOR SPIN CASES ABAB AND BABA
C AND INITIALIZES THE INTERMEDIATE LIST WITH  <MN||IF> -  Q(MNIF)
C NOTE MINUS SIGN).
C THE CONTRACTION IS A RATHER DIFFICULT ONE AND REQUIRES FANCY I/O
C TRICKS TO WRITE AS A MATRIX MULTIPLICATION.
C
C FOR THE SPECIFIC SPIN CASES, THE CONTRACTIONS ARE
C
C Q(IfMn) = SUM T(Ef,Mo) * <Io|En> + SUM T(Ef,On) * <IO||EM> 
C           o,E                      O,E
C
C         - SUM T(ef,on) * <Io|Me>        (ABAB)
C           o,e
C
C Q(iFmN) = SUM T(eF,mO) * <iO|eN> + SUM T(eF,oN) * <io||em> 
C           O,e                      o,e
C
C         - SUM T(EF,ON) * <iO|mE>        (BABA)
C           O,E
C
C ONLY THE ABAB SPIN CASE IS PERFORMED FOR RHF
C
C THE INTERMEDIATE LIST NUMBER IS 109 FOR SPIN CASE BABA
C                                 110 FOR SPIN CASE ABAB
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,ALPHA,BETA
      CHARACTER*4 SPCASE
      LOGICAL RHF,TERM1,TERM2
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      DATA ONE   /1.0/
      DATA ONEM  /-1.0/
      DATA ZILCH /0.0/
      BETA =ZILCH
      ALPHA=ONE
      RHF=.FALSE.
      IF(IUHF.EQ.0)RHF=.TRUE.
      IF(SPCASE.EQ.'ABAB')THEN
C 
C CODE FOR SPIN CASE ABAB.
C
       LSTTAR=110+IOFFLIST
       IF(.NOT.TERM2)GOTO 202
       LISTI1=9
       LISTI2=7
       LISTI3=10
       LISTT1=38
       LISTT2=37
       LISTT3=35
       IF(RHF)THEN
        LISTT1=39
        LISTI1=10
        LISTT3=34
       ENDIF
       ISCRSZ1=-1
       ISCRSZ2=-1
       DO 99 IRREP=1,NIRREP 
        ISCRSZ1=MAX(ISCRSZ1,IRPDPD(IRREP,ISYTYP(1,LISTI1)),
     &              IRPDPD(IRREP,ISYTYP(1,LSTTAR)))
        ISCRSZ2=MAX(ISCRSZ2,IRPDPD(IRREP,ISYTYP(1,LISTI2)),
     &              IRPDPD(IRREP,ISYTYP(1,LSTTAR)))
99     CONTINUE
       DO 100 IRREPDO=1,NIRREP
        TARDSZ=IRPDPD(IRREPDO,ISYTYP(1,LSTTAR))
        TARDIS=IRPDPD(IRREPDO,ISYTYP(2,LSTTAR))
C
C DO FIRST CONTRACTION
C
C        Q(IfMn) = SUM T(Ef,Mo) * <Io|En> 
C                  o,E                   
C
        INTDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTI1))
        INTDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTI1))
        NT2DSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTT1))
        NT2DIS=IRPDPD(IRREPDO,ISYTYP(2,LISTT1))
C
C TARGET MATRIX AS USED IN FIRST CONTRACTION
C
        TP1ADZ=IRPDPD(IRREPDO,18)
        TP1ADS=IRPDPD(IRREPDO,14)
C
C I MATRIX AS USED IN FIRST CONTRACTION
C
        IF(RHF)THEN
         TP1BDZ=IRPDPD(IRREPDO,18)
         TP1BDS=IRPDPD(IRREPDO,14)
        ELSE
         TP1BDZ=IRPDPD(IRREPDO,14)
         TP1BDS=IRPDPD(IRREPDO,11)
        ENDIF
C
C TARGET MATRIX AS USED IN SECOND CONTRACTION
C
        TP2ADZ=IRPDPD(IRREPDO,21)
        TP2ADS=IRPDPD(IRREPDO,10)
C
C I MATRIX AS USED IN SECOND CONTRACTION
C
        TP2BDZ=IRPDPD(IRREPDO,9)
        TP2BDS=IRPDPD(IRREPDO,21)
C
C TARGET MATRIX AS USED IN THIRD CONTRACTION
C
        TP3ADZ=TP2ADZ
        TP3ADS=TP2ADS
C
C INTEGRAL MATRIX AS USED IN THIRD CONTRACTION
C
        TP3BDZ=IRPDPD(IRREPDO,21)
        TP3BDS=IRPDPD(IRREPDO,10)

        I000=1
        I010=I000+IINTFP*TP1ADZ*TP1ADS
        CALL IZERO(ICORE,IINTFP*TP1ADZ*TP1ADS)
        I020=I010+IINTFP*TP1BDZ*TP1BDS
        I030=I020+IINTFP*ISCRSZ1
        IF(.NOT.RHF)THEN
C
C FOR UHF, READ <Io|En> INTEGRALS INTO AN I(In,Eo) MATRIX
C
         CALL FANCYGET(ICORE(I010),ICORE(I020),'FF','1432',
     &                 POP(1,1),POP(1,2),VRT(1,1),POP(1,2),
     &                 1,TP1BDS,2,IRREPDO,LISTI1)
C
C READ T2 AMPLITUDES INTO A T2(fM,Eo) MATRIX
C
         I030=I020+IINTFP*NT2DSZ*NT2DIS
         CALL GETLST(ICORE(I020),1,NT2DIS,1,IRREPDO,LISTT1)
C
C TRANSPOSE BRA INDICES TO FORM T2(Mf,Eo)
C
         I040=I030+IINTFP*NT2DIS
         I050=I040+IINTFP*NT2DIS
         I060=I050+IINTFP*NT2DIS
         CALL SYMTR3(IRREPDO,VRT(1,2),POP(1,1),NT2DSZ,NT2DIS,
     &               ICORE(I020),ICORE(I030),ICORE(I040),ICORE(I050))
C
C PERFORM THE MATRIX MULTIPLICATION
C                                           +
C             Z(Mf,In) = T(Mf,Eo) * I(In,Eo)
C
         BETA=ZILCH
         CALL XGEMM('N','T',NT2DSZ,TP1BDZ,NT2DIS,ALPHA,ICORE(I020),
     &              NT2DSZ,ICORE(I010),TP1BDZ,BETA,ICORE(I000),NT2DSZ)
C
C
C
        ELSEIF(RHF)THEN
C
C FOR RHF, READ THE <Oi|Ne> (= <oI|nE>) INTEGRALS INTO AN I(oE,nI)
C  MATRIX
C
         CALL FANCYGET(ICORE(I010),ICORE(I020),'FF','1432',
     &                 POP(1,1),POP(1,1),POP(1,1),VRT(1,1),
     &                 1,TP1BDS,2,IRREPDO,LISTI1)
C
C TRANSPOSE KET INDICES TO FORM I(oE,In)
C
         I030=I020+IINTFP*TP1BDZ
         I040=I030+IINTFP*TP1BDZ
         I050=I040+IINTFP*TP1BDZ
         CALL SYMTR1(IRREPDO,POP(1,1),POP(1,1),TP1BDZ,ICORE(I010),
     &               ICORE(I020),ICORE(I030),ICORE(I040))
C
C READ T2 AMPLITUDES INTO A T2(Eo,fM) MATRIX
C
         I030=I020+IINTFP*NT2DSZ*NT2DIS
         CALL GETLST(ICORE(I020),1,NT2DIS,1,IRREPDO,LISTT1)
C
C TRANSPOSE BRA AND KET INDICES TO FORM T2(oE,Mf)
C
         I040=I030+IINTFP*MAX(NT2DSZ,NT2DIS)
         I050=I040+IINTFP*MAX(NT2DSZ,NT2DIS)
         I060=I050+IINTFP*MAX(NT2DSZ,NT2DIS)
         CALL SYMTR1(IRREPDO,VRT(1,1),POP(1,1),NT2DSZ,ICORE(I020),
     &               ICORE(I030),ICORE(I040),ICORE(I050))
         CALL SYMTR3(IRREPDO,VRT(1,1),POP(1,1),NT2DSZ,NT2DIS,
     &               ICORE(I020),ICORE(I030),ICORE(I040),ICORE(I050))
C
C PERFORM THE MATRIX MULTIPLICATION
C                                +           
C             Z(Mf,In) = T(oE,Mf) * I(oE,In)
C
         BETA=ZILCH
         CALL XGEMM('T','N',NT2DIS,TP1BDS,NT2DSZ,ALPHA,ICORE(I020),
     &              NT2DSZ,ICORE(I010),TP1BDZ,BETA,ICORE(I000),NT2DIS)
C
        ENDIF
C
C
C NOW WRITE THE PRODUCT TO DISK AS Z(Mn,If)
C
        CALL FANCYPUT(ICORE(I000),ICORE(I010),'FF','1432','N','S',
     &                 ONE,POP(1,1),POP(1,2),POP(1,1),VRT(1,2),
     &                 1,TP1ADS,2,IRREPDO,LSTTAR)
C
C
c
C
C
C NOW PERFORM SECOND CONTRACTION
C
C              Q(IfMn) = SUM T(Ef,On) * <OI||ME> (ABAB)
C                        O,E
C
        INTDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTI2))
        INTDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTI2))
        NT2DSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTT2))
        NT2DIS=IRPDPD(IRREPDO,ISYTYP(2,LISTT2))
        I000=1
        I010=I000+IINTFP*TP2ADZ*TP2ADS
        CALL IZERO(ICORE,IINTFP*TP2ADZ*TP2ADS)
        I020=I010+IINTFP*TP2BDZ*TP2BDS
        I030=I020+IINTFP*ISCRSZ2
C
C READ <OI||ME> INTEGRALS INTO AN I(OE,MI) MATRIX AND TRANSPOSE
C  BRA INDICES GIVING I(EO,MI).
C
        CALL FANCYGET(ICORE(I010),ICORE(I020),'PF','1432',
     &                POP(1,1),POP(1,1),POP(1,1),VRT(1,1),
     &                1,TP2BDS,2,IRREPDO,LISTI2)
        I030=I020+IINTFP*TP2BDS
        I040=I030+IINTFP*TP2BDS
        I050=I040+IINTFP*TP2BDS
        CALL SYMTR3(IRREPDO,POP(1,1),VRT(1,1),TP2BDZ,TP2BDS,
     &              ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040))
C
C READ T2 AMPLITUDES INTO A T2(EO,fn) MATRIX
C
        I030=I020+IINTFP*NT2DSZ*NT2DIS
        CALL GETLST(ICORE(I020),1,NT2DIS,1,IRREPDO,LISTT2)
C
C NOW PERFORM THE MATRIX MULTIPLICATION
C                             +
C          Z(MI,fn) = I(EO,MI) * T(EO,fn)
C
        BETA=ZILCH
        CALL XGEMM('T','N',TP2BDS,NT2DIS,TP2BDZ,ALPHA,ICORE(I010),
     &             TP2BDZ,ICORE(I020),NT2DSZ,BETA,ICORE(I000),
     &             TP2BDS)
C
C
C
C NOW PERFORM THIRD CONTRACTION AND ACCUMULATE INTO BOTTOM OF ICORE.
C
C              Q(IfMn) = SUM T(ef,on) * <Io|Me> (ABAB)
C                        o,e
C
        INTDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTI3))
        INTDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTI3))
        NT2DSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTT3))
        NT2DIS=IRPDPD(IRREPDO,ISYTYP(2,LISTT3))
        I000=1
        I010=I000+IINTFP*TP3ADZ*TP3ADS
        I020=I010+IINTFP*TP3BDZ*TP3BDS
        I030=I020+IINTFP*ISCRSZ2
C
C READ <Io|Me> INTEGRALS INTO AN I(IM,oe) MATRIX AND TRANSPOSE BRA AND KET 
C  INDICES TO GIVE I(MI,eo)
C
        CALL FANCYGET(ICORE(I010),ICORE(I020),'FF','1324',
     &                POP(1,1),POP(1,2),POP(1,1),VRT(1,2),
     &                1,TP3BDS,2,IRREPDO,LISTI3)
        I030=I020+IINTFP*MAX(TP3BDZ,TP3BDS)
        I040=I030+IINTFP*MAX(TP3BDZ,TP3BDS)
        I050=I040+IINTFP*MAX(TP3BDZ,TP3BDS)
        CALL SYMTR1(IRREPDO,POP(1,2),VRT(1,2),TP3BDZ,ICORE(I010),
     &              ICORE(I020),ICORE(I030),ICORE(I040))
        CALL SYMTR3(IRREPDO,POP(1,1),POP(1,1),TP3BDZ,TP3BDS,
     &              ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040))
C
C READ T2 AMPLITUDES INTO A T2(eo,fn) MATRIX.
C        
        I030=I020+IINTFP*NT2DSZ*NT2DIS
        CALL GETLST(ICORE(I020),1,NT2DIS,1,IRREPDO,LISTT3)
C
C PERFORM MATRIX MULTIPLICATION
C
C         Z(MI,fn) = I(MI,eo) * T(eo,fn)
C
        BETA=ONE
        CALL XGEMM('N','N',TP3BDZ,NT2DIS,TP3BDS,ALPHA,ICORE(I010),
     &             TP3BDZ,ICORE(I020),NT2DSZ,BETA,ICORE(I000),
     &             TP3BDZ)
C
C NOW AUGMENT THE LISTS ON DISK, WRITING THIS AND PRECEDING 
C   CONTRIBUTION AS Z(Mn,If)
C
        CALL FANCYPUT(ICORE(I000),ICORE(I010),'FF','1423','N','S',
     &                ONE,POP(1,1),POP(1,2),POP(1,1),VRT(1,2),
     &                1,TP3ADS,2,IRREPDO,LSTTAR)
100    CONTINUE
C
C
      ELSEIF(SPCASE.EQ.'BABA')THEN
C
C CODE FOR SPIN CASE BABA.
C
       LSTTAR=109+IOFFLIST
       IF(.NOT.TERM2)GOTO 202
       LISTI1=10
       LISTI2=8
       LISTI3=9
       LISTT1=39
       LISTT2=36
       LISTT3=34
       ISCRSZ1=-1
       ISCRSZ2=-1
       DO 199 IRREP=1,NIRREP 
        ISCRSZ1=MAX(ISCRSZ1,IRPDPD(IRREP,ISYTYP(1,LISTI1)),
     &              IRPDPD(IRREP,ISYTYP(1,LSTTAR)))
        ISCRSZ2=MAX(ISCRSZ2,IRPDPD(IRREP,ISYTYP(1,LISTI2)),
     &              IRPDPD(IRREP,ISYTYP(1,LSTTAR)))
199     CONTINUE
       DO 200 IRREPDO=1,NIRREP
        TARDSZ=IRPDPD(IRREPDO,ISYTYP(1,LSTTAR))
        TARDIS=IRPDPD(IRREPDO,ISYTYP(2,LSTTAR))
C
C DO FIRST CONTRACTION
C
C        Q(iFmN) = SUM T(eF,mO) * <Oi|Ne>
C                  O,e                   
C
        INTDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTI1))
        INTDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTI1))
        NT2DSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTT1))
        NT2DIS=IRPDPD(IRREPDO,ISYTYP(2,LISTT1))
C
C TARGET MATRIX AS USED IN FIRST CONTRACTION
C
        TP1ADZ=IRPDPD(IRREPDO,14)
        TP1ADS=IRPDPD(IRREPDO,11)
C
C I MATRIX AS USED IN FIRST CONTRACTION
C
        TP1BDZ=IRPDPD(IRREPDO,18)
        TP1BDS=IRPDPD(IRREPDO,14)
C
C TARGET MATRIX AS USED IN SECOND CONTRACTION
C
        TP2ADZ=IRPDPD(IRREPDO,16)
        TP2ADS=IRPDPD(IRREPDO,22)
C
C I MATRIX AS USED IN SECOND CONTRACTION
C
        TP2BDZ=IRPDPD(IRREPDO,17)
        TP2BDS=IRPDPD(IRREPDO,22)
C
C TARGET MATRIX AS USED IN THIRD CONTRACTION
C
        TP3ADZ=TP2ADZ
        TP3ADS=TP2ADS
C
C I MATRIX AS USED IN THIRD CONTRACTION
C
        TP3BDZ=TP3ADZ
        TP3BDS=TP2BDS
C
        I000=1
        I010=I000+IINTFP*TP1ADZ*TP1ADS
        I020=I010+IINTFP*TP1BDZ*TP1BDS
        I030=I020+IINTFP*ISCRSZ1
C
C READ <Oi|Ne> INTEGRALS INTO AN I(Oe,Ni) MATRIX
C
        CALL FANCYGET(ICORE(I010),ICORE(I020),'FF','1432',
     &                POP(1,1),POP(1,2),POP(1,1),VRT(1,2),
     &                1,TP1BDS,2,IRREPDO,LISTI1)
C
C READ T2 AMPLITUDES INTO A T2(Fm,eO) MATRIX
C
        I030=I020+IINTFP*NT2DSZ*NT2DIS
        CALL GETLST(ICORE(I020),1,NT2DIS,1,IRREPDO,LISTT1)
C
C TRANSPOSE KET INDICES TO FORM T2(Fm,Oe)
C
        I040=I030+IINTFP*NT2DSZ
        I050=I040+IINTFP*NT2DSZ
        I060=I050+IINTFP*NT2DSZ
        CALL SYMTR1(IRREPDO,VRT(1,2),POP(1,1),NT2DSZ,ICORE(I020),
     &              ICORE(I030),ICORE(I040),ICORE(I050))
C
C PERFORM THE MATRIX MULTIPLICATION
C                                +          +
C             Z(Ni,Fm) = I(Oe,Ni) * T(Fm,Oe)
C
        BETA=ZILCH
        CALL XGEMM('T','T',TP1BDS,NT2DSZ,TP1BDZ,ALPHA,ICORE(I010),
     &              TP1BDZ,ICORE(I020),NT2DSZ,BETA,ICORE(I000),TP1BDS)
C
C
C NOW WRITE THE PRODUCT TO DISK AS Z(Nm,Fi)
C
         CALL FANCYPUT(ICORE(I000),ICORE(I010),'FF','1432','N','S',
     &                 ONE,POP(1,1),POP(1,2),VRT(1,1),POP(1,2),
     &                 1,TP1ADS,2,IRREPDO,LSTTAR)
C
C
C NOW PERFORM SECOND CONTRACTION
C
C              Q(iFmN) = SUM  T(eF,oN) * <oi||me> (BABA)
C                        o,e
C
        INTDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTI2))
        INTDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTI2))
        NT2DSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTT2))
        NT2DIS=IRPDPD(IRREPDO,ISYTYP(2,LISTT2))
        I000=1
        I010=I000+IINTFP*TP2ADZ*TP2ADS
        I020=I010+IINTFP*TP2BDZ*TP2BDS
        I030=I020+IINTFP*ISCRSZ2
C
C READ <oi||me> INTEGRALS INTO AN I(oe,mi) MATRIX AND TRANSPOSE
C  BRA AND KET INDICES TO FORM I(eo,im).
C
        CALL FANCYGET(ICORE(I010),ICORE(I020),'PF','1432',
     &                POP(1,2),POP(1,2),POP(1,2),VRT(1,2),
     &                1,TP2BDS,2,IRREPDO,LISTI2)
        I030=I020+IINTFP*MAX(TP2BDS,TP2BDZ)
        I040=I030+IINTFP*MAX(TP2BDS,TP2BDZ)
        I050=I040+IINTFP*MAX(TP2BDS,TP2BDZ)
        CALL SYMTR3(IRREPDO,POP(1,2),VRT(1,2),TP2BDZ,TP2BDS,
     &              ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040))
        CALL SYMTR1(IRREPDO,POP(1,2),POP(1,2),TP2BDZ,
     &              ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040))
C
C READ T2 AMPLITUDES INTO A T2(eo,FN) MATRIX AND TRANSPOSE KET INDICES
C  TO FORM T2(eo,NF)
C
        I030=I020+IINTFP*NT2DSZ*NT2DIS
        CALL GETLST(ICORE(I020),1,NT2DIS,1,IRREPDO,LISTT2)
        I040=I030+IINTFP*NT2DSZ
        I050=I040+IINTFP*NT2DSZ
        I060=I050+IINTFP*NT2DSZ
        CALL SYMTR1(IRREPDO,VRT(1,1),POP(1,1),NT2DSZ,ICORE(I020),
     &              ICORE(I030),ICORE(I040),ICORE(I050))
C
C NOW PERFORM THE MATRIX MULTIPLICATION
C                             +         
C          Z(NF,im) = T(eo,NF) * I(eo,im)
C
        BETA=ZILCH
        CALL XGEMM('T','N',NT2DIS,TP2BDS,NT2DSZ,ALPHA,ICORE(I020),
     &             NT2DSZ,ICORE(I010),TP2BDZ,BETA,ICORE(I000),
     &             NT2DIS)
C
C
C NOW PERFORM THIRD CONTRACTION
C
C              Q(iFmN) = SUM T(EF,ON) * <Oi|Em>        (BABA)
C                        O,E
C
        INTDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTI3))
        INTDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTI3))
        NT2DSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTT3))
        NT2DIS=IRPDPD(IRREPDO,ISYTYP(2,LISTT3))
        I000=1
        I010=I000+IINTFP*TP3ADZ*TP3ADS
        I020=I010+IINTFP*TP3BDZ*TP3BDS
        I030=I020+IINTFP*ISCRSZ2
C
C READ <Oi|Em> INTEGRALS INTO AN I(OE,im) MATRIX AND TRANSPOSE
C  BRA INDICES TO FORM I(EO,im)
C
        CALL FANCYGET(ICORE(I010),ICORE(I020),'FF','1324',
     &                POP(1,1),POP(1,2),VRT(1,1),POP(1,2),
     &                1,TP3BDS,2,IRREPDO,LISTI3)
        I030=I020+IINTFP*TP3BDS
        I040=I030+IINTFP*TP3BDS
        I050=I040+IINTFP*TP3BDS
        CALL SYMTR3(IRREPDO,POP(1,1),VRT(1,1),TP3BDZ,TP3BDS,
     &              ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040))
C
C READ T2 AMPLITUDES INTO A T2(EO,FN) MATRIX.
C
        I030=I020+IINTFP*NT2DIS*NT2DSZ
        CALL GETLST(ICORE(I020),1,NT2DIS,1,IRREPDO,LISTT3)
C
C TRANSPOSE KET INDICES TO FORM T(EO,NF)
C
        I040=I030+IINTFP*NT2DSZ
        I050=I040+IINTFP*NT2DSZ
        I060=I050+IINTFP*NT2DSZ
        CALL SYMTR1(IRREPDO,VRT(1,1),POP(1,1),NT2DSZ,ICORE(I020),
     &              ICORE(I030),ICORE(I040),ICORE(I050))
C
C NOW FORM THE MATRIX PRODUCT AND ACCUMULATE WITH LAST TERM
C                                     +
C                  Z(NF,im) = T(EO,NF) * I(EO,im)
C
        BETA=ONE
        CALL XGEMM('T','N',NT2DIS,TP3BDS,NT2DSZ,ALPHA,ICORE(I020),
     &             NT2DSZ,ICORE(I010),TP3BDZ,BETA,ICORE(I000),
     &             NT2DIS)
C
C NOW AUGMENT THE LISTS ON DISK, WRITING THIS AND PREVIOUS TERM
C   AS Z(Nm,Fi)
C
        CALL FANCYPUT(ICORE(I000),ICORE(I010),'FF','1423','N','S',
     &                ONE,POP(1,1),POP(1,2),VRT(1,1),POP(1,2),
     &                1,TP2ADS,2,IRREPDO,LSTTAR)
C
200    CONTINUE
      ENDIF
C
C NOW AUGMENT THE TARGET LIST WITH THE CORRESPONDING INTEGRALS
C
202   DO 1000 IRREP=1,NIRREP
       TARDSZ=IRPDPD(IRREP,ISYTYP(1,LSTTAR))
       TARDIS=IRPDPD(IRREP,ISYTYP(2,LSTTAR))
       I000=1
       I010=I000+IINTFP*TARDIS*TARDSZ
       I020=I010+IINTFP*TARDIS*TARDSZ
       IF(TERM1)THEN
        CALL GETLST(ICORE(I000),1,TARDIS,1,IRREP,LSTTAR-100-IOFFLIST)
       ELSE
        CALL ZERO(ICORE(I000),TARDIS*TARDSZ)
       ENDIF
       IF(TERM2)THEN
        CALL GETLST(ICORE(I010),1,TARDIS,1,IRREP,LSTTAR)
       ELSE
        CALL ZERO(ICORE(I010),TARDIS*TARDSZ)
       ENDIF
       CALL SAXPY(TARDIS*TARDSZ,ONEM,ICORE(I010),1,ICORE(I000),1)
       CALL PUTLST(ICORE(I000),1,TARDIS,1,IRREP,LSTTAR)
1000  CONTINUE
C
C WHEW!  WE'RE FINALLY DONE WITH ALL OF THIS AWFUL CRAP!
C
      RETURN
      END
