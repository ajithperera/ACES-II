      SUBROUTINE W5ABIN(ICORE,MAXCOR,IUHF,SPCASE,TERM1,TERM2)
C
C
C THIS ROUTINE CALCULATES THE CONTRIBUTIONS
C
C           Z(efam) =  - P(ef) SUM T(fg,mn) * <ga||en>
C                              n,g
C
C TO THE W(efam) INTERMEDIATE FOR SPIN CASES ABAB AND BABA
C AND INITIALIZES THE INTERMEDIATE LIST WITH  <EF||AM> -  Z(EFAM)
C NOTE MINUS SIGN).
C THE CONTRACTION IS A RATHER DIFFICULT ONE AND REQUIRES FANCY I/O
C TRICKS TO WRITE AS A MATRIX MULTIPLICATION.
C
C FOR THE SPECIFIC SPIN CASES, THE CONTRACTIONS ARE
C
C Z(EfAm) =  SUM T(fg,mn) * <gA|nE> - SUM T(fG,mN) * <GA||EN>
C            n,g                      N,G
C
C          - SUM T(Eg,Nm) * <Ag|Nf>            (ABAB)
C            N,g
C
C Z(eFaM) = SUM T(FG,MN) * <Ga|Ne> - SUM T(Fg|Mn> * <ga||en>
C           N,G                      n,g
C
C         - SUM T(Ge,Mn) * <Ga|Fn>             (BABA)
C           n,G
C
C ONLY THE ABAB SPIN CASE IS PERFORMED FOR RHF
C
C SPIN ADAPTED CODE IS USED FOR RHF
C
C Z(eF,Am) = SUM (2 T(Fg,Mn) - T(Fg,Nm) (<gA|nE> - 1/2 <aG|nE>
C
C                - 1/2 t(Fg,Nm) <Ag|nE> - T(Eg,Nm) <Ag|Nf>
C
C THE INTERMEDIATE LIST NUMBER IS 129 FOR SPIN CASE BABA
C                                 130 FOR SPIN CASE ABAB
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,ALPHA,BETA,FACT,HALF,HALFM
      CHARACTER*4 SPCASE
      LOGICAL RHF,TERM1,TERM2
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      DATA ONE   /1.0D0/
      DATA ONEM  /-1.0D0/
      DATA HALF,HALFM /.5D0,-.5D0/
      DATA ZILCH /0.0D0/
      BETA =ZILCH
      ALPHA=ONE
      RHF=.FALSE.
      IF(IUHF.EQ.0)RHF=.TRUE.
C
C DETERMINE SCRATCH SPACE
C
       ISCRSZ1=0
       ISCRSZ2=0
       DO 10 IRREP=1,NIRREP
        ISCRSZ1=MAX(ISCRSZ1,IRPDPD(IRREP,ISYTYP(1,30)))
        ISCRSZ1=MAX(ISCRSZ1,IRPDPD(IRREP,ISYTYP(1,30-IUHF)))
        ISCRSZ2=MAX(ISCRSZ2,IRPDPD(IRREP,19),
     &              IRPDPD(IRREP,19+IUHF))
10     CONTINUE
       ITMPSIZ=MAX(ISCRSZ1,ISCRSZ2)*10
C
      IF(SPCASE.EQ.'ABAB')THEN
C 
C Z(EfAm) =  SUM T(fg,mn) * <Ag|En> - SUM T(Gf,Nm) * <GA||EN>
C            n,g                      N,G
C
C          - SUM T(Eg,Nm) * <Ag|Nf>            (ABAB)
C            N,g
C
C CODE FOR SPIN CASE ABAB.
C
       LSTTAR=130
       IF(.NOT.TERM2)GOTO 202
       LISTI1=30
       LISTI2=27
       LISTI3=30-IUHF
       LISTT1=34+IUHF
       LISTT1A=37
       LISTT2=37
       LISTT3=39-IUHF
C
       DO 100 IRREPDO=1,NIRREP
C
C TARGET MATRIX AS USED IN FIRST CONTRACTION
C
        TP1ADZ=IRPDPD(IRREPDO,19)
        TP1ADS=IRPDPD(IRREPDO,10)
C
C I MATRIX AS USED IN FIRST CONTRACTION
C
        TP1BDZ=IRPDPD(IRREPDO,19)
        TP1BDS=IRPDPD(IRREPDO,10)
C
C TARGET MATRIX AS USED IN SECOND CONTRACTION
C
        TP2ADZ=TP1ADZ
        TP2ADS=TP1ADS
C
C I MATRIX AS USED IN SECOND CONTRACTION
C
        TP2BDZ=IRPDPD(IRREPDO,9)
        TP2BDS=IRPDPD(IRREPDO,19)
C
C TARGET MATRIX AS USED IN THIRD CONTRACTION
C
        TP3ADZ=IRPDPD(IRREPDO,11)
        TP3ADS=IRPDPD(IRREPDO,13)
C
C INTEGRAL MATRIX AS USED IN THIRD CONTRACTION
C
        IF(RHF)THEN
         TP3BDZ=IRPDPD(IRREPDO,12)
         TP3BDS=IRPDPD(IRREPDO,13)
        ELSE
         TP3BDZ=IRPDPD(IRREPDO,13)
         TP3BDS=IRPDPD(IRREPDO,12)
        ENDIF
C
C DO FIRST CONTRACTION
C
C        Z(EfAm) = SUM T(fg,mn) * <Ag|En>
C                  n,g
C
C ALPHA SET TO -1 BECAUSE OF T2(AA) AND T2(BB) RING ORDERED STORAGE
C  MODE
C
        ALPHA=ONEM
        BETA=ZILCH
        INTDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTI1))
        INTDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTI1))
        NT2DSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTT1))
        NT2DIS=IRPDPD(IRREPDO,ISYTYP(2,LISTT1))
        MAXTP =MAX(TP1BDZ,TP1BDS,ITMPSIZ)
        I000=1
        I010=I000+IINTFP*MAX(TP1ADZ*TP1ADS,NT2DIS*NT2DSZ*(1-IUHF))
        ITMP1=I010+IINTFP*NT2DSZ*NT2DIS
        ITMP2=ITMP1+IINTFP*MAXTP
        ITMP3=ITMP2+IINTFP*MAXTP
        I020=ITMP3+IINTFP*MAXTP
        I030=I020+IINTFP*MAXTP
        IF(I030.GE.MAXCOR) CALL INSMEM('W5ABIN',I030,MAXCOR)
        MXCOR=MAXCOR-I030+1
C
C READ T2 AMPLITUDES INTO A T2(gn,fm) MATRIX
C
        CALL GETLST(ICORE(I010),1,NT2DIS,1,IRREPDO,LISTT1)
C
C SPIN ADAPT THE AMPLITUDES FOR RHF
C
        IF(RHF) THEN
C
C    T2(gn,fm) + T2(gn,FM) (34 - 37), THERE IS A NEGATIVE SIGN IN 34 AND OVERALL WE
C                                     KEEP A NEGATIVE SIGN
C
         CALL GETLST(ICORE(I000),1,NT2DIS,1,IRREPDO,LISTT1A)
         CALL SSCAL(NT2DIS*NT2DSZ,HALF,ICORE(I010),1) 
         CALL SAXPY(NT2DIS*NT2DSZ,HALFM,ICORE(I000),1,ICORE(I010),1)
        ENDIF 
        CALL IZERO(ICORE,IINTFP*TP1ADS*TP1ADZ)
C  
C  READ THE <Ag|En> INTEGRALS AS AN I(AE,gn) MATRIX AND
C   TRANSPOSE BRA INDICES GIVING I(EA,gn)
C
        IF(TP1BDZ.NE.0)THEN
         NINCOR=MXCOR/(TP1BDZ*IINTFP)
        ELSE
         NINCOR=TP1BDS
        ENDIF
        NLEFT =TP1BDS
        NFIRST=1
        FACT=ONE
        IOFT=I010
1       NREAD=MIN(NLEFT,NINCOR)
        CALL FANCYGET1(ICORE(I020),ICORE(ITMP1),'FF','1324',
     &                VRT(1,1),VRT(1,2),VRT(1,1),POP(1,2),
     &                NFIRST,NREAD,2,IRREPDO,LISTI1,RHF,ICORE(ITMP2),
     &                ICORE(ITMP3))
        CALL SYMTR3(IRREPDO,VRT(1,1),VRT(1,1),TP1BDZ,NREAD,
     &              ICORE(I020),ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
        CALL XGEMM('N','N',TP1BDZ,NT2DIS,NREAD,ONEM,ICORE(I020),
     &             TP1BDZ,ICORE(IOFT),NT2DSZ,FACT,ICORE(I000),TP1BDZ)
        NFIRST=NFIRST+NREAD
        NLEFT =NLEFT -NREAD
        IOFT  =IOFT  +NREAD*IINTFP
        IF(NLEFT.NE.0)GOTO 1
c        CALL GETLST(ICORE(I020),NFIRST,NREAD,1,IRREP,LISTWA)
c        CALL SYMTR3(IRREPDO,VRT(1,1),VRT(1,1),TP1BDZ,TP1DIS,
c     &              ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040))
C
C PERFORM THE MATRIX MULTIPLICATION
C                                          
C             Z(EA,fm) = I(EA,gn) * T(gn,fm)
C
c         CALL XGEMM('N','N',TP1BDZ,NT2DIS,TP1BDS,ALPHA,ICORE(I010),
c     &              TP1BDZ,ICORE(I020),NT2DSZ,BETA,ICORE(I000),
c     &              TP1BDZ)
C
C DON'T WRITE PRODUCT, BUT RATHER ACCUMULATE IN NEXT CONTRACTION 
C
C NOW PERFORM SECOND CONTRACTION
C
C Z(EfAm) =   - SUM T(fG,mN) * <GA||EN>
C               N,G
C
C (UHF ONLY)
C
        IF(.NOT.RHF) THEN
C
        ALPHA=ONEM
        BETA=ONE
        INTDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTI2))
        INTDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTI2))
        NT2DSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTT2))
        NT2DIS=IRPDPD(IRREPDO,ISYTYP(2,LISTT2))
        ITMP1=I010+IINTFP*NT2DSZ*NT2DIS
        MAXTP=MAX(TP2BDS,TP2BDZ,ITMPSIZ)
        I020=ITMP1+IINTFP*MAXTP
        I030=I020+IINTFP*MAXTP
        IF(I030.GE.MAXCOR) CALL INSMEM('W5ABIN',I030,MAXCOR)
        MXCOR=MAXCOR-I030+1
C
C READ T2 AMPLITUDES INTO A T2(GN,fm) MATRIX
C
        CALL GETLST(ICORE(I010),1,NT2DIS,1,IRREPDO,LISTT2)
C  
C  READ THE <Ag|En> INTEGRALS AS AN I(AE,gn) MATRIX AND
C   TRANSPOSE BRA INDICES GIVING I(EA,gn)
C
C READ <GA||EN> INTEGRALS INTO AN I(GN,EA) MATRIX.
C
c        CALL FANCYGET(ICORE(I010),ICORE(I020),'PF','1432',
c     &                VRT(1,1),VRT(1,1),VRT(1,1),POP(1,1) ,
c     &                1,TP2BDS,2,IRREPDO,LISTI2)
C
C NOW PERFORM THE MATRIX MULTIPLICATION
C                             +
C          Z(EA,fm) = I(GN,EA) * T(GN,fm)
C
c        CALL XGEMM('T','N',TP2BDS,NT2DIS,TP2BDZ,ALPHA,ICORE(I010),
c     &             TP2BDZ,ICORE(I020),NT2DSZ,BETA,ICORE(I000),
c     &             TP2BDS)
        IF(TP2BDZ.NE.0)THEN
         NINCOR=MXCOR/(TP2BDZ*IINTFP)
        ELSE
         NINCOR=TP2BDS
        ENDIF
        NLEFT =TP2BDS
        NFIRST=1
        FACT=BETA
        IOFZ=I000
2       NREAD=MIN(NLEFT,NINCOR)
        CALL FANCYGET(ICORE(I020),ICORE(ITMP1),'PF','1432',
     &                VRT(1,1),VRT(1,1),VRT(1,1),POP(1,1),
     &                NFIRST,NREAD,2,IRREPDO,LISTI2)
        CALL XGEMM('T','N',NREAD,NT2DIS,TP2BDZ,ONEM,ICORE(I020),
     &             TP2BDZ,ICORE(I010),NT2DSZ,FACT,ICORE(IOFZ),TP2BDS)
        NFIRST=NFIRST+NREAD
        NLEFT =NLEFT-NREAD
        IOFZ=IOFZ+NREAD*IINTFP
        IF(NLEFT.NE.0)GOTO 2
        ENDIF
C
C NOW WRITE PRODUCT + PREVIOUS TERM TO DISK AS Z(Ef,Am)
C
        CALL FANCYPUT(ICORE(I000),ICORE(I010),'FF','1324','N','S',
     &                ONE,VRT(1,1),VRT(1,2),VRT(1,1),POP(1,2),
     &                1,TP2ADS,2,IRREPDO,LSTTAR)
C
C
C NOW PERFORM THIRD CONTRACTION. 
C
C  Z(EfAm) = - SUM T(Eg,Nm) * <Ag|Nf>            (ABAB)
C              N,g
C
        BETA=ZILCH
        ALPHA=ONEM
        INTDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTI3))
        INTDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTI3))
        NT2DSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTT3))
        NT2DIS=IRPDPD(IRREPDO,ISYTYP(2,LISTT3))
        I000=1
        I010=I000+IINTFP*TP3ADZ*TP3ADS
        ITMP1=I010+IINTFP*NT2DSZ*NT2DIS
        MAXTP=MAX(TP3BDS,TP3BDZ,ITMPSIZ)
        I020=ITMP1+IINTFP*MAXTP
        I030=I020+IINTFP*MAXTP
        IF(I030.GE.MAXCOR) CALL INSMEM('W5ABIN',I030,MAXCOR)
        MXCOR=MAXCOR-I030+1
        CALL IZERO(ICORE,IINTFP*TP3ADZ*TP3ADS)
C
C READ T2 AMPLITUDES INTO A T2(gN,Em) MATRIX.
C        
        CALL GETLST(ICORE(I010),1,NT2DIS,1,IRREPDO,LISTT3)
C
C READ <Ag|Nf> INTEGRALS INTO A I(Af,gN) MATRIX.
C
        IF(.NOT.RHF)THEN
C                            +          +
C         Z(Em,Af) = T(gN,Em) * I(Af,gN) [UHF]
C
c         CALL FANCYGET(ICORE(I010),ICORE(I020),'FF','1423',
c     &                 VRT(1,1),VRT(1,2),POP(1,1),VRT(1,2),
c     &                 1,TP3BDS,2,IRREPDO,LISTI3)
c         CALL XGEMM('T','T',NT2DIS,TP3BDZ,NT2DSZ,ALPHA,ICORE(I020),
c     &              NT2DSZ,ICORE(I010),TP3BDZ,BETA,ICORE(I000),
c     &              NT2DIS)
c
         IF(TP3BDZ.NE.0)THEN
          NINCOR=MXCOR/(TP3BDZ*IINTFP)
         ELSE
          NINCOR=TP3BDS
         ENDIF
         NLEFT =TP3BDS
         NFIRST=1
         FACT=BETA
         IOFT2=I010
3        NREAD=MIN(NLEFT,NINCOR)
         CALL FANCYGET(ICORE(I020),ICORE(ITMP1),'FF','1423',
     &                 VRT(1,1),VRT(1,2),POP(1,1),VRT(1,2),
     &                 NFIRST,NREAD,2,IRREPDO,LISTI3)
         CALL XGEMM('T','T',NT2DIS,TP3BDZ,NREAD,ALPHA,ICORE(IOFT2),
     &              NT2DSZ,ICORE(I020),TP3BDZ,FACT,ICORE(I000),NT2DIS)
         FACT=ONE
         NFIRST=NFIRST+NREAD
         NLEFT =NLEFT-NREAD
         IOFT2=IOFT2+NREAD*IINTFP
         IF(NLEFT.NE.0)GOTO 3
        ELSE
C
C FOR RHF, <Ag|Nf> = <Ga|Fn>, SO THIS READ GIVES US I(gN,Af).
C
C                            +        
C         Z(Em,Af) = T(gN,Em) * I(gN,Af) [RHF]
C
c         CALL FANCYGET(ICORE(I010),ICORE(I020),'FF','1423',
c     &                 VRT(1,1),VRT(1,1),VRT(1,1),POP(1,1),
c     &                 1,TP3BDS,2,IRREPDO,LISTI3)
c         CALL XGEMM('T','N',NT2DIS,TP3BDS,NT2DSZ,ALPHA,ICORE(I020),
c     &              NT2DSZ,ICORE(I010),TP3BDZ,BETA,ICORE(I000),
c     &              NT2DIS)
         IF(TP3BDZ.NE.0)THEN
          NINCOR=MXCOR/(TP3BDZ*IINTFP)
         ELSE
          NINCOR=TP3BDS
         ENDIF
         NLEFT =TP3BDS
         NFIRST=1
         FACT=BETA
         IOFZ=I000
4        NREAD=MIN(NLEFT,NINCOR)
         CALL FANCYGET(ICORE(I020),ICORE(ITMP1),'FF','1423',
     &                 VRT(1,1),VRT(1,1),VRT(1,1),POP(1,1),
     &                 NFIRST,NREAD,2,IRREPDO,LISTI3)
         CALL XGEMM('T','N',NT2DIS,NREAD,NT2DSZ,ALPHA,ICORE(I010),
     &              NT2DSZ,ICORE(I020),TP3BDZ,FACT,ICORE(IOFZ),NT2DIS)
         IOFZ=IOFZ+NREAD*NT2DIS*IINTFP
         NFIRST=NFIRST+NREAD
         NLEFT =NLEFT-NREAD
         IF(NLEFT.NE.0)GOTO 4
        ENDIF
C
C NOW AUGMENT THE LISTS ON DISK, WRITING THIS
C   CONTRIBUTION AS Z(Ef,Am)
C
        CALL FANCYPUT1(ICORE(I000),ICORE(I010),'FF','1432','N','S',
     &                ONE,VRT(1,1),VRT(1,2),VRT(1,1),POP(1,2),
     &                1,TP3ADS,2,IRREPDO,LSTTAR,RHF)
100    CONTINUE
C
      ELSEIF(SPCASE.EQ.'BABA')THEN
C
C CODE FOR SPIN CASE BABA.
C
       LSTTAR=129
       IF(.NOT.TERM2)GOTO 202
       LISTI1=29
       LISTI2=28
       LISTI3=30
       LISTT1=34
       LISTT2=36
       LISTT3=39
       DO 200 IRREPDO=1,NIRREP
C
C TARGET MATRIX AS USED IN FIRST CONTRACTION
C
        TP1ADZ=IRPDPD(IRREPDO,9)
        TP1ADS=IRPDPD(IRREPDO,20)
C
C I MATRIX AS USED IN FIRST CONTRACTION
C
        TP1BDZ=IRPDPD(IRREPDO,9)
        TP1BDS=IRPDPD(IRREPDO,20)
C
C TARGET MATRIX AS USED IN SECOND CONTRACTION
C
        TP2ADZ=TP1ADZ
        TP2ADS=TP1ADS
C
C I MATRIX AS USED IN SECOND CONTRACTION
C
        TP2BDZ=IRPDPD(IRREPDO,10)
        TP2BDS=IRPDPD(IRREPDO,20)
C
C TARGET MATRIX AS USED IN THIRD CONTRACTION
C
        TP3ADZ=IRPDPD(IRREPDO,13)
        TP3ADS=IRPDPD(IRREPDO,12)
C
C I MATRIX AS USED IN THIRD CONTRACTION
C
        TP3BDZ=IRPDPD(IRREPDO,11)
        TP3BDS=IRPDPD(IRREPDO,13)
C
        MAXTP=ITMPSIZ
C
C
C DO FIRST CONTRACTION
C
C        Z(eFaM) = SUM T(FG,MN) * <Ga|Ne>
C                  N,G
C
C ALPHA SET TO -1 BECAUSE OF T2(AA) AND T2(BB) RING ORDERED STORAGE MODE
C
        ALPHA=ONEM
        BETA=ZILCH
        INTDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTI1))
        INTDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTI1))
        NT2DSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTT1))
        NT2DIS=IRPDPD(IRREPDO,ISYTYP(2,LISTT1))
        I000=1
        I010=I000+IINTFP*TP1ADZ*TP1ADS
        ITMP1=I010+IINTFP*NT2DSZ*NT2DIS
        I020=ITMP1+IINTFP*MAXTP
        I030=I020+IINTFP*MAXTP
        IF(I030.GE.MAXCOR) CALL INSMEM('W5ABIN',I030,MAXCOR)
        MXCOR=MAXCOR-I030+1
        CALL IZERO(ICORE,IINTFP*TP1ADZ*TP1ADS)
C
C READ T2 AMPLITUDES INTO A T2(FM,GN) MATRIX
C
        CALL GETLST(ICORE(I010),1,NT2DIS,1,IRREPDO,LISTT1)
C
        IF(TP1BDZ.NE.0)THEN
         NINCOR=MXCOR/(TP1BDZ*IINTFP)
        ELSE
         NINCOR=TP1BDS
        ENDIF
        NLEFT =TP1BDS
        NFIRST=1
        IOFZ=I000
11      NREAD=MIN(NLEFT,NINCOR)
        CALL FANCYGET(ICORE(I020),ICORE(ITMP1),'FF','1342',
     &                VRT(1,1),VRT(1,2),POP(1,1),VRT(1,2),
     &                NFIRST,NREAD,2,IRREPDO,LISTI1)
        CALL XGEMM('N','N',NT2DSZ,NREAD,NT2DIS,ALPHA,ICORE(I010),
     &             NT2DSZ,ICORE(I020),TP1BDZ,BETA,ICORE(IOFZ),NT2DSZ)
        NFIRST=NFIRST+NREAD
        NLEFT =NLEFT -NREAD
        IOFZ  =IOFZ  +NREAD*NT2DSZ*IINTFP
        IF(NLEFT.NE.0)GOTO 11
C
C READ <Ga|Ne> INTEGRALS INTO AN I(GN,ea) MATRIX
C
c        CALL FANCYGET(ICORE(I010),ICORE(I020),'FF','1342',
c     &                VRT(1,1),VRT(1,2),POP(1,1),VRT(1,2),
c     &                1,TP1BDS,2,IRREPDO,LISTI1)
C
C PERFORM THE MATRIX MULTIPLICATION
C                              
C             Z(FM,ea) = T(FM,GN) * I(GN,ea) 
C
c        CALL XGEMM('N','N',NT2DSZ,TP1BDS,NT2DIS,ALPHA,ICORE(I020),
c     &             NT2DSZ,ICORE(I010),TP1BDZ,BETA,ICORE(I000),NT2DSZ)
C
C
C NOW PERFORM SECOND CONTRACTION
C
C              Z(eFaM) =  - SUM T(Fg|Mn> * <ga||en>
C                           n,g
C
C THE RESULT WILL BE ADDED TO THE PREVIOUS TERM BY ACCUMULATION IN XGEMM
C
        ALPHA=ONEM
        BETA=ONE
        INTDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTI2))
        INTDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTI2))
        NT2DSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTT2))
        NT2DIS=IRPDPD(IRREPDO,ISYTYP(2,LISTT2))
        ITMP1=I010+IINTFP*NT2DSZ*NT2DIS
        I020=ITMP1+IINTFP*MAXTP
        I030=I020+IINTFP*MAXTP
        IF(I030.GE.MAXCOR) CALL INSMEM('W5ABIN',I030,MAXCOR)
        MXCOR=MAXCOR-I030+1
C
C READ T2 AMPLITUDES INTO A T2(gn,FM) MATRIX.
C
        CALL GETLST(ICORE(I010),1,NT2DIS,1,IRREPDO,LISTT2)
C
        IF(TP2BDZ.NE.0)THEN
         NINCOR=MXCOR/(TP2BDZ*IINTFP)
        ELSE
         NINCOR=TP2BDS
        ENDIF
        NLEFT =TP2BDS
        NFIRST=1
        IOFZ=I000
12      NREAD=MIN(NLEFT,NINCOR)
        CALL FANCYGET(ICORE(I020),ICORE(ITMP1),'PF','1432',
     &                VRT(1,2),VRT(1,2),VRT(1,2),POP(1,2),
     &                NFIRST,NREAD,2,IRREPDO,LISTI2)
        CALL XGEMM('T','N',NT2DIS,NREAD,NT2DSZ,ALPHA,ICORE(I010),
     &             NT2DSZ,ICORE(I020),TP2BDZ,BETA,ICORE(IOFZ),NT2DIS)
        NFIRST=NFIRST+NREAD
        NLEFT =NLEFT -NREAD
        IOFZ  =IOFZ  +NREAD*NT2DIS*IINTFP
        IF(NLEFT.NE.0)GOTO 12
C
C READ <ga||en> INTEGRALS INTO AN I(gn,ea) MATRIX.
C
c        CALL FANCYGET(ICORE(I010),ICORE(I020),'PF','1432',
c     &                VRT(1,2),VRT(1,2),VRT(1,2),POP(1,2),
c     &                1,TP2BDS,2,IRREPDO,LISTI2)
C
C NOW PERFORM THE MATRIX MULTIPLICATION
C                             +         
C          Z(FM,ea) = T(gn,FM)  *  I(gn,ea)  
C
c       CALL XGEMM('T','N',NT2DIS,TP2BDS,NT2DSZ,ALPHA,ICORE(I020),
c     &             NT2DSZ,ICORE(I010),TP2BDZ,BETA,ICORE(I000),NT2DIS)
C
C NOW WRITE THIS TERM AND THE PREVIOUS CONTRIBUTION TO DISK
C  AS Z(Fe,Ma)
C
        CALL FANCYPUT(ICORE(I000),ICORE(I010),'FF','1324','N','S',
     &                ONE,VRT(1,1),VRT(1,2),POP(1,1),VRT(1,2),
     &                1,TP2ADS,2,IRREPDO,LSTTAR)
C
C
C NOW PERFORM THIRD CONTRACTION
C
C              Z(eFaM) =-SUM T(Ge,Mn) * <Ga|Fn>             (BABA)
C                        n,G
C
        ALPHA=ONEM
        BETA=ZILCH
        INTDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTI3))
        INTDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTI3))
        NT2DSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTT3))
        NT2DIS=IRPDPD(IRREPDO,ISYTYP(2,LISTT3))
        I000=1
        I010=I000+IINTFP*TP3ADZ*TP3ADS
        ITMP1=I010+IINTFP*NT2DSZ*NT2DIS
        I020=ITMP1+IINTFP*MAXTP
        I030=I020+IINTFP*MAXTP
        IF(I030.GE.MAXCOR) CALL INSMEM('W5ABIN',I030,MAXCOR)
        MXCOR=MAXCOR-I030+1
        CALL IZERO(ICORE,IINTFP*TP3ADZ*TP3ADS)
C
C READ T2 AMPLITUDES INTO A T2(Gn,eM) MATRIX.
C
        CALL GETLST(ICORE(I010),1,NT2DIS,1,IRREPDO,LISTT3)
C
        IF(TP3BDZ.NE.0)THEN
         NINCOR=MXCOR/(TP3BDZ*IINTFP)
        ELSE
         NINCOR=TP3BDS
        ENDIF
        NLEFT =TP3BDS
        NFIRST=1
        IOFZ=I000
13      NREAD=MIN(NLEFT,NINCOR)
        CALL FANCYGET(ICORE(I020),ICORE(ITMP1),'FF','1432',
     &                VRT(1,1),VRT(1,2),VRT(1,1),POP(1,2),
     &                NFIRST,NREAD,2,IRREPDO,LISTI3)
        CALL XGEMM('T','N',NREAD,NT2DIS,NT2DSZ,ALPHA,ICORE(I020),
     &             TP3BDZ,ICORE(I010),NT2DSZ,BETA,ICORE(IOFZ),TP3BDS)
        NFIRST=NFIRST+NREAD
        NLEFT =NLEFT -NREAD
        IOFZ  =IOFZ  +NREAD*IINTFP
        IF(NLEFT.NE.0)GOTO 13
C
C READ <Ga|Fn> INTEGRALS INTO AN I(Gn,Fa) MATRIX.
C
c        CALL FANCYGET(ICORE(I010),ICORE(I020),'FF','1432',
c     &                VRT(1,1),VRT(1,2),VRT(1,1),POP(1,2),
c     &                1,TP3BDS,2,IRREPDO,LISTI3)
C
C NOW FORM THE MATRIX PRODUCT.
C                                     +
C                  Z(Fa,eM) = I(Gn,Fa)  * T(Gn,eM)
C
c        CALL XGEMM('T','N',TP3BDS,NT2DIS,NT2DSZ,ALPHA,ICORE(I010),
c     &             TP3BDZ,ICORE(I020),NT2DSZ,BETA,ICORE(I000),
c     &             TP3BDS)
C
C NOW AUGMENT THE LISTS ON DISK, WRITING THIS AND PREVIOUS TERM
C   AS Z(Fe,Ma)
C
        CALL FANCYPUT(ICORE(I000),ICORE(I010),'FF','1342','N','S',
     &                ONE,VRT(1,1),VRT(1,2),POP(1,1),VRT(1,2),
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
        CALL GETLST(ICORE(I000),1,TARDIS,1,IRREP,LSTTAR-100)
       ELSE
        CALL ZERO(ICORE(I000),TARDIS*TARDSZ)
       ENDIF
       IF(I020.LE.MAXCOR)THEN
        IF(TERM2)THEN
         CALL GETLST(ICORE(I010),1,TARDIS,1,IRREP,LSTTAR)
        ELSE
         CALL ZERO(ICORE(I010),TARDIS*TARDSZ)
        ENDIF
        CALL SAXPY(TARDIS*TARDSZ,ONE,ICORE(I010),1,ICORE(I000),1)
       ELSE
        I020=I010+IINTFP*TARDSZ
        IOFF=I000
        DO 1001 IDIS=1,TARDIS
         IF(TERM2)THEN
          CALL GETLST(ICORE(I010),IDIS,1,1,IRREP,LSTTAR)
         ELSE
          CALL ZERO(ICORE(I010),TARDSZ)
         ENDIF
         CALL SAXPY(TARDSZ,ONE,ICORE(I010),1,ICORE(IOFF),1)
         IOFF=IOFF+TARDSZ*IINTFP
1001    CONTINUE
       ENDIF
       CALL PUTLST(ICORE(I000),1,TARDIS,1,IRREP,LSTTAR)
1000  CONTINUE
C
C WHEW!  WE'RE FINALLY DONE WITH ALL OF THIS AWFUL CRAP!
C
      RETURN
      END
