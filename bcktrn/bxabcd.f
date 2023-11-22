      SUBROUTINE BXABCD(EVEC,BUF1,BUF2,BUF3,BUFTMP,LEN1,LEN2,LEN3,
     &                  LENTMP,IUHF)
C
C THIS ROUTINE TRANSFORMS T2 AMPLITUDES OF SYMMETRY TYPE
C  ABCD FROM THE MO TO THE AO BASIS.  THE AMPLITUDES ARE
C  DEALT WITH IN MULLIKEN ORDER (AI,BJ) AND THE TRANSFORMATION
C  PASSES THROUGH TWO STAGES:
C
C    T(AI,BJ) -> T(XX,BJ) -> T(XX,XX)
C
C EVEC IS THE SYMMETRY PACKED EIGENVECTOR MATRIX AS PREPARED BY SYPKEV.
C
C REQUIRED SPIN CASES:
C             AIBJ
C             ----
C             AAAA (UHF ONLY, COMBINED WITH BBAA)
C             BBBB (UHF ONLY, COMBINED WITh AABB)
C             AABB (RHF AND UHF)
C             BBAA (UHF ONLY, COMBINED WITH AAAA)
C
C MEMORY REQUIREMENTS:
C       
C       BUF1  = MAX(NTAIBJ,NTXXBJ)                        (RHF AND UHF)
C       BUF2  = MAX(NTAIBJ,NXX)                           (RHF AND UHF)
C       BUF3  = NTAIBJ                                    (UHF ONLY)
C       BUFTMP= 2*NAO*NAO                                 (RHF AND UHF)
C 
C  WHERE: 
C        NTAIBJ - LENGTH OF THE IRREP 1 DPD FOR T2.
C        NTXXBJ - MAX[NAO(IR1)*NAO(IR1)*NOC(IR2)*NVRT(IR2)] OVER ALL 
C                 DIRECT PRODUCT PAIRS.
C        NXX    - MAX[NAO(IRR)*NAO(IRR)] FOR ALL IRREPS.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BUF1(LEN1),BUF2(LEN2),BUF3(LEN3),BUFTMP(LENTMP),EVEC(1)
      INTEGER POP,VRT,DIRPRD,DISSIZ,AOPOP,IOFF(8)
      INTEGER IROFF2(8,8),IROFF3(8,8),IDID(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /AOPOPS/ AOPOP(8),MOPOP(8),NAO,NAOSPH,NMO
      COMMON /AOOFST/ INDOCC(8,2),INDVRT(8,2)
C
      DATA ONE   /1.0D0/
      DATA ZILCH /0.0D0/
      DATA FOURTH  /0.25D0/
      DATA FOUR  /4.D0/
      DATA HALF /0.5D0/
C
C LOOP OVER THE DIFFERENT TYPES
C
      DO 1000 ITYPE=1,IUHF+1
C
       LISTW=118
C
       IF(ITYPE.EQ.1) THEN
C
C AABB AND BBBB SPIN CASES
C
        LISTW2=120
        ISL=1
        ISR=2
C
       ELSE
C
C BBAA AND AAAA SPIN CASES
C
        LISTW2=119
        ISL=2
        ISR=1
C
       ENDIF
C
C READ IN SECOND THROUGH LAST IRREPS OF THE SYMMETRY-PACKED T2 VECTOR
C  AND PROCESS EACH ON THE FLY.  ALL ELEMENTS OF THIS CLASS BELONG TO
C  SYMMETRY TYPE ABAB, BAAB AND ABCD.  FROM THESE LISTS, THE ABCD 
C  AMPLITUDES WILL BE EXTRACTED AND ORDERED.
C
       ITHRU=0
       DO 100 IRREPDO=2,NIRREP
C
C READ IN ONE DPD IRREP OF THE T2 VECTOR
C
        NUMDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTW))
        DISSIZ=IRPDPD(IRREPDO,ISYTYP(1,LISTW))
        IF(ITYPE.EQ.1) THEN
         CALL GETLST(BUF1,1,NUMDIS,1,IRREPDO,LISTW)
        ELSE
         CALL GETLST(BUF2,1,NUMDIS,1,IRREPDO,LISTW)
         CALL TRANSP(BUF2,BUF1,NUMDIS,DISSIZ)
         ITMP=NUMDIS
         NUMDIS=DISSIZ
         DISSIZ=ITMP
        ENDIF
C
C COMPUTE OFFSETS FOR THE BEGINNING OF ALL AB SYMMETRY 
C  TYPES WITHIN A LOGICAL RECORD OF T AMPLITUDES.
C
        IPOS=1
        DO 101 IRREPO=1,NIRREP
         IRREPV=DIRPRD(IRREPO,IRREPDO)
         IOFF(IRREPO)=IPOS
         IPOS=IPOS+VRT(IRREPV,ISL)*POP(IRREPO,ISL)
101     CONTINUE
C
C NOW ORDER ELEMENTS AS FOLLOWS:
C
C           X1C X1B X1A 1
C           X1B X1C X1A 1
C           X2C X2B X2A 2
C           X2B X2C X2A 2
C           X3C X3B X3A 3
C           X3B X3C X3A 3
C             .   .   . .
C           XhC XhB XhA h
C           XhB XhC XhA h
C
C WHERE THE SPECIFIC IRREPS ARE GIVEN BY A NUMBER AND
C
C
C    DIRPRD(XiA,i)  = IRREPDO
C    DIRPRD(XiC,XiB)= IRREPDO
C    XiC.NE.XiB.NE.XiA.NE.i
C
C ALSO, WE DON'T WANT TO DO STUFF LIKE 3 4 1 2 AND THEN 1 2 3 4
C  SINCE THESE ARE THE SAME FOR EVERYTHING EXCEPT UHF AB.
C
        IOFF2=1
        IOFFRCB=0
        DO 110 IRREPI=1,NIRREP
         IRREPXIA=DIRPRD(IRREPI,IRREPDO)
         IBOT=MAX(IRREPI,IRREPXIA)+1
         CALL IZERO(IDID,NIRREP)
         DO 111 IRREPTMP=IBOT,NIRREP
          IRREPXIC=DIRPRD(IRREPTMP,IRREPDO)
          IF(IRREPI.NE.IRREPTMP.AND.IRREPI.NE.IRREPXIC.AND.
     &       IDID(IRREPTMP).EQ.0)THEN
           IRREPXIB=MIN(IRREPTMP,IRREPXIC)
           IRREPXIC=MAX(IRREPTMP,IRREPXIC)
           IDID(IRREPXIB)=1
           IDID(IRREPXIC)=1
           IOFFRC=IOFFRCB
           NSIZVR=VRT(IRREPXIA,ISR)
           NSIZOR=POP(IRREPI,ISR)
           IROFF2(IRREPXIB,IRREPI)=IOFF2
           NSIZVL=VRT(IRREPXIC,ISL)
           NSIZOL=POP(IRREPXIB,ISL)
           LEN=NSIZVL*NSIZOL
           DO 112 INDEXJ=1,NSIZOR
            DO 113 INDEXB=1,NSIZVR
C
C COPY IN THE (XiC,XiB) PIECE FROM THIS RECORD.
C
             IOFF1=IOFFRC+IOFF(IRREPXIB)
             CALL SCOPY(LEN,BUF1(IOFF1),1,BUF2(IOFF2),1)
             IOFF2=IOFF2+LEN
C
             IOFFRC=IOFFRC+DISSIZ
C
113         CONTINUE
112        CONTINUE
           IROFF2(IRREPXIC,IRREPI)=IOFF2
           NSIZVL=VRT(IRREPXIB,ISL)
           NSIZOL=POP(IRREPXIC,ISL)
           LEN=NSIZVL*NSIZOL
           IOFFRC=IOFFRCB
           DO 212 INDEXJ=1,NSIZOR
            DO 213 INDEXB=1,NSIZVR
C
C COPY IN THE (XiB,XiC) PIECE FROM THIS RECORD.
C
             IOFF1=IOFFRC+IOFF(IRREPXIC)
             CALL SCOPY(LEN,BUF1(IOFF1),1,BUF2(IOFF2),1)
             IOFF2=IOFF2+LEN
C
             IOFFRC=IOFFRC+DISSIZ
C
213         CONTINUE
212        CONTINUE
          ENDIF
111      CONTINUE
         IOFFRCB=IOFFRCB+DISSIZ*VRT(IRREPXIA,ISR)*POP(IRREPI,ISR)
110     CONTINUE
C
C FOR UHF TREAT HERE AS WELL THE AAAA AND BBBB LIST, SCALE ALSO
C THE AABB AND BBAA LIST WITh A FACTOR OF 0.5 ( SAME ARGUMENT AS IN
C BXAABB AND BXABAB)
C
        IF(IUHF.EQ.1) THEN
C
         CALL SSCAL(IOFF2-1,HALF,BUF2,1)
C
         DISSIZ=IRPDPD(IRREPDO,ISYTYP(1,LISTW2))
         NUMIDS=IRPDPD(IRREPDO,ISYTYP(2,LISTW2))
         CALL GETLST(BUF1,1,NUMDIS,1,IRREPDO,LISTW2)
C
C COMPUTE OFFSETS FOR THE BEGINNING OF ALL AB SYMMETRY 
C  TYPES WITHIN A LOGICAL RECORD OF T AMPLITUDES.
C
         IPOS=1
         DO 201 IRREPO=1,NIRREP
          IRREPV=DIRPRD(IRREPO,IRREPDO)
          IOFF(IRREPO)=IPOS
          IPOS=IPOS+VRT(IRREPV,ISR)*POP(IRREPO,ISR)
201      CONTINUE
C
C NOW ORDER ELEMENTS AS FOLLOWS:
C
C           X1C X1B X1A 1
C           X1B X1C X1A 1
C           X2C X2B X2A 2
C           X2B X2C X2A 2
C           X3C X3B X3A 3
C           X3B X3C X3A 3
C             .   .   . .
C           XhC XhB XhA h
C           XhB XhC XhA h
C
C WHERE THE SPECIFIC IRREPS ARE GIVEN BY A NUMBER AND
C
C
C    DIRPRD(XiA,i)  = IRREPDO
C    DIRPRD(XiC,XiB)= IRREPDO
C    XiC.NE.XiB.NE.XiA.NE.i
C
C ALSO, WE DON'T WANT TO DO STUFF LIKE 3 4 1 2 AND THEN 1 2 3 4
C  SINCE THESE ARE THE SAME FOR EVERYTHING EXCEPT UHF AB.
C
         IOFF2=1
         IOFFRCB=0
         DO 210 IRREPI=1,NIRREP
          IRREPXIA=DIRPRD(IRREPI,IRREPDO)
          IBOT=MAX(IRREPI,IRREPXIA)+1
          CALL IZERO(IDID,NIRREP)
          DO 211 IRREPTMP=IBOT,NIRREP
           IRREPXIC=DIRPRD(IRREPTMP,IRREPDO)
           IF(IRREPI.NE.IRREPTMP.AND.IRREPI.NE.IRREPXIC.AND.
     &        IDID(IRREPTMP).EQ.0)THEN
            IRREPXIB=MIN(IRREPTMP,IRREPXIC)
            IRREPXIC=MAX(IRREPTMP,IRREPXIC)
            IDID(IRREPXIB)=1
            IDID(IRREPXIC)=1
            IOFFRC=IOFFRCB
            NSIZVR=VRT(IRREPXIA,ISR)
            NSIZOR=POP(IRREPI,ISR)
            IROFF3(IRREPXIB,IRREPI)=IOFF2
            NSIZVL=VRT(IRREPXIC,ISR)
            NSIZOL=POP(IRREPXIB,ISR)
            LEN=NSIZVL*NSIZOL
            DO 312 INDEXJ=1,NSIZOR
             DO 313 INDEXB=1,NSIZVR
C
C COPY IN THE (XiC,XiB) PIECE FROM THIS RECORD.
C
              IOFF1=IOFFRC+IOFF(IRREPXIB)
              CALL SCOPY(LEN,BUF1(IOFF1),1,BUF3(IOFF2),1)
              IOFF2=IOFF2+LEN
C
              IOFFRC=IOFFRC+DISSIZ
C
313          CONTINUE
312         CONTINUE
            IROFF3(IRREPXIC,IRREPI)=IOFF2
            NSIZVL=VRT(IRREPXIB,ISR)
            NSIZOL=POP(IRREPXIC,ISR)
            LEN=NSIZVL*NSIZOL
            IOFFRC=IOFFRCB
            DO 412 INDEXJ=1,NSIZOR
             DO 413 INDEXB=1,NSIZVR
C
C COPY IN THE (XiB,XiC) PIECE FROM THIS RECORD.
C
              IOFF1=IOFFRC+IOFF(IRREPXIC)
              CALL SCOPY(LEN,BUF1(IOFF1),1,BUF3(IOFF2),1)
              IOFF2=IOFF2+LEN
C 
              IOFFRC=IOFFRC+DISSIZ
C
413          CONTINUE
412         CONTINUE
           ENDIF
211       CONTINUE
          IOFFRCB=IOFFRCB+DISSIZ*VRT(IRREPXIA,ISR)*POP(IRREPI,ISR)
210      CONTINUE
C
        ENDIF
C 
C DRIVE TRANSFORMATION BY LOOPING OVER APPROPRIATE IRREP COMBINATIONS.
C
C KEEP IRREPS CANONICAL!
C 
       DO 120 IRREPI=1,NIRREP
        IRREPXIA=DIRPRD(IRREPI,IRREPDO)
        IF(IRREPXIA.LT.IRREPI)GOTO 120
        IBOT=MAX(IRREPI,IRREPXIA)+1
        CALL IZERO(IDID,8)
        DO 130 ITMP=IBOT,NIRREP
         CALL ZERO(BUF1,LEN1)
         IRREPXIC=DIRPRD(ITMP,IRREPDO)
         ITMP2=MIN(IRREPXIC,ITMP)
         IRREPXIC=MAX(ITMP,IRREPXIC)
         IRREPXIB=ITMP2
         IF(MAX(IDID(IRREPXIC),IDID(IRREPXIB)).NE.0)GOTO 130
         IDID(IRREPXIC)=1
         IDID(IRREPXIB)=1
         ITHRU=ITHRU+1
C
C NOW THEY ARE CANONICAL AND UNIQUE.  PROCEED WITH THE TRANSFORMATION.
C
         NAOXIC=AOPOP(IRREPXIC)
         NAOXIB=AOPOP(IRREPXIB)
         NAOXIA=AOPOP(IRREPXIA)
         NAOI  =AOPOP(IRREPI)
         NSIZFL=NAOXIC*NAOXIB
         NSIZFR=NAOXIA*NAOI
C
         NAXIC =VRT(IRREPXIC,ISL)
         NAXIB =VRT(IRREPXIB,ISL)
         NAXIA =VRT(IRREPXIA,ISL)
         NAI   =VRT(IRREPI,ISL)
         NIXIC =POP(IRREPXIC,ISL)
         NIXIB =POP(IRREPXIB,ISL)
         NIXIA =POP(IRREPXIA,ISL)
         NII   =POP(IRREPI,ISL)
         NBXIC =VRT(IRREPXIC,ISR)
         NBXIB =VRT(IRREPXIB,ISR)
         NBXIA =VRT(IRREPXIA,ISR)
         NBI   =VRT(IRREPI,ISR)
         NJXIC =POP(IRREPXIC,ISR)
         NJXIB =POP(IRREPXIB,ISR)
         NJXIA =POP(IRREPXIA,ISR)
         NJI   =POP(IRREPI,ISR)
C
         NSZAI1=NAXIC*NIXIB
         NSZBJ1=NBXIA*NJI
         NSZAI2=NAXIB*NIXIC
         NSZBJ2=NBI*NJXIA
C
         IOFF1  =1
         IOFFB  =NSIZFL+1
         IOFF2A1=IROFF2(IRREPXIB,IRREPI)
         IOFF2B1=IROFF2(IRREPXIC,IRREPI)
         IOFF2A2=IROFF2(IRREPXIB,IRREPXIA)
         IOFF2B2=IROFF2(IRREPXIC,IRREPXIA)
         IOFFEVA=INDVRT(IRREPXIC,ISL)
         IOFFEVB=INDVRT(IRREPXIB,ISL)
         IOFFEOA=INDOCC(IRREPXIC,ISL)
         IOFFEOB=INDOCC(IRREPXIB,ISL)
         IBLKSZA=NAXIC*NIXIB
         IBLKSZB=NAXIB*NIXIC
C
C DO THE SYMMETRIZED FIRST HALF TRANSFORMATION:
C
C                                          +  
CT2(XX,BxiaJi) =  Cxic T2(AxicIib,BxiaJi) C xib
C
C                                             +    + 
C               + [Cxib T2(AxibIxic,Bxia,Ji) C xic] 
C
         DO 140 INDEXJ=1,NJI
          DO 150 INDEXB=1,NBXIA
           CALL SYTRAB(BUF2(IOFF2A1),BUF2(IOFF2B1),EVEC(IOFFEVA),
     &                 EVEC(IOFFEVB),EVEC(IOFFEOA),EVEC(IOFFEOB),
     &                 NAOXIC,NAOXIB,NAXIC,NAXIB,NIXIC,NIXIB,BUFTMP,
     &                 BUF1(IOFF1),0)
           IOFF2A1=IOFF2A1+IBLKSZA
           IOFF2B1=IOFF2B1+IBLKSZB
           IOFF1=IOFF1+NSIZFL
150       CONTINUE
140      CONTINUE
C
C NOW DO THE SYMMETRIZED FIRST HALF TRANSFORMATION:
C
C                                            +   
C T2(XX,BxiaJi) =  Cxic T2(AxicIxib,BiJxia) C xib 
C
C                                           +    +
C               + [Cxib T2(AibIxic,BiJxia) C xic]
C
         ILCIXI=IOFF1
         DO 160 INDEXJ=1,NJXIA
          DO 170 INDEXB=1,NBI
           CALL SYTRAB(BUF2(IOFF2A2),BUF2(IOFF2B2),EVEC(IOFFEVA),
     &                 EVEC(IOFFEVB),EVEC(IOFFEOA),EVEC(IOFFEOB),
     &                 NAOXIC,NAOXIB,NAXIC,NAXIB,NIXIC,NIXIB,BUFTMP,
     &                 BUF1(IOFF1),0)
           IOFF2A2=IOFF2A2+IBLKSZA
           IOFF2B2=IOFF2B2+IBLKSZB
           IOFF1=IOFF1+NSIZFL
170       CONTINUE
160      CONTINUE
C
C FOR UHF SDOME MORE WORK HAS TO BE DONE 
C
         IF(IUHF.EQ.1) THEN
C
         NAXIC =VRT(IRREPXIC,ISR)
         NAXIB =VRT(IRREPXIB,ISR)
         NAXIA =VRT(IRREPXIA,ISR)
         NAI   =VRT(IRREPI,ISR)
         NIXIC =POP(IRREPXIC,ISR)
         NIXIB =POP(IRREPXIB,ISR)
         NIXIA =POP(IRREPXIA,ISR)
         NII   =POP(IRREPI,ISR)
C
         NSZAI1=NAXIC*NIXIB
         NSZAI2=NAXIB*NIXIC
C
         IOFF1  =1
         IOFFB  =NSIZFL+1
         IOFF2A1=IROFF3(IRREPXIB,IRREPI)
         IOFF2B1=IROFF3(IRREPXIC,IRREPI)
         IOFF2A2=IROFF3(IRREPXIB,IRREPXIA)
         IOFF2B2=IROFF3(IRREPXIC,IRREPXIA)
         IOFFEVA=INDVRT(IRREPXIC,ISR)
         IOFFEVB=INDVRT(IRREPXIB,ISR)
         IOFFEOA=INDOCC(IRREPXIC,ISR)
         IOFFEOB=INDOCC(IRREPXIB,ISR)
         IBLKSZA=NAXIC*NIXIB
         IBLKSZB=NAXIB*NIXIC
C
C DO THE SYMMETRIZED FIRST HALF TRANSFORMATION:
C
C                                          +  
CT2(XX,BxiaJi) =  Cxic T2(AxicIib,BxiaJi) C xib
C
C                                             +    + 
C               + [Cxib T2(AxibIxic,Bxia,Ji) C xic] 
C
         DO 240 INDEXJ=1,NJI
          DO 250 INDEXB=1,NBXIA
           CALL SYTRAB(BUF3(IOFF2A1),BUF3(IOFF2B1),EVEC(IOFFEVA),
     &                 EVEC(IOFFEVB),EVEC(IOFFEOA),EVEC(IOFFEOB),
     &                 NAOXIC,NAOXIB,NAXIC,NAXIB,NIXIC,NIXIB,BUFTMP,
     &                 BUF1(IOFF1),1)
           IOFF2A1=IOFF2A1+IBLKSZA
           IOFF2B1=IOFF2B1+IBLKSZB
           IOFF1=IOFF1+NSIZFL
250       CONTINUE
240      CONTINUE
C
C NOW DO THE SYMMETRIZED FIRST HALF TRANSFORMATION:
C
C                                            +   
C T2(XX,BxiaJi) =  Cxic T2(AxicIxib,BiJxia) C xib 
C
C                                           +    +
C               + [Cxib T2(AibIxic,BiJxia) C xic]
C
         ILCIXI=IOFF1
         DO 260 INDEXJ=1,NJXIA
          DO 270 INDEXB=1,NBI
           CALL SYTRAB(BUF3(IOFF2A2),BUF3(IOFF2B2),EVEC(IOFFEVA),
     &                 EVEC(IOFFEVB),EVEC(IOFFEOA),EVEC(IOFFEOB),
     &                 NAOXIC,NAOXIB,NAXIC,NAXIB,NIXIC,NIXIB,BUFTMP,
     &                 BUF1(IOFF1),1)
           IOFF2A2=IOFF2A2+IBLKSZA
           IOFF2B2=IOFF2B2+IBLKSZB
           IOFF1=IOFF1+NSIZFL
270       CONTINUE
260      CONTINUE
C
       ENDIF
C
C
C 
C NOW BOTH THE T2(XX,BxiaJi) AND T2(XX,BiJxia) HALF TRANSFORMED 
C  AMPLITUDES ARE ON BUF1.  PICK THESE OFF AND DO THE SYMMETRIZED 
C  SECOND HALF TRANSFORMS.
C
C DO THE SYMMETRIZED HALF TRANSFORMATION
C
C                                       +                       +
C       T2(XX,XX) = Cxia T2(XX,BxiaJi) C i + [Ci T2(XX,BiJxia) C xia]
C
         IF(ITYPE.EQ.1) THEN
         IOFFEVA=INDVRT(IRREPXIA,ISR)
         IOFFEVB=INDVRT(IRREPI  ,ISR)
         IOFFEOA=INDOCC(IRREPXIA,ISR)
         IOFFEOB=INDOCC(IRREPI  ,ISR)
         IOFF1A=1
         IOFF1B=ILCIXI
         IOFFBA=NSZBJ1+1
         IOFFBB=NSZBJ2+IOFFBA
         IOFFCC=NAOXIA*NAOI+IOFFBB
C
         DO 180 INDXNU=1,NAOXIB
          DO 190 INDXMU=1,NAOXIC
           CALL SCOPY(NSZBJ1,BUF1(IOFF1A),NSIZFL,BUFTMP,1)
           CALL SCOPY(NSZBJ2,BUF1(IOFF1B),NSIZFL,
     &                BUFTMP(IOFFBA),1)
           CALL SYTRAB(BUFTMP,BUFTMP(IOFFBA),EVEC(IOFFEVA),
     &                 EVEC(IOFFEVB),EVEC(IOFFEOA),EVEC(IOFFEOB),
     &                 NAOXIA,NAOI,NAXIA,NAI,NIXIA,NII,
     &                 BUFTMP(IOFFBB),BUFTMP(IOFFCC),0)
           CALL SSCAL(NAOXIA*NAOI,FOURTH,BUFTMP(IOFFCC),1)
           CALL PUTLST(BUFTMP(IOFFCC),IOFF1A,1,1,4,ITHRU)
           IOFF1A=IOFF1A+1
           IOFF1B=IOFF1B+1
190       CONTINUE
180      CONTINUE
C
       ELSE
         IOFFEVA=INDVRT(IRREPXIA,ISR)
         IOFFEVB=INDVRT(IRREPI  ,ISR)
         IOFFEOA=INDOCC(IRREPXIA,ISR)
         IOFFEOB=INDOCC(IRREPI  ,ISR)
         IOFF1A=1
         IOFF1B=ILCIXI
         IOFFBA=NSZBJ1+1
         IOFFBB=NSZBJ2+IOFFBA
         IOFFCC=NAOXIA*NAOI+IOFFBB
C
         DO 280 INDXNU=1,NAOXIB
          DO 290 INDXMU=1,NAOXIC
           CALL SCOPY(NSZBJ1,BUF1(IOFF1A),NSIZFL,BUFTMP,1)
           CALL SCOPY(NSZBJ2,BUF1(IOFF1B),NSIZFL,
     &                BUFTMP(IOFFBA),1)
           CALL GETLST(BUFTMP(IOFFCC),IOFF1A,1,1,4,ITHRU)
           CALL SSCAL(NAOXIA*NAOI,FOUR,BUFTMP(IOFFCC),1)       
           CALL SYTRAB(BUFTMP,BUFTMP(IOFFBA),EVEC(IOFFEVA),
     &                 EVEC(IOFFEVB),EVEC(IOFFEOA),EVEC(IOFFEOB),
     &                 NAOXIA,NAOI,NAXIA,NAI,NIXIA,NII,
     &                 BUFTMP(IOFFBB),BUFTMP(IOFFCC),1)
           CALL SSCAL(NAOXIA*NAOI,FOURTH,BUFTMP(IOFFCC),1)
           CALL PUTLST(BUFTMP(IOFFCC),IOFF1A,1,1,4,ITHRU)
           IOFF1A=IOFF1A+1
           IOFF1B=IOFF1B+1
290       CONTINUE
280      CONTINUE
C
          ENDIF
130     CONTINUE
120    CONTINUE
100   CONTINUE
C
1000  CONTINUE
      RETURN
      END
