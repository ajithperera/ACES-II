      SUBROUTINE BXABAB1(EVEC,BUF1,BUF2,BUF3,BUFTMP,L1MAX,LEN2,LEN3,
     &                   LENTMP,IUHF)
C
C THIS ROUTINE TRANSFORMS T2 AMPLITUDES OF SYMMETRY TYPE
C  ABAB FROM THE MO TO THE AO BASIS.  THE AMPLITUDES ARE
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
C             BBBB (UHF ONLY, COMBINED WITH AABB)
C             AABB (RHF AND UHF)
C             BBAA (UHF ONLY, COMBINED WITH AAAA)
C
C MEMORY REQUIREMENTS:
C       
C       BUF1  = MAX(NTAIBJ,NTXXBJ)                  (RHF AND UHF)
C       BUF2  = MAX(NTAIBJ,NXX)                     (RHF AND UHF)
C       BUF3  = NTAIBJ                              (UHF ONLY)
C       BUFTMP= 2*NAO*NAO                           (RHF AND UHF)
C 
C  WHERE: 
C        NTAIBJ - LENGTH OF THE IRREP 1 DPD FOR T2.
C        NTXXBJ - MAX[NAO(IR1)*NAO(IR1)*NOC(IR2)*NVRT(IR2)] OVER ALL 
C                 DIRECT PRODUCT PAIRS.
C        NXX    - MAX[NAO(IRR)*NAO(IRR)] FOR ALL IRREPS.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BUF1(L1MAX),BUF2(LEN2),BUF3(LEN3),BUFTMP(LENTMP),EVEC(1)
      INTEGER POP,VRT,DIRPRD,DISSIZ,AOPOP,ISTART2(8),ISTART3(8),IOFF(8)
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
      DATA FOURTH/0.25D0/
      DATA HALF/0.5D0/
      DATA FOUR /4.0D0/
C
C  LOOP OVER THE DIFFERENT TYPES
C
      DO 1000 ITYPE=1,IUHF+1
C
       LISTW=118
C
       IF(ITYPE.EQ.1) THEN
C
C AABB ABD BBBB SPIN CASES
C
        LISTW2=120
        ISL=1
        ISR=2
C
       ELSE
C
        LISTW2=119
        ISL=2
        ISR=1
C
       ENDIF
C
C READ IN SECOND THROUGH LAST IRREPS OF THE SYMMETRY-PACKED T2 VECTOR
C  AND PROCESS EACH ON THE FLY.  ALL ELEMENTS OF THIS CLASS BELONG TO
C  SYMMETRY TYPE ABAB, BAAB AND ABCD.  FROM THESE LISTS, THE ABAB AND 
C  BAAB AMPLITUDES WILL BE EXTRACTED AND TRANSFORMED.
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
C           X1  1 X1  1
C            1 X1 X1  1
C           X2  2 X2  2
C            2 X2 X2  2
C            .  .  .  .
C           Xh  h Xh  h
C            h Xh Xh  h
C
C WHERE THE SPECIFIC IRREPS ARE GIVEN BY A NUMBER AND
C    DIRPRD(Xi,i) = IRREPDO
C
        IOFF2=1
        IOFFRC=0
        DO 110 IRREPO=1,NIRREP
         IRREPV=DIRPRD(IRREPO,IRREPDO)
         NSIZVL=VRT(IRREPV,ISL)
         NSIZVR=VRT(IRREPV,ISR)
         NSIZOL=POP(IRREPO,ISL)
         NSIZOR=POP(IRREPO,ISR) 
         NSIZEL=NSIZVL*NSIZOL
         NSIZER=NSIZVR*NSIZOR
         ISTART2(IRREPO)=IOFF2
         LEN11=VRT(IRREPV,ISL)*POP(IRREPO,ISL)
         LEN22=VRT(IRREPO,ISL)*POP(IRREPV,ISL)
         DO 111 INDEXJ=1,NSIZOR
          DO 112 INDEXB=1,NSIZVR
C
C COPY IN THE (Xi,i) PIECE FIRST, FOLLOWED BY (i,Xi) PIECE.
C
           IOFF1=IOFFRC+IOFF(IRREPO)
           CALL SCOPY(LEN11,BUF1(IOFF1),1,BUF2(IOFF2),1)
           IOFF2=IOFF2+LEN11
C
           IOFF1=IOFFRC+IOFF(IRREPV)
           CALL SCOPY(LEN22,BUF1(IOFF1),1,BUF2(IOFF2),1)
           IOFF2=IOFF2+LEN22
C
           IOFFRC=IOFFRC+DISSIZ
C
112       CONTINUE
111      CONTINUE
110     CONTINUE
C
C FOR UHF TREAT HERE AS WELL THE AAAA AND BBBB LIST, SCALE ALSO
C THE AABB AND BBAA LIST WITh A FACTOR OF 0.5 ( SAME ARGUMENT AS IN
C BXAABB)
        IF(IUHF.EQ.1) THEN
C
         CALL SSCAL(IOFF2-1,HALF,BUF2,1) 
C
         DISSIZ=IRPDPD(IRREPDO,ISYTYP(1,LISTW2))
         NUMDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTW2))
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
C           X1  1 X1  1
C            1 X1 X1  1
C           X2  2 X2  2
C            2 X2 X2  2
C            .  .  .  .
C           Xh  h Xh  h
C            h Xh Xh  h
C
C WHERE THE SPECIFIC IRREPS ARE GIVEN BY A NUMBER AND
C    DIRPRD(Xi,i) = IRREPDO
C
         IOFF2=1
         IOFFRC=0
         DO 210 IRREPO=1,NIRREP
          IRREPV=DIRPRD(IRREPO,IRREPDO)
          NSIZVL=VRT(IRREPV,ISR)
          NSIZVR=VRT(IRREPV,ISR)
          NSIZOL=POP(IRREPO,ISR)
          NSIZOR=POP(IRREPO,ISR) 
          NSIZEL=NSIZVL*NSIZOL
          NSIZER=NSIZVR*NSIZOR
          ISTART3(IRREPO)=IOFF2
          LEN11=VRT(IRREPV,ISR)*POP(IRREPO,ISR)
          LEN22=VRT(IRREPO,ISR)*POP(IRREPV,ISR)
          DO 211 INDEXJ=1,NSIZOR
           DO 212 INDEXB=1,NSIZVR
C 
C COPY IN THE (Xi,i) PIECE FIRST, FOLLOWED BY (i,Xi) PIECE.
C
            IOFF1=IOFFRC+IOFF(IRREPO)
            CALL SCOPY(LEN11,BUF1(IOFF1),1,BUF3(IOFF2),1)
            IOFF2=IOFF2+LEN11
C
            IOFF1=IOFFRC+IOFF(IRREPV)
            CALL SCOPY(LEN22,BUF1(IOFF1),1,BUF3(IOFF2),1)
            IOFF2=IOFF2+LEN22
C
            IOFFRC=IOFFRC+DISSIZ
C
212        CONTINUE
211       CONTINUE
210      CONTINUE
C
        ENDIF
C 
C DRIVE TRANSFORMATION BY LOOPING OVER IRREP PAIRS. SKIP LOOP IF
C  IRREPXI < IRREPI 
C 
        DO 120 IRREPI=1,NIRREP
         IRREPXI=DIRPRD(IRREPI,IRREPDO)
         IF(IRREPXI.LT.IRREPI)GOTO 120
         ITHRU=ITHRU+1
         NAOXI =AOPOP(IRREPXI)
         NAOI  =AOPOP(IRREPI)
         NSIZF =NAOXI*NAOI
C
         NAXI  =VRT(IRREPXI,ISL)
         NAI   =VRT(IRREPI,ISL)
         NIXI  =POP(IRREPXI,ISL)
         NII   =POP(IRREPI,ISL)
         NBXI  =VRT(IRREPXI,ISR)
         NBI   =VRT(IRREPI,ISR)
         NJXI  =POP(IRREPXI,ISR)
         NJI   =POP(IRREPI,ISR)
C
         NSZAI1=NAXI*NII
         NSZBJ1=NBXI*NJI
         NSZAI2=NAI*NIXI
         NSZBJ2=NBI*NJXI
C
C
C  FOR OUT OF CORE 
C
         ISTART=1
C 
         LEN1=(NJI*NBXI+NBI*NJXI)*NAOXI
C 
         IF(LEN1.NE.0) THEN
          MAXAO=MIN((L1MAX/LEN1),NAOI)
         ELSE
          MAXAO=NAOI
         ENDIF
C
C REENTRY POINT FOR OUT OF CORE
C
10       CONTINUE         
C
         CALL ZERO(BUF1,LEN1*MAXAO)
         IEND=MIN(NAOI,ISTART+MAXAO-1)
         NAOI1=IEND+1-ISTART
         NSIZF1=NAOXI*NAOI1
C
         IOFF1  =1
         IOFF2A =ISTART2(IRREPI)
         IOFF2B =ISTART2(IRREPXI)
         IOFFEVA=INDVRT(IRREPXI,ISL)
         IOFFEVB1=INDVRT(IRREPI,ISL)+ISTART-1
         IOFFEVB2=INDVRT(IRREPI,ISL)
         IOFFEOA=INDOCC(IRREPXI,ISL)
         IOFFEOB1=INDOCC(IRREPI,ISL)+ISTART-1
         IOFFEOB2=INDOCC(IRREPI,ISL)
         IOFFG1A=IOFF2A
         IOFFG2A=IOFFG1A+NAXI*NII
         IOFFG1B=IOFF2B
         IOFFG2B=IOFFG1B+NAI*NIXI
         IBLKSIZ=NAXI*NII+NAI*NIXI
C
C DO THE SYMMETRIZED FIRST HALF TRANSFORMATION:
C
C                                      +                         +   +
C T2(XX,BxiJi) =  Cxi T2(AxiIi,BxiJi) C i + [Ci T2(AiIxi,BxiJi) C xi]
C
         DO 130 INDEXJ=1,NJI
          DO 140 INDEXB=1,NBXI
c           call checksum('sytrab2',buf2(ioffg1a),naxi*nii)
           CALL SYTRAB2(BUF2(IOFFG1A),BUF2(IOFFG2A),EVEC(IOFFEVA),
     &                  EVEC(IOFFEVB1),EVEC(IOFFEOA),EVEC(IOFFEOB1),
     &                  NAOXI,NAOI1,NAOI,NAXI,NAI,NIXI,NII,BUFTMP,
     &                  BUF1(IOFF1),0)
c           write(*,*) 'incdices',indexj,indexb
c           write(*,1234) (buf1(ioff1+i-1),i=1,nsizf1)
           IOFFG1A=IOFFG1A+IBLKSIZ
           IOFFG2A=IOFFG2A+IBLKSIZ
           IOFF1=IOFF1+NSIZF1
140       CONTINUE
130      CONTINUE
C
C NOW DO THE SYMMETRIZED FIRST HALF TRANSFORMATION:
C
C                                      +                         +   +
C T2(XX,BxiJi) =  Cxi T2(AxiIi,BiJxi) C i + [Ci T2(AiIxi,BiJxi) C xi]
C
         ILCIXI=IOFF1
         DO 150 INDEXJ=1,NJXI
          DO 160 INDEXB=1,NBI
           CALL SYTRAB2(BUF2(IOFFG2B),BUF2(IOFFG1B),EVEC(IOFFEVA),
     &                  EVEC(IOFFEVB1),EVEC(IOFFEOA),EVEC(IOFFEOB1),
     &                  NAOXI,NAOI1,NAOI,NAXI,NAI,NIXI,NII,BUFTMP,
     &                  BUF1(IOFF1),0)
c           write(*,*) 'incdices',indexj,indexb
c           write(*,1234) (buf1(ioff1+i-1),i=1,nsizf1)
           IOFFG1B=IOFFG1B+IBLKSIZ
           IOFFG2B=IOFFG2B+IBLKSIZ
           IOFF1=IOFF1+NSIZF1
160       CONTINUE
150      CONTINUE
C
C  FOR UHF SOME MORE WORK HAS TO BE DONE
C
         IF(IUHF.EQ.1) THEN
          NAXI  =VRT(IRREPXI,ISR)
          NAI   =VRT(IRREPI,ISR)
          NIXI  =POP(IRREPXI,ISR)
          NII   =POP(IRREPI,ISR)
C
          NSZAI1=NAXI*NII
          NSZAI2=NAI*NIXI
C
          IOFF1  =1
          IOFF3A =ISTART3(IRREPI)
          IOFF3B =ISTART3(IRREPXI)
          IOFFEVA=INDVRT(IRREPXI,ISR)
          IOFFEVB1=INDVRT(IRREPI,ISR)+ISTART-1
          IOFFEVB2=INDVRT(IRREPI,ISR)
          IOFFEOA=INDOCC(IRREPXI,ISR)
          IOFFEOB1=INDOCC(IRREPI,ISR)+ISTART-1
          IOFFEOB2=INDOCC(IRREPI,ISR)
          IOFFG1A=IOFF3A
          IOFFG2A=IOFFG1A+NAXI*NII
          IOFFG1B=IOFF3B
          IOFFG2B=IOFFG1B+NAI*NIXI
          IBLKSIZ=NAXI*NII+NAI*NIXI
C
C DO THE SYMMETRIZED FIRST HALF TRANSFORMATION:
C
C                                      +                         +   +
C T2(XX,BxiJi) =  Cxi T2(AxiIi,BxiJi) C i + [Ci T2(AiIxi,BxiJi) C xi]
C
          DO 230 INDEXJ=1,NJI
           DO 240 INDEXB=1,NBXI
            CALL SYTRAB2(BUF3(IOFFG1A),BUF3(IOFFG2A),EVEC(IOFFEVA),
     &                   EVEC(IOFFEVB1),EVEC(IOFFEOA),EVEC(IOFFEOB1),
     &                   NAOXI,NAOI1,NAOI,NAXI,NAI,NIXI,NII,BUFTMP,
     &                   BUF1(IOFF1),1)
            IOFFG1A=IOFFG1A+IBLKSIZ
            IOFFG2A=IOFFG2A+IBLKSIZ
            IOFF1=IOFF1+NSIZF1
240       CONTINUE
230      CONTINUE
C
C NOW DO THE SYMMETRIZED FIRST HALF TRANSFORMATION:
C
C                                      +                         +   +
C T2(XX,BxiJi) =  Cxi T2(AxiIi,BiJxi) C i + [Ci T2(AiIxi,BiJxi) C xi]
C
          ILCIXI=IOFF1
          DO 250 INDEXJ=1,NJXI
           DO 260 INDEXB=1,NBI
            CALL SYTRAB2(BUF3(IOFFG2B),BUF3(IOFFG1B),EVEC(IOFFEVA),
     &                   EVEC(IOFFEVB1),EVEC(IOFFEOA),EVEC(IOFFEOB1),
     &                   NAOXI,NAOI1,NAOI,NAXI,NAI,NIXI,NII,BUFTMP,
     &                   BUF1(IOFF1),1)
            IOFFG1B=IOFFG1B+IBLKSIZ
            IOFFG2B=IOFFG2B+IBLKSIZ
            IOFF1=IOFF1+NSIZF1
260        CONTINUE
250       CONTINUE
C
         ENDIF
C
C 
C NOW BOTH THE T2(XX,BxiJi) AND T2(XX,BiJxi) HALF TRANSFORMED 
C  AMPLITUDES ARE ON BUF1.  PICK THESE OFF AND DO THE SYMMETRIZED 
C  SECOND HALF TRANSFORMS.
C
C DO THE SYMMETRIZED HALF TRANSFORMATION
C
C                                     +                      +
C       T2(XX,XX) = Cxi T2(XX,BxiJi) C i + [Ci T2(XX,BiJxi) C xi]
C
         IF(ITYPE.EQ.1) THEN 
C
          IOFF1A=1
          IOFF1AA=1+(ISTART-1)*NAOXI
          IOFF1B=ILCIXI
          IOFFBA=NSZBJ1+1
          IOFFBB=NSZBJ2+IOFFBA
          IOFFCC=IOFFBB+NAOXI*NAOI
C
          DO 170 INDXNU=ISTART,IEND
           DO 180 INDXMU=1,NAOXI
            CALL SCOPY(NSZBJ1,BUF1(IOFF1A),NSIZF1,BUFTMP,1)
c            write(*,*) 'before '
c            write(*,1234) (buftmp(i),i=1,nszbj1)
            CALL SCOPY(NSZBJ2,BUF1(IOFF1B),NSIZF1,
     &               BUFTMP(IOFFBA),1)
c            write(*,1234) (buftmp(IOFFBA-1+i),i=1,nszbj2)
            CALL SYTRAB(BUFTMP,BUFTMP(IOFFBA),EVEC(IOFFEVA),
     &                EVEC(IOFFEVB2),EVEC(IOFFEOA),EVEC(IOFFEOB2),
     &                NAOXI,NAOI,NAXI,NAI,NIXI,NII,
     &                BUFTMP(IOFFBB),BUFTMP(IOFFCC),0)
            CALL SSCAL(NAOXI*NAOI,FOURTH,BUFTMP(IOFFCC),1)
c            write(*,*) 'aos',indxnu,indxmu
c            write(*,*)
c            write(*,1234) (buftmp(IOFFCC+I-1),I=1,NAOXI*NAOI)
c1234        format(8(F10.7,1x))
            CALL PUTLST(BUFTMP(IOFFCC),IOFF1AA,1,1,3,ITHRU)
            IOFF1A=IOFF1A+1
            IOFF1AA=IOFF1AA+1
            IOFF1B=IOFF1B+1
180        CONTINUE
170       CONTINUE
C
         ELSE
C
          IOFF1A=1
          IOFF1AA=1+(ISTART-1)*NAOXI
          IOFF1B=ILCIXI
          IOFFBA=NSZBJ1+1
          IOFFBB=NSZBJ2+IOFFBA
          IOFFCC=IOFFBB+NAOXI*NAOI
C
          DO 270 INDXNU=ISTART,IEND
           DO 280 INDXMU=1,NAOXI
            CALL SCOPY(NSZBJ1,BUF1(IOFF1A),NSIZF1,BUFTMP,1)
            CALL SCOPY(NSZBJ2,BUF1(IOFF1B),NSIZF1,
     &               BUFTMP(IOFFBA),1)
            CALL GETLST(BUFTMP(IOFFCC),IOFF1AA,1,1,3,ITHRU)
            CALL SSCAL(NAOXI*NAOI,FOUR,BUFTMP(IOFFCC),1)
            CALL SYTRAB(BUFTMP,BUFTMP(IOFFBA),EVEC(IOFFEVA),
     &                EVEC(IOFFEVB2),EVEC(IOFFEOA),EVEC(IOFFEOB2),
     &                NAOXI,NAOI,NAXI,NAI,NIXI,NII,
     &                BUFTMP(IOFFBB),BUFTMP(IOFFCC),1)
            CALL SSCAL(NAOXI*NAOI,FOURTH,BUFTMP(IOFFCC),1)
            CALL PUTLST(BUFTMP(IOFFCC),IOFF1AA,1,1,3,ITHRU)
            IOFF1A=IOFF1A+1
            IOFF1AA=IOFF1AA+1
            IOFF1B=IOFF1B+1
280        CONTINUE
270       CONTINUE
C
         ENDIF
C
         write(*,*) 'istart, iend ',istart,iend 

         ISTART=ISTART+NAOI1
C
        IF(IEND.NE.NAOI) GO TO 10
C
120     CONTINUE
100    CONTINUE
C
1000  CONTINUE
      RETURN
      END
