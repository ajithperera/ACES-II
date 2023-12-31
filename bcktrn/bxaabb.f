      SUBROUTINE BXAABB(EVEC,BUF1,BUF2,BUF3,BUFTMP,LEN1,LEN2,LEN3,
     &                  LENTMP,IUHF)
C
C THIS ROUTINE TRANSFORMS T2 AMPLITUDES OF SYMMETRY TYPE
C  AABB FROM THE MO TO THE AO BASIS.  THE AMPLITUDES ARE
C  DEALT WITH IN MULLIKEN ORDER (AI,BJ) AND THE TRANSFORMATION
C  PASSES THROUGH TWO STAGES:
C
C    T(AI,BJ) -> T(XX,BJ) -> T(X<X,BJ) -> T(X<X,XX) -> T(X<X,X<X)
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
C       BUF1  = MAX(NTAIBJ,NTXXBJ)          (RHF AND UHF)
C       BUF2  = MAX(NTAIBJ,NXX)             (RHF AND UHF)
C       BUF3  = NTAIBJ                      (UHF ONLY)
C       BUFTMP= 2*NAO*NAO                   (RHF AND UHF)
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
      INTEGER POP,VRT,DIRPRD,DISSIZ,AOPOP 
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
C
C  LOOP OVER THE DIFFERENT TYPES
C
      DO 1000 ITYPE=1,IUHF+1 
C
       LISTW=118
C
       IF(ITYPE.EQ.1) THEN
C
C  AABB AND BBBB SPIN CASES
C 
        LISTW2=120
        ISL=1
        ISR=2
C
       ELSE
C
C  AAAA AND BBAA SPIN CASES
C
        LISTW2=119
        ISL=2
        ISR=1
C
       ENDIF
C
C READ IN FIRST IRREP OF SYMMETRY-PACKED T2 VECTOR.  THESE AMPLITUDES
C  BELONG TO ONE OF THE SYMMETRY TYPES AAAA OR AABB.  
C  EACH LOGICAL RECORD CONTAINS - ALL AA, FOLLOWED BY ALL BB, CC, ETC
C  FOR A SPECIFIC BJ PAIR WHICH IS ALSO AA,BB,CC ETC.  FROM THIS LIST,
C  THE AABB AMPLITUDES WILL BE EXTRACTED AND ORDERED
C
C                  1 1 2 2
C                  1 1 3 3
C                  2 2 3 3
C                  1 1 4 4 
C                  2 2 4 4 
C                  3 3 4 4 
C                  . . . .
C                h-1 h-1 h h
C
C WHERE EACH SYMMETRY BLOCK IS DESIGNATED BY A NUMBER.  THE AABB
C  AMPLITUDES WILL THEN BE TRANSFORMED ONE SET AT A TIME.
C
       DISSIZ=IRPDPD(1,ISYTYP(1,LISTW))
       NUMDIS=IRPDPD(1,ISYTYP(2,LISTW))
       IF(ITYPE.EQ.1) THEN
        CALL GETLST(BUF1,1,NUMDIS,1,1,LISTW)
       ELSE
C
C  FOR ITYPE EQ 2 WE WANT TO HAVE THE BBAA LIST ( TRANSPOSE OF AABB)
C
        CALL GETLST(BUF2,1,NUMDIS,1,1,LISTW)
        CALL TRANSP(BUF2,BUF1,NUMDIS,DISSIZ)
        ITMP=NUMDIS
        NUMDIS=DISSIZ
        DISSIZ=ITMP
       ENDIF
C
C NOW RESORT THE ELEMENTS IN THE LIST SO THAT THEY ARE ORDERED
C  ACCORDING TO THE SCHEME ABOVE.  (TRICKY STUFF)
C
       IOFF1S=1+DISSIZ*POP(1,ISR)*VRT(1,ISR)
       IOFF2=1
       DO 10 IRREPR=2,NIRREP
        IOFF0=0
        DO 20 IRREPL=1,IRREPR-1
         NSIZVL=VRT(IRREPL,ISL)
         NSIZVR=VRT(IRREPR,ISR)
         NSIZOL=POP(IRREPL,ISL)
         NSIZOR=POP(IRREPR,ISR)
         NSIZEL=NSIZVL*NSIZOL
         NSIZER=NSIZVR*NSIZOR
         IOFF1=IOFF1S+IOFF0
         DO 30 INDEXJ=1,NSIZOR
          DO 40 INDEXB=1,NSIZVR
           CALL SCOPY(NSIZEL,BUF1(IOFF1),1,BUF2(IOFF2),1)
           IOFF1=IOFF1+DISSIZ
           IOFF2=IOFF2+NSIZEL
40        CONTINUE
30       CONTINUE
         IOFF0=IOFF0+NSIZEL
20      CONTINUE
        IOFF1S=IOFF1S+NSIZER*DISSIZ
10     CONTINUE
C
C  FOR UHF TREAT AS WELL THE AAAA AND BBBB LIST, SCALE ALSO THE AABB
C  AND BBAA LIST WITH A FACTOR OF HALF ( IN RHF WE HAVE ONLY AABB 
C  WITH A FACTOR OF ONE, HERE WE HAVE AABB AND BBAA THUS A FACTOR OF HALF  
C
       IF(IUHF.EQ.1) THEN
C
        CALL SSCAL(IOFF2-1,HALF,BUF2,1)
C
        DISSIZ=IRPDPD(1,ISYTYP(1,LISTW2))
        NUMDIS=IRPDPD(1,ISYTYP(2,LISTW2))
        CALL GETLST(BUF1,1,NUMDIS,1,1,LISTW2)
C
C NOW RESORT THE ELEMENTS IN THE LIST SO THAT THEY ARE ORDERED
C  ACCORDING TO THE SCHEME ABOVE.  (TRICKY STUFF)
C
        IOFF1S=1+DISSIZ*POP(1,ISR)*VRT(1,ISR)
        IOFF2=1
        DO 15 IRREPR=2,NIRREP
         IOFF0=0
         DO 25 IRREPL=1,IRREPR-1
          NSIZVL=VRT(IRREPL,ISR)
          NSIZVR=VRT(IRREPR,ISR)
          NSIZOL=POP(IRREPL,ISR)
          NSIZOR=POP(IRREPR,ISR)
          NSIZEL=NSIZVL*NSIZOL
          NSIZER=NSIZVR*NSIZOR
          IOFF1=IOFF1S+IOFF0
          DO 35 INDEXJ=1,NSIZOR
           DO 45 INDEXB=1,NSIZVR
            CALL SCOPY(NSIZEL,BUF1(IOFF1),1,BUF3(IOFF2),1)
            IOFF1=IOFF1+DISSIZ
            IOFF2=IOFF2+NSIZEL
45         CONTINUE
35        CONTINUE
          IOFF0=IOFF0+NSIZEL
25       CONTINUE
         IOFF1S=IOFF1S+NSIZER*DISSIZ
15      CONTINUE
C
       ENDIF
C
C  NOW RESORT ELEMENTS
C DRIVE TRANSFORMATION BY LOOPING OVER IRREP PAIRS
C
       IOFF2=1
       IOFF3=1
       ITHRU=0
       DO 110 IRREPR=2,NIRREP
        DO 120 IRREPL=1,IRREPR-1
         ITHRU=ITHRU+1
         CALL ZERO(BUF1,LEN1)
C
C FIRST HALF TRANSFORM THE AMPLITUDES TO T2(XX,BJ).  HALF TRANSFORMED
C  AMPLITUDES ARE HELD ON BUF1.
C
         NSIZAO=AOPOP(IRREPL)
         NSIZF=NSIZAO*NSIZAO
         NSIZT=(NSIZAO*(NSIZAO+1))/2
         NSIZVL=VRT(IRREPL,ISL)
         NSIZVR=VRT(IRREPR,ISR)
         NSIZOL=POP(IRREPL,ISL)
         NSIZOR=POP(IRREPR,ISR)
         NSIZEL=NSIZVL*NSIZOL
         NSIZER=NSIZVR*NSIZOR
         IOFFB=NSIZAO*NSIZAO+1
         IOFF1=1
         IOFFV=INDVRT(IRREPL,ISL)
         IOFFO=INDOCC(IRREPL,ISL)
         DO 130 INDEXJ=1,NSIZOR
          DO 140 INDEXB=1,NSIZVR
           CALL SYTRAA(BUF2(IOFF2),EVEC(IOFFV),EVEC(IOFFO),NSIZAO,
     &                 NSIZVL,NSIZOL,BUFTMP(IOFFB),BUFTMP,0)
           CALL MPMT(BUFTMP,NSIZAO)
C
C NOW COMPRESS MATRIX TO T2(X.LE.X;BJ) AND PUT DIRECTLY ON TARGET LIST
C
           CALL SQUEZ2(BUFTMP,BUF1(IOFF1),NSIZAO)
           IOFF1=IOFF1+NSIZT
           IOFF2=IOFF2+NSIZEL
140       CONTINUE
130      CONTINUE
C
C  FOR UHF WE HAVE TO SOME EXTRA WORK
C
         IF(IUHF.EQ.1) THEN
C
C
C FIRST HALF TRANSFORM THE AMPLITUDES TO T2(XX,BJ).  HALF TRANSFORMED
C  AMPLITUDES ARE HELD ON BUF1.
C
          NSIZVL=VRT(IRREPL,ISR)
          NSIZOL=POP(IRREPL,ISR)
          NSIZEL=NSIZVL*NSIZOL
          IOFF1=1
          IOFFV=INDVRT(IRREPL,ISR)
          IOFFO=INDOCC(IRREPL,ISR)
          DO 230 INDEXJ=1,NSIZOR
           DO 240 INDEXB=1,NSIZVR
            CALL SYTRAA(BUF3(IOFF3),EVEC(IOFFV),EVEC(IOFFO),NSIZAO,
     &                  NSIZVL,NSIZOL,BUFTMP(IOFFB),BUFTMP,0)
            CALL MPMT(BUFTMP,NSIZAO)
C
C NOW COMPRESS MATRIX TO T2(X.LE.X;BJ) AND PUT DIRECTLY ON TARGET LIST
C
            CALL SQUEZ3(BUFTMP,BUF1(IOFF1),NSIZAO)
            IOFF1=IOFF1+NSIZT
            IOFF3=IOFF3+NSIZEL
240        CONTINUE
230       CONTINUE
C
         ENDIF
C
C NOW DO SECOND HALF TRANSFORMATION 
C
C FIRST PART IS
C
C        T2(XX;XJ) = C(X,B) * T2(XX,BJ) 
C
C HALF-TRANSFORMED INTEGRALS MUST BE PLUCKED OFF OF BUF1
C
         IF(ITYPE.EQ.1) THEN
          NTOP=NSIZAO
          NSIZAO=AOPOP(IRREPR)
          IOFFB=NSIZAO*NSIZAO+1
          IOFFC=NSIZAO*NSIZAO+IOFFB
          IOFF1=1
          IOFFV=INDVRT(IRREPR,ISR)
          IOFFO=INDOCC(IRREPR,ISR)
          DO 150 INDXMU=1,NTOP
           DO 160 INDXNU=1,INDXMU
            CALL SCOPY(NSIZER,BUF1(IOFF1),NSIZT,BUFTMP(IOFFC),1)
            CALL SYTRAA(BUFTMP(IOFFC),EVEC(IOFFV),EVEC(IOFFO),NSIZAO,
     &                  NSIZVR,NSIZOR,BUFTMP(IOFFB),BUFTMP,0)
            CALL MPMT(BUFTMP,NSIZAO) 
c
c to compare with grnfnc
c
            CALL SSCAL(NSIZAO*NSIZAO,FOURTH,BUFTMP,1)
            CALL SQUEZ2(BUFTMP,BUFTMP(IOFFB),NSIZAO)
            CALL PUTLST(BUFTMP(IOFFB),IOFF1,1,1,2,ITHRU)
            IOFF1=IOFF1+1
160        CONTINUE
150       CONTINUE
C
         ELSE
          NTOP=NSIZAO
          NSIZAO=AOPOP(IRREPR)
          IOFFB=NSIZAO*NSIZAO+1
          IOFFC=NSIZAO*NSIZAO+IOFFB
          IOFF1=1
          IOFFV=INDVRT(IRREPR,ISR)
          IOFFO=INDOCC(IRREPR,ISR)
          DO 250 INDXMU=1,NTOP
           DO 260 INDXNU=1,INDXMU
            CALL SCOPY(NSIZER,BUF1(IOFF1),NSIZT,BUFTMP(IOFFC),1)
            CALL SYTRAA(BUFTMP(IOFFC),EVEC(IOFFV),EVEC(IOFFO),NSIZAO,
     &                  NSIZVR,NSIZOR,BUFTMP(IOFFB),BUFTMP,0)
            CALL MPMT(BUFTMP,NSIZAO) 
c
c to compare with grnfnc
c
            CALL SSCAL(NSIZAO*NSIZAO,FOURTH,BUFTMP,1)
            CALL GETLST(BUFTMP(IOFFB),IOFF1,1,1,2,ITHRU)
            CALL SQUEZ3(BUFTMP,BUFTMP(IOFFB),NSIZAO)
            CALL PUTLST(BUFTMP(IOFFB),IOFF1,1,1,2,ITHRU)
            IOFF1=IOFF1+1
260        CONTINUE
250       CONTINUE
C
         ENDIF
C
120     CONTINUE
110    CONTINUE
C
1000  CONTINUE
      RETURN
      END
