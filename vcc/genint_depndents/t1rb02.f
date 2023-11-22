      SUBROUTINE T1RB02(ICORE,MAXCOR,IUHF,LAMBDA)
C
C THIS SUBROUTINE COMPUTES TWO T1*W CONTRIBUTIONS TO THE
C  W(mbej) INTERMEDIATE.  
C
C     W(mBEj) =  SUM T(j,f) * <fE|mB> - SUM T(N,B) * <mN|jE>
C     W(MbeJ) =  SUM T(J,F) * <Fe|Mb> - SUM T(n,b) * <Mn|Je> (UHF only)
C
C ALSO COMPUTE ONE OF THE CONTRIBUTIONS TO THE F(ea) INTERMEDIATE:
C
C     F(EA)   = SUM T(m,f) * <fE|mA> 
C     F(EA)   = SUM T(M,F) * <Fe|Ma>  (UHF only)
C
C BY TAKING A GENERALIZED TRACE OVER THE FIRST TERM IN THE FIRST TWO EQUATIONS
C  ABOVE.
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL LAMBDA
      DOUBLE PRECISION ONE,ONEM,ZILCH,ALPHA,BETA,FACTOR
      CHARACTER*4 SSTSPN(2)
      LOGICAL INCORE,RHF
      DIMENSION ICORE(MAXCOR),IOFFT1(8,2),IOFFZL(8,4),SIZEVV(2)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     &             NF1BB,NF2BB
      COMMON /INFO/ NOCCO(2),NVRTO(2)
C
      DATA ONE /1.0/
      DATA ZILCH /0.0/
      DATA ONEM /-1.0/
      DATA SSTSPN /'ABAB','BABA'/
C
C FIRST PICK UP T1 VECTOR.
C
      CALL GETT1(ICORE,MAXCOR,MXCOR,IUHF,IOFFT1)
      RHF=.FALSE.
      IF(IUHF.EQ.0)RHF=.TRUE.
      LSTFEA=92
C
C COMPUTE SIZES OF A FULL VV DISTRIBUTION WHICH TRANSFORMS AS THE
C  TOTALLY SYMMETRIC REP.
C
      CALL IZERO(SIZEVV,2)
      DO 3000 ISPIN=1,2
         DO 3001 IRREP=1,NIRREP
            SIZEVV(ISPIN)=SIZEVV(ISPIN)+VRT(IRREP,ISPIN)*
     &                    VRT(IRREP,ISPIN)
 3001    CONTINUE
 3000 CONTINUE

C
C SPIN CASES BAAB AND ABBA, RESPECTIVELY.
C
      DO 10 ISPIN=1,1+IUHF
         LSTTAR=57+ISPIN
         IF(LAMBDA)THEN
            LSTOUT=124+ISPIN
            FACTOR=ONEM
         ELSE
            LSTOUT=LSTTAR
            FACTOR=ONE
         ENDIF
C
C LOOP OVER IRREPS - IN THE FIRST BLOCK OF CODE THIS CORRESPONDS TO mB,
C                    WHILE IT IS jE IN THE SECOND BLOCK.
C
         DO 20 IRREPDO=1,NIRREP

C
C COMPUTE OFFSETS INTO AN i,A DISTRIBUTION FOR THIS IRREP (RIGHT INDEX - 1)
C                         I,a DISTRIBUTION FOR THIS IRREP (RIGHT INDEX - 2)
C                         b,M DISTRIBUTION FOR THIS IRREP (RIGHT INDEX - 3)
C                         B,m DISTRIBUTION FOR THIS IRREP (RIGHT INDEX - 4)
C
            IOFF1=0
            IOFF2=0
            IOFF3=0
            IOFF4=0
            DO 1001 IRREPA=1,NIRREP
               IRREPM=IRREPA
               IRREPI=DIRPRD(IRREPA,IRREPDO)
               IRREPB=IRREPI
               IOFFZL(IRREPA,1)=IOFF1
               IOFFZL(IRREPA,2)=IOFF2
               IOFFZL(IRREPM,3)=IOFF3
               IOFFZL(IRREPM,4)=IOFF4
               IOFF1=IOFF1+POP(IRREPI,2)*VRT(IRREPA,1)
               IOFF2=IOFF2+POP(IRREPI,1)  *VRT(IRREPA,2)
               IOFF3=IOFF3+VRT(IRREPB,2)*POP(IRREPM,1)
               IOFF4=IOFF4+VRT(IRREPB,1)  *POP(IRREPM,2)
 1001       CONTINUE
C
C COMPUTE DIMENSIONS OF TARGET MATRIX.
C
            DISTAR=IRPDPD(IRREPDO,ISYTYP(2,LSTTAR))
            DSZTAR=IRPDPD(IRREPDO,ISYTYP(1,LSTTAR))
            TARSIZ=DISTAR*DSZTAR
C
C FIRST DO  W(mBEj) =   SUM T(j,f) * <Ef|Bm>  
C           W(MbeJ) =   SUM T(J,F) * <Fe|Mb>  (UHF only)
C
C THIS PRODUCT IS INITIALLY PACKED jE-Bm [Je-Mb].
C
            LISTW1=31-ISPIN
            DISW  =IRPDPD(IRREPDO,ISYTYP(2,LISTW1))
            DSZW  =IRPDPD(IRREPDO,ISYTYP(1,LISTW1))
C
C I000 HOLDS THE W(jEmB) TARGET FOR GAMMA(mB).
C I010 HOLDS AREA EVENTUALLY USED AS SCRATCH IN SYMTR1a
C I011 AND I012 SCRATCH ARRAYS FOR SYMTR1.
C I020 HOLDS THE <Ef||Bm> INTEGRALS.
C
            I000  =1
            I010  =I000+IINTFP*DISTAR*DSZTAR
            I011  =I010+IINTFP*MAX(DSZTAR,DISTAR,DSZW,SIZEVV(ISPIN))
            I012  =I011+IINTFP*MAX(DSZTAR,DISTAR)  
            I020  =I012+IINTFP*MAX(DSZTAR,DISTAR)
            I030  =I020+IINTFP*DSZW*DISW
            CALL IZERO(ICORE,IINTFP*DISTAR*DSZTAR)
C
            IF(I030.LE.MXCOR)THEN
               INCORE=.TRUE.
               CALL GETLST(ICORE(I020),1,DISW,2,IRREPDO,LISTW1)
            ELSE
               INCORE=.FALSE.
               I030=I020+IINTFP*DSZW
            ENDIF
C
            DO 30 INUMBM=1,DISW
               IF(INCORE)THEN
                  IOFFW1R=(INUMBM-1)*DSZW*IINTFP+I020
               ELSE
                  CALL GETLST(ICORE(I020),INUMBM,1,2,IRREPDO,LISTW1)
                  IOFFW1R=I020
               ENDIF
C
               IOFFW1L=0
               IOFFZR=(INUMBM-1)*DSZTAR*IINTFP+I000
               IF(ISPIN.EQ.1)THEN
C
                  DO 40 IRREPF=1,NIRREP
                     IRREPE=DIRPRD(IRREPF,IRREPDO)
                     IRREPT=IRREPF
                     IOFFT =IOFFT1(IRREPT,2)
                     IOFFW1=IOFFW1R+IOFFW1L
                     IOFFZ =IOFFZR+IOFFZL(IRREPE,1)*IINTFP
                     NROW=POP(IRREPT,2)
                     NCOL=VRT(IRREPE,1)
                     NSUM=VRT(IRREPT,2)
                     ALPHA=ONE*FACTOR
                     BETA=ZILCH
                     IF(MIN(NROW,NCOL,NSUM).GT.0)THEN
                        CALL XGEMM('T','T',NROW,NCOL,NSUM,ALPHA,
     &                              ICORE(IOFFT),NSUM,ICORE(IOFFW1),
     &                              NCOL,BETA,ICORE(IOFFZ),NROW)
                     ENDIF
                     IOFFW1L=IOFFW1L+NCOL*NSUM*IINTFP
 40               CONTINUE
               ELSE
C
                  DO 41 IRREPE=1,NIRREP
                     IRREPF=DIRPRD(IRREPE,IRREPDO)
                     IRREPT=IRREPF
                     IOFFT =IOFFT1(IRREPT,1)
                     IOFFW1=IOFFW1R+IOFFW1L
                     IOFFZ =IOFFZR+IOFFZL(IRREPE,2)*IINTFP
                     NROW=POP(IRREPT,1)
                     NCOL=VRT(IRREPE,2)
                     NSUM=VRT(IRREPT,1)
                     ALPHA=ONE*FACTOR
                     BETA=ZILCH
                     IF(MIN(NROW,NCOL,NSUM).GT.0)THEN
                        CALL XGEMM('T','N',NROW,NCOL,NSUM,ALPHA,
     &                              ICORE(IOFFT),NSUM,ICORE(IOFFW1),
     &                              NSUM,BETA,ICORE(IOFFZ),NROW)
                     ENDIF
                     IOFFW1L=IOFFW1L+NCOL*NSUM*IINTFP
 41               CONTINUE
               ENDIF
 30         CONTINUE
C
C COMPUTE CONTRIBUTION TO F(EA) AND F(ea) INTERMEDIATES.
C
            IF(.NOT.LAMBDA) THEN
               IF(ISPIN.EQ.1)THEN
                  CALL TRACEOO('OVVO',IRREPDO,POP(1,2),VRT(1,1),DISTAR,
     &                          SIZEVV(1),ICORE(I000),ICORE(I010))
                  CALL SUMSYM3(ICORE(I010),ICORE(I020),SIZEVV(1),
     &                         1,1,LSTFEA)
               ELSE
                  CALL TRACEOO('OVOV',IRREPDO,POP(1,1),VRT(1,2),DISTAR,
     &                          SIZEVV(2),ICORE(I000),ICORE(I010))
                  CALL SUMSYM3(ICORE(I010),ICORE(I020),SIZEVV(2),
     &                         1,2,LSTFEA)
               ENDIF
            ENDIF
C
C THE CODE ABOVE PRODUCES mbej INTERMEDIATES IN THE ORDER:
C
C                   jE-Bm (ISPIN=1)
C                   Je-Mb (ISPIN=2)
C
C THE NEXT PIECE WILL BE EVALUATED AS
C
C                   Bm-Ej (ISPIN=1, UHF)
C                   bM-Je (ISPIN=1, RHF)
C                   bM-Je (ISPIN=2)
C
C WE NEED TO REORDER WHAT WE HAVE SO THAT IT WILL MATCH UP AND
C  CAN BE ACCUMULATED BY XGEMM. 
C
            IF(ISPIN.EQ.1.AND..NOT.RHF)THEN
               CALL SYMTR3(IRREPDO,POP(1,2),VRT(1,1),DSZTAR,
     &                     DISTAR,ICORE(I000),ICORE(I010),  
     &                    ICORE(I011),ICORE(I012))
               CALL MTRAN2(ICORE(I000),DSZTAR)
            ELSEIF(ISPIN.EQ.2)THEN
               CALL SYMTR1(IRREPDO,POP(1,1),VRT(1,2),DSZTAR,
     &                     ICORE(I000),ICORE(I010),ICORE(I011),
     &                     ICORE(I012))
               CALL MTRAN2(ICORE(I000),DSZTAR)
            ELSEIF(RHF)THEN
               CALL MTRAN2(ICORE(I000),DSZTAR)
            ENDIF
C  
C
C NOW DO    W(mBEj) = W(mBEj) - SUM T(N,B) * <Nm|Ej> (ISPIN=1)
C           W(MbeJ) = W(MbeJ) - SUM T(n,b) * <Mn|Je> (ISPIN=2 OR RHF)
C
C THIS PRODUCT IS INITIALLY PACKED Bm-Ej [bM-Je].
C 
            LISTW2=8+ISPIN
            IF(RHF)LISTW2=10
            DISW  =IRPDPD(IRREPDO,ISYTYP(2,LISTW2))
            DSZW  =IRPDPD(IRREPDO,ISYTYP(1,LISTW2))
C
C I020 NOW HOLDS THE <Nm|Ej> [<Mn|Je>] INTEGRALS.
C
            I020  =I010+IINTFP*DSZW
            I030  =I020+IINTFP*DSZW*DISW
            IF(I030.LE.MXCOR)THEN
               INCORE=.TRUE.
               CALL GETLST(ICORE(I020),1,DISW,2,IRREPDO,LISTW2)
            ELSE
               INCORE=.FALSE.
               I030=I020+IINTFP*DSZW
            ENDIF
C
            DO 50 INUMEJ=1,DISW
               IF(INCORE)THEN
                  IOFFW2R=(INUMEJ-1)*DSZW*IINTFP+I020
               ELSE
                  CALL GETLST(ICORE(I020),INUMEJ,1,2,IRREPDO,LISTW2)
                  IOFFW2R=I020
               ENDIF
C
               IOFFW2L=0
               IOFFZR=(INUMEJ-1)*DSZTAR*IINTFP+I000
               IF(ISPIN.EQ.1.AND..NOT.RHF)THEN
C
                  DO 60 IRREPM=1,NIRREP
                     IRREPN=DIRPRD(IRREPM,IRREPDO)
                     IRREPT=IRREPN
                     IOFFT=IOFFT1(IRREPT,1)
                     IOFFW2=IOFFW2R+IOFFW2L
                     IOFFZ=IOFFZR+IOFFZL(IRREPM,4)*IINTFP
                     NROW=VRT(IRREPT,1)
                     NCOL=POP(IRREPM,2)
                     NSUM=POP(IRREPT,1)
                     ALPHA=ONEM*FACTOR
                     BETA=ONE
                     IF(MIN(NROW,NCOL,NSUM).GT.0)THEN
                        CALL XGEMM('N','N',NROW,NCOL,NSUM,ALPHA,
     &                              ICORE(IOFFT),NROW,ICORE(IOFFW2),
     &                              NSUM,BETA,ICORE(IOFFZ),NROW)
                     ENDIF
                     IOFFW2L=IOFFW2L+NCOL*NSUM*IINTFP
 60               CONTINUE
               ELSEIF(ISPIN.EQ.2.OR.RHF)THEN
C
                  DO 61 IRREPN=1,NIRREP
                     IRREPM=DIRPRD(IRREPN,IRREPDO)
                     IRREPT=IRREPN 
                     IOFFT=IOFFT1(IRREPT,2)
                     IOFFW2=IOFFW2R+IOFFW2L
                     IOFFZ=IOFFZR+IOFFZL(IRREPM,3)*IINTFP
                     NROW=VRT(IRREPT,2)
                     NCOL=POP(IRREPM,1)
                     NSUM=POP(IRREPT,2)
                     ALPHA=ONEM*FACTOR
                     BETA=ONE
                     IF(MIN(NROW,NCOL,NSUM).GT.0)THEN
                        CALL XGEMM('N','T',NROW,NCOL,NSUM,ALPHA,
     &                              ICORE(IOFFT),NROW,ICORE(IOFFW2),
     &                              NCOL,BETA,ICORE(IOFFZ),NROW)
                     ENDIF
                     IOFFW2L=IOFFW2L+NCOL*NSUM*IINTFP
 61               CONTINUE
               ENDIF
 50         CONTINUE
C
            IF(ISPIN.EQ.2.OR.RHF)THEN
               CALL SYMTR1(IRREPDO,POP(1,1),VRT(1,2),DSZTAR,
     &                     ICORE(I000),ICORE(I010),ICORE(I011),
     &                     ICORE(I012))
            ENDIF
C
C NOW WRITE THESE TO DISK FOR EACH IRREP.  THESE ARE ORDERED 
C
C                   Bm-Ej (ISPIN=1, UHF)
C                   bM-eJ (ISPIN=2 OR RHF)
C
            CALL PUTLST(ICORE(I000),1,DISTAR,1,IRREPDO,LSTOUT)
 20      CONTINUE
C
C NOW SWITCH ORDERING FROM Bm-Ej TO Bj-Em
C
         ISCSIZ=(NVRTO(1)+NVRTO(2))*(NOCCO(1)+NOCCO(2))
         TARSIZ=ISYMSZ(ISYTYP(1,LSTTAR),ISYTYP(2,LSTTAR))
         I000=1
         I010=I000+TARSIZ*IINTFP
         I020=I010+TARSIZ*IINTFP
         I030=I020+ISCSIZ
         IF(I030.GT.MXCOR)CALL INSMEM('T1RABBA',I030,MXCOR)
         CALL GETALL(ICORE(I010),TARSIZ,1,LSTOUT)
         CALL SSTRNG(ICORE(I010),ICORE(I000),TARSIZ,TARSIZ,ICORE(I020),
     &               SSTSPN(ISPIN))
C
C NOW DO TRANSPOSITION TO GET E,m-B,j ORDERING IRREP BY IRREP.
C
         IOFF=1
         DO 5000 IRREP=1,NIRREP
            DSZTAR=IRPDPD(IRREP,ISYTYP(1,LSTTAR))
            DISTAR=DSZTAR
            CALL MTRAN2(ICORE(IOFF),DSZTAR)
c CALL PUTLST(ICORE(IOFF),1,DISTAR,2,IRREP,LSTTAR)
            IOFF=IOFF+IINTFP*DISTAR*DSZTAR
 5000    CONTINUE
C
C NOW WRITE THE (PROPERLY ORDERED) PIECE TO THE TARGET LIST.
C
C FOR LAMBDA UPDATE THE TARGET LIST AND COPY ORIGINAL INTERMEDIATES
C  TO LSTOUT
C
         IF(LAMBDA) THEN
            CALL GETALL(ICORE(I010),TARSIZ,1,LSTTAR)
            CALL PUTALL(ICORE(I010),TARSIZ,1,LSTOUT)
            CALL SAXPY(TARSIZ,ONE,ICORE(I010),1,ICORE(I000),1)
         ENDIF
C
         CALL PUTALL(ICORE(I000),TARSIZ,1,LSTTAR)
 10   CONTINUE
C
      RETURN
      END
