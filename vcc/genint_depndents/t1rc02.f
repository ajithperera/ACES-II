      SUBROUTINE T1RC02(ICORE,MAXCOR,IUHF,LAMBDA)
C
C THIS SUBROUTINE COMPUTES TWO T1*W CONTRIBUTIONS TO THE
C  W(mbej) INTERMEDIATE.  
C
C     W(mBeJ) =   SUM T(J,F) * <Fe|Bm> - T(N,B) * <Nm|Je>
C     W(MbEj) =   SUM T(j,f) * <fE|bM> - T(n,b) * <nM|jE>
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL LAMBDA,CHANGE
      DOUBLE PRECISION ONE,ONEM,ZILCH,ALPHA,BETA,FACTOR
      CHARACTER*4 SSTSPN
      LOGICAL INCORE,RHF
      DIMENSION ICORE(MAXCOR),IOFFT1(8,2),IOFFZL(8,4)
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
      CHANGE=.TRUE.
      RHF=.FALSE.
      IF(IUHF.EQ.0)RHF=.TRUE.
C
C FIRST PICK UP T1 VECTOR.
C
      CALL GETT1(ICORE,MAXCOR,MXCOR,IUHF,IOFFT1)
C
C SPIN CASES ABAB AND BABA, RESPECTIVELY.
C
      DO 10 ISPIN=1,1+IUHF
         LSTTAR=55+ISPIN
         IF(LAMBDA)THEN
            LSTOUT=119-ISPIN
            FACTOR=ONEM
         ELSE
            LSTOUT=LSTTAR
            FACTOR=ONE
         ENDIF
 
C
C RESET THE TARGET LIST PARAMETERS BECAUSE IT MUST FIRST BE USED AS A
C SCRATCH LIST FOR A LIST WHICH IS PACKED DIFFERENTLY.
C
         IF(IUHF.NE.0)THEN
            LSTSCR=37+ISPIN
            CALL NEWTYP(LSTOUT,ISYTYP(1,LSTSCR),ISYTYP(2,LSTSCR),
     &                  CHANGE)
         ENDIF
C
C LOOP OVER IRREPS - IN THE FIRST BLOCK OF CODE THIS CORRESPONDS TO mB,
C                    WHILE IT IS jE IN THE SECOND BLOCK.
C
C     W(MbEj) =   SUM T(j,f) * <fE|bM> - T(n,b) * <nM|jE>
C     W(mBeJ) =   SUM T(J,F) * <Fe|Bm> - T(N,B) * <Nm|Je>
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
               IOFF2=IOFF2+POP(IRREPI,1)*VRT(IRREPA,2)
               IOFF3=IOFF3+VRT(IRREPB,2)*POP(IRREPM,1)
               IOFF4=IOFF4+VRT(IRREPB,1)*POP(IRREPM,2)
 1001       CONTINUE
C
C COMPUTE DIMENSIONS OF TARGET MATRIX.
C
            LSTTMP=40-ISPIN
            DSZTAR=IRPDPD(IRREPDO,ISYTYP(1,LSTOUT))
            DISTAR=IRPDPD(IRREPDO,ISYTYP(2,LSTOUT))
            DSZTMP=IRPDPD(IRREPDO,ISYTYP(1,LSTTMP))
            DISTMP=IRPDPD(IRREPDO,ISYTYP(2,LSTTMP))
C
C FIRST DO       W(MbEj) =   SUM T(j,f) * <Ef|Mb> (ISPIN=1)
C                W(mBeJ) =   SUM T(J,F) * <Fe|Bm> (ISPIN=2 OR RHF).
C
C THIS PRODUCT IS INITIALLY PACKED j,E-M,b [J,e;B,m].
C
            LISTW1=28+ISPIN
            IF(RHF)LISTW1=30
            DISW  =IRPDPD(IRREPDO,ISYTYP(2,LISTW1))
            DSZW  =IRPDPD(IRREPDO,ISYTYP(1,LISTW1))
C
C I000 HOLDS THE W(jEMb) [W(JemB)] TARGET FOR GAMMA(Mb [Bm]).
C I010 HOLDS AREA EVENTUALLY USED AS SCRATCH IN SYMTR.
C I020 HOLDS THE <Ef|Mb> [<Fe|Bm>] INTEGRALS.
C
            I000  =1
            I010  =I000+IINTFP*DISTMP*DSZTMP
            CALL IZERO(ICORE,IINTFP*DSZTMP*DISTMP)
            I011  =I010+IINTFP*MAX(DISTMP,DSZTMP)
            I012  =I011+IINTFP*MAX(DSZTMP,DISTMP)
            I020  =I012+IINTFP*MAX(DSZTMP,DISTMP)
            I030  =I020+IINTFP*DSZW*DISW
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
               IOFFZR=(INUMBM-1)*DSZTMP*IINTFP+I000
C
C FIRST DO       W(MbEj) =   SUM T(j,f) * <Ef|Mb> (ISPIN=1)
C                W(mBeJ) =   SUM T(J,F) * <Fe|Bm> (ISPIN=2 OR RHF).
C
               IF(ISPIN.EQ.1.AND..NOT.RHF)THEN
                  DO 40 IRREPF=1,NIRREP
                     IRREPT=IRREPF
                     IRREPE=DIRPRD(IRREPF,IRREPDO) 
                     IOFFT=IOFFT1(IRREPT,2)
                     IOFFW1=IOFFW1R+IOFFW1L
                     IOFFZ=IOFFZR+IOFFZL(IRREPE,1)*IINTFP
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
C
               ELSEIF(ISPIN.EQ.E2.OR.RHF)THEN
C
                  DO 41 IRREPE=1,NIRREP
                     IRREPT=DIRPRD(IRREPE,IRREPDO)
                     IRREPF=IRREPT
                     IOFFT=IOFFT1(IRREPT,1)
                     IOFFW1=IOFFW1R+IOFFW1L
                     IOFFZ=IOFFZR+IOFFZL(IRREPE,2)*IINTFP
                     NROW=POP(IRREPT,1)
                     NCOL=VRT(IRREPE,2)
                     NSUM=VRT(IRREPF,1)
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
C NOW WE HAVE A jE-Mb (ISPIN=1) OR Je-Bm (ISPIN=2 OR RHF) ORDERED
C QUANTITY.  THE NEXT PIECE WILL BE ORDERED bM-Ej (ISPIN=1) OR
C Bm-Je (ISPIN=2 OR RHF), SO WE NEED TO REORDER WHAT WE HAVE TO MATCH
C THIS, THEREBY ALLOWING ACCUMULATION IN MATRIX MULTIPLY OPERATIONS.
C
            IF(ISPIN.EQ.1.AND..NOT.RHF)THEN
               CALL SYMTR1(IRREPDO,POP(1,1),VRT(1,2),DSZTMP,ICORE(I000),
     &                     ICORE(I010),ICORE(I011),ICORE(I012))
               CALL SYMTR3(IRREPDO,POP(1,2),VRT(1,1),DSZTMP,DISTMP,
     &                     ICORE(I000),ICORE(I010),ICORE(I011),
     &                     ICORE(I012))
            ENDIF
            CALL TRANSP(ICORE(I000),ICORE(I010),DISTMP,DSZTMP)
c YAU : old
c       CALL ICOPY(DISTMP*DSZTMP*IINTFP,ICORE(I010),1,ICORE(I000),1)
c YAU : new
            CALL DCOPY(DISTMP*DSZTMP,ICORE(I010),1,ICORE(I000),1)
c YAU : end
C
C NOW DO    W(MbEj) =   - T(n,b) * <Mn|Ej>   (ISPIN=1)
C           W(mBeJ) =   - T(N,B) * <Nm|Je>   (ISPIN=2 OR RHF)
C
            LISTW2=8+ISPIN
            IF(RHF)LISTW2=10
            LSTTMP=37+ISPIN
            IF(RHF)LSTTMP=39 
            DSZTMP=IRPDPD(IRREPDO,ISYTYP(1,LSTTMP))
            DISTMP=IRPDPD(IRREPDO,ISYTYP(2,LSTTMP))
            DISW  =IRPDPD(IRREPDO,ISYTYP(2,LISTW2))
            DSZW  =IRPDPD(IRREPDO,ISYTYP(1,LISTW2))
C
C I020 NOW HOLDS THE <Mn|Ej> [<Nm|jE>] INTEGRALS.
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
               IOFFZR=(INUMEJ-1)*DSZTMP*IINTFP+I000
C
C NOW DO    W(MbEj) =   - T(n,b) * <Mn|Ej>   (ISPIN=1)
C           W(mBeJ) =   - T(N,B) * <Nm|Je>   (ISPIN=2 OR RHF)
C
               IF(ISPIN.EQ.1.AND..NOT.RHF)THEN
                  DO 60 IRREPN=1,NIRREP
                     IRREPM=DIRPRD(IRREPN,IRREPDO)
                     IRREPT=IRREPN
                     IOFFT=IOFFT1(IRREPT,2)
                     IOFFW2=IOFFW2R+IOFFW2L
                     IOFFZ=IOFFZR+IOFFZL(IRREPM,3)*IINTFP
                     NROW=VRT(IRREPT,2)
                     NCOL=POP(IRREPM,1)
                     NSUM=POP(IRREPN,2)
                     ALPHA=ONEM*FACTOR
                     BETA=ONE
                     IF(MIN(NROW,NCOL,NSUM).GT.0)THEN
                        CALL XGEMM('N','T',NROW,NCOL,NSUM,ALPHA,
     &                              ICORE(IOFFT),NROW,ICORE(IOFFW2),
     &                              NCOL,BETA,ICORE(IOFFZ),NROW)
                     ENDIF
                     IOFFW2L=IOFFW2L+NCOL*NSUM*IINTFP
 60               CONTINUE
C
               ELSEIF(ISPIN.EQ.2.OR.RHF)THEN
                  DO 61 IRREPM=1,NIRREP
                     IRREPN=DIRPRD(IRREPM,IRREPDO)
                     IRREPT=IRREPN
                     IOFFT=IOFFT1(IRREPT,1)
                     IOFFW2=IOFFW2R+IOFFW2L
                     IOFFZ=IOFFZR+IOFFZL(IRREPM,4)*IINTFP
                     NROW=VRT(IRREPT,1)
                     NCOL=POP(IRREPM,2)
                     NSUM=POP(IRREPN,1)
                     ALPHA=ONEM*FACTOR
                     BETA=ONE
                     IF(MIN(NROW,NCOL,NSUM).GT.0)THEN
                        CALL XGEMM('N','N',NROW,NCOL,NSUM,ALPHA,
     &                              ICORE(IOFFT),NROW,ICORE(IOFFW2),
     &                              NSUM,BETA,ICORE(IOFFZ),NROW)
                     ENDIF
                     IOFFW2L=IOFFW2L+NCOL*NSUM*IINTFP
 61               CONTINUE
               ENDIF
 50         CONTINUE
C
C REORDER TO 
C           bM-Ej ->  bM-Ej (ISPIN=1)
C           Bm-Je ->  Bm-eJ (ISPIN=2 OR RHF)
C     
            IF(ISPIN.EQ.1.AND..NOT.RHF)THEN
               CONTINUE
               SSTSPN='BAAB'
            ELSE
               CALL SYMTR1(IRREPDO,POP(1,1),VRT(1,2),DSZTMP,
     &                     ICORE(I000),ICORE(I010),ICORE(I011),
     &                     ICORE(I012))
               SSTSPN='ABBA'
            ENDIF
C
C NOW WRITE THESE TO DISK FOR EACH IRREP.  
C
            CALL PUTLST(ICORE(I000),1,DISTMP,1,IRREPDO,LSTOUT)
 20      CONTINUE
C
C NOW SWITCH ORDERING 
C
C                   bM-Ej -> bj-EM (ISPIN=1)
C                   Bm-eJ -> BJ-em (ISPIN=2 OR RHF)
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
     &               SSTSPN)
C
C NOW TRANSPOSE AND WRITE TO DISK IRREP BY IRREP.
C
         IOFF=1
         IF(IUHF.NE.0)CALL NEWTYP(LSTOUT,8+ISPIN,11-ISPIN,CHANGE)
         DO 5000 IRREP=1,NIRREP
            DSZTAR=IRPDPD(IRREP,ISYTYP(1,LSTTAR))
            DISTAR=IRPDPD(IRREP,ISYTYP(2,LSTTAR))
            CALL TRANSP(ICORE(IOFF),ICORE(I010),DSZTAR,DISTAR)
C
C FOR LAMBDA UPDATE TARGET LIST AND COPY ORIGINAL LIST TO LSTOUT
C
            IF(LAMBDA) THEN
               CALL GETLST(ICORE(IOFF),1,DISTAR,2,IRREP,LSTTAR)
               CALL PUTLST(ICORE(IOFF),1,DISTAR,2,IRREP,LSTOUT)
               CALL SAXPY(DSZTAR*DISTAR,ONE,ICORE(IOFF),1,ICORE(I010),1)
            ENDIF
C
            CALL PUTLST(ICORE(I010),1,DISTAR,2,IRREP,LSTTAR)
            IOFF=IOFF+IINTFP*DISTAR*DSZTAR
 5000    CONTINUE 
 10   CONTINUE
C
      RETURN
      END
