      SUBROUTINE T1RC01_R(ICORE,MAXCOR,IUHF,LAMBDA,IOFFLIST,LISTT,
     &           LISTTOFF)
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
      LOGICAL RHF
      DIMENSION ICORE(MAXCOR),IOFFT1(8,2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     &             NF1BB,NF2BB
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      DATA ONE /1.0/
      DATA ZILCH /0.0/
      DATA ONEM /-1.0/
      integer  aces_list_rows, aces_list_cols
      external aces_list_rows, aces_list_cols
      CHANGE=.TRUE.
      RHF=.FALSE.
      IF(IUHF.EQ.0)RHF=.TRUE.
C
C FIRST PICK UP T1 VECTOR.
C
      CALL GETR1(ICORE,MAXCOR,MXCOR,IUHF,IOFFT1,LISTT,LISTTOFF)
C
C SPIN CASES ABAB AND BABA, RESPECTIVELY.
C
      DO 10 ISPIN=1,1+IUHF
       LSTTAR=55+ISPIN+IOFFLIST
       IF(LAMBDA)THEN
        LSTOUT=119-ISPIN
        FACTOR=ONEM
       ELSE
        LSTOUT=LSTTAR
        FACTOR=ONE
       ENDIF
C
C RESET THE TARGET LIST PARAMETERS BECAUSE IT MUST FIRST BE USED AS A
C  SCRATCH LIST FOR A LIST WHICH IS PACKED DIFFERENTLY.
C
       IF(IUHF.NE.0)THEN
        LSTSCR=40-ISPIN
        CALL NEWTYP(LSTOUT,ISYTYP(1,LSTSCR),ISYTYP(2,LSTSCR),CHANGE)
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
C COMPUTE DIMENSIONS OF TARGET MATRIX.
C
        DSZTAR=IRPDPD(IRREPDO,ISYTYP(1,LSTOUT))
        DISTAR=IRPDPD(IRREPDO,ISYTYP(2,LSTOUT))
        DSZTMP=IRPDPD(IRREPDO,13-ISPIN)
        DISTMP=IRPDPD(IRREPDO,10+ISPIN)
C
C FIRST DO       W(MbEj) =   SUM T(j,f) * <Ef|Mb> (ISPIN=1)
C                W(mBeJ) =   SUM T(J,F) * <Fe|Bm> (ISPIN=2 OR RHF).
C
C THIS PRODUCT IS INITIALLY PACKED Mb,Ej [Bm,eJ].
C
       CALL IZERO(ICORE,IINTFP*DSZTMP*DISTMP)
        ALPHA=ONE*FACTOR
        BETA=ONE
        LISTW1=28+ISPIN
        IF(RHF)LISTW1=30
        DISW  =aces_list_cols(IRREPDO,LISTW1)
        DSZW  =aces_list_rows(IRREPDO,LISTW1)
        I000  =1
        I010  =I000+IINTFP*DISTMP*DSZTMP
        I020  =I010+IINTFP*DSZW*DISW
        I030  =I020+IINTFP*DSZW
C
C READ INTEGRALS INTO AN I(Mb,Ef) [I(Bm,Fe)] MATRIX
C
        CALL GETTRN(ICORE(I010),ICORE(I020),DSZW,DISW,2,IRREPDO,
     &              LISTW1)
        ITMP=DISW
        DISW=DSZW
        DSZW=ITMP
C
C TRANSPOSE KET INDICES IF ISPIN=2
C
        IF(ISPIN.EQ.2.OR.RHF)THEN
         I030=I020+IINTFP*DSZW
         I040=I030+IINTFP*DSZW
         I050=I040+DSZW
         CALL SYMTR1(IRREPDO,VRT(1,1),VRT(1,2),DSZW,ICORE(I010),
     &               ICORE(I020),ICORE(I030),ICORE(I040))
        ENDIF
C
C WE NOW HAVE
C
C                   I(Mb,Ef)    [ISPIN=1]
C                   I(Bm,eF)    [ISPIN=2 OR RHF]
C
C SETUP AND PERFORM THE MATRIX MULTIPLY
C
C                 I(MbE,f) * T(f,j) [ISPIN=1]
C                 I(Bme,F) * T(F,J) [ISPIN=2 OR RHF]
C
C
         IOFFI=I010
         IOFFZ=I000
         DO 100 IRREPF=1,NIRREP
          IRREPE=DIRPRD(IRREPF,IRREPDO)
          IRREPJ=IRREPF
          NROWI=DSZW*VRT(IRREPE,ISPIN)
          NCOLI=VRT(IRREPF,3-ISPIN)
          NROWT=VRT(IRREPF,3-ISPIN)
          NCOLT=POP(IRREPJ,3-ISPIN)
          NROWZ=NROWI
          NCOLZ=NCOLT
          IOFFT=IOFFT1(IRREPF,3-ISPIN)
          IF(MIN(NROWZ,NCOLZ,NCOLI).GT.0)THEN
           CALL XGEMM('N','N',NROWZ,NCOLZ,NCOLI,ALPHA,ICORE(IOFFI),
     &                NROWI,ICORE(IOFFT),NROWT,BETA,ICORE(IOFFZ),NROWZ)
          ENDIF
          IOFFI=IOFFI+IINTFP*NROWI*NCOLI
          IOFFZ=IOFFZ+IINTFP*NROWZ*NCOLZ
100      CONTINUE
C
C NOW WE HAVE A Mb-Ej (ISPIN=1) OR Bm-eJ (ISPIN=2 OR RHF) ORDERED
C  QUANTITY.  THE NEXT PIECE WILL BE ORDERED Ej-Mb (ISPIN=1) OR
C  Je-mB (ISPIN=2 OR RHF), SO WE NEED TO REORDER WHAT WE HAVE TO MATCH
C  THIS, THEREBY ALLOWING ACCUMULATION IN MATRIX MULTIPLY OPERATIONS.
C
        IF(ISPIN.EQ.2.OR.RHF)THEN
         I020=I010+IINTFP*MAX(DSZTMP,DISTMP)
         I030=I020+IINTFP*MAX(DSZTMP,DISTMP)
         I040=I030+IINTFP*MAX(DSZTMP,DISTMP)
         CALL SYMTR1(IRREPDO,VRT(1,2),POP(1,1),DSZTMP,ICORE(I000),
     &               ICORE(I010),ICORE(I020),ICORE(I030))
         CALL SYMTR3(IRREPDO,VRT(1,1),POP(1,2),DSZTMP,DISTMP,
     &               ICORE(I000),ICORE(I010),ICORE(I020),ICORE(I030))
        ENDIF
        I020=I010+IINTFP*DSZTMP*DISTMP
        CALL TRANSP(ICORE(I000),ICORE(I010),DISTMP,DSZTMP)
c YAU : old
c       CALL ICOPY(IINTFP*DISTMP*DSZTMP,ICORE(I010),1,ICORE(I000),1)
c YAU : new
        CALL DCOPY(DISTMP*DSZTMP,ICORE(I010),1,ICORE(I000),1)
c YAU : end
C
C NOW DO    W(MbEj) =   - T(n,b) * <Mn|Ej>   (ISPIN=1)
C           W(mBeJ) =   - T(N,B) * <Nm|Je>   (ISPIN=2 OR RHF)
C
        ALPHA=ONEM*FACTOR
        BETA=ONE
        LISTW2=8+ISPIN
        IF(RHF)LISTW2=10
        DSZTMP=IRPDPD(IRREPDO,10+ISPIN)
        DISTMP=IRPDPD(IRREPDO,13-ISPIN)
        DISW  =aces_list_cols(IRREPDO,LISTW2)
        DSZW  =aces_list_rows(IRREPDO,LISTW2)
        I020  =I010+IINTFP*DSZW*DISW
        I030  =I020+IINTFP*DSZW
        CALL GETTRN(ICORE(I010),ICORE(I020),DSZW,DISW,2,IRREPDO,LISTW2)
        ITMP=DISW
        DISW=DSZW
        DSZW=ITMP
C
C WE NOW HAVE
C                I(Ej,Mn) [ISPIN=1]
C                I(Je,Nm) [ISPIN=2 OR RHF]
C
C WE NEED TO HAVE N AS THE SLOWEST INDEX, SO TRANSPOSE KET INDICES FOR ISPIN=2.
C
        IF(ISPIN.EQ.2.OR.RHF)THEN
         I030=I020+IINTFP*DSZW
         I040=I030+IINTFP*DSZW
         I050=I040+IINTFP*DSZW
         CALL SYMTR1(IRREPDO,POP(1,1),POP(1,2),DSZW,ICORE(I010),
     &               ICORE(I020),ICORE(I030),ICORE(I040))
        ENDIF
C
C WE NOW HAVE
C                I(Ej,Mn) [ISPIN=1]
C                I(Je,mN) [ISPIN=2 OR RHF]
C
C
C AND WE CAN DO THE MATRIX MULTIPLY
C                                 t
C                I(Ej,Mn) * T(b,n)  [ISPIN=1]
C                                 t
C                I(Je,mN) * T(B,N)  [ISPIN=2 OR RHF]
C
        IOFFI=I010
        IOFFZ=I000
        DO 200 IRREPN=1,NIRREP
         IRREPM=DIRPRD(IRREPN,IRREPDO)
         IRREPB=IRREPN
         NROWI=DSZW*POP(IRREPM,ISPIN)
         NCOLI=POP(IRREPN,3-ISPIN)
         NROWT=VRT(IRREPB,3-ISPIN)
         NCOLT=POP(IRREPN,3-ISPIN)
         NROWZ=NROWI
         NCOLZ=NROWT
         IOFFT=IOFFT1(IRREPB,3-ISPIN)
         IF(MIN(NROWZ,NCOLZ,NCOLI).GT.0)THEN
          CALL XGEMM('N','T',NROWZ,NCOLZ,NCOLI,ALPHA,ICORE(IOFFI),
     &               NROWI,ICORE(IOFFT),NROWT,BETA,ICORE(IOFFZ),
     &               NROWZ)
         ENDIF
         IOFFI=IOFFI+IINTFP*NROWI*NCOLI
         IOFFZ=IOFFZ+IINTFP*NROWZ*NCOLZ
200     CONTINUE
C
C REORDER TO
C           Ej-Mb ->  Ej-bM [ISPIN=1]
C           Je-mB ->  eJ-Bm [ISPIN=2 OR RHF]
C
        I020=I010+IINTFP*DSZTMP
        I030=I020+IINTFP*DSZTMP
        I040=I030+IINTFP*DSZTMP
        CALL SYMTR1(IRREPDO,POP(1,ISPIN),VRT(1,3-ISPIN),DSZTMP,
     &              ICORE(I000),ICORE(I010),ICORE(I020),
     &              ICORE(I030))
        IF(ISPIN.EQ.2.OR.RHF)THEN
         I020=I010+IINTFP*MAX(DSZTMP,DISTMP)
         I030=I020+IINTFP*MAX(DSZTMP,DISTMP)
         I040=I030+IINTFP*MAX(DSZTMP,DISTMP)
         CALL SYMTR3(IRREPDO,POP(1,1),VRT(1,2),DSZTMP,DISTMP,
     &               ICORE(I000),ICORE(I010),ICORE(I020),
     &               ICORE(I030))
        ENDIF
        IF(ISPIN.EQ.1.AND..NOT.RHF)THEN
         SSTSPN='ABBA'
        ELSE
         SSTSPN='BAAB'
        ENDIF
C
C NOW WRITE THESE TO DISK FOR EACH IRREP.
C
        CALL PUTLST(ICORE(I000),1,DISTMP,1,IRREPDO,LSTOUT)
20     CONTINUE
C
C NOW SWITCH ORDERING
C
C           Ej-bM ->  EM-bj [ISPIN=1]
C           eJ-Bm ->  em-BJ [ISPIN=2 OR RHF]
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
     &             SSTSPN)
       IF(IUHF.NE.0)CALL NEWTYP(LSTOUT,8+ISPIN,11-ISPIN,CHANGE)
C
C  FOR LAMBDA UPDATE LSTTAR AND COPY ORIGINAL INTERMEDIATES TO LSTOUT
C
       IF(LAMBDA) THEN
        CALL GETALL(ICORE(I010),TARSIZ,1,LSTTAR)
        CALL PUTALL(ICORE(I010),TARSIZ,1,LSTOUT)
        CALL SAXPY(TARSIZ,ONE,ICORE(I010),1,ICORE(I000),1)
       ENDIF
C
       CALL PUTALL(ICORE(I000),TARSIZ,1,LSTTAR)
10    CONTINUE
      RETURN
      END
