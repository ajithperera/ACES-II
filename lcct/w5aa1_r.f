      SUBROUTINE W5AA1_R(ICORE,MAXCOR,IUHF,TERM1,TERM2,TERM3,TAU,
     &                   IOFFLIST,LISTT,LISTTOFF)
C
C           Z(amef) = T1(g,m) * W (efag)
C
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL CCGRAD,TERM1,TERM2,TERM3,TAU
      DOUBLE PRECISION ONE,ONEM,ZILCH,ALPHA,BETA,HALF
      DIMENSION ICORE(MAXCOR),MNFULL(8),ABFULL(8),IOFFT1(8,2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /FLAGS/ IFLAGS(100)
C
      DATA ONE   /1.0/
      DATA ONEM  /-1.0/
      DATA ZILCH /0.0/
C

      CALL GETR1(ICORE,MAXCOR,MXCOR,IUHF,IOFFT1,LISTT,LISTTOFF)
      BETA=ONE
      DO 10 ISPIN=1,2
       LSTTAR=26+ISPIN+IOFFLIST
       DO 5 IRREP=1,NIRREP
        ABFULL(IRREP)=IRPDPD(IRREP,18+ISPIN)
        MNFULL(IRREP)=IRPDPD(IRREP,20+ISPIN)
5      CONTINUE
       DO 20 IRREPDO=1,NIRREP
        TARDSZ=IRPDPD(IRREPDO,ISYTYP(1,26+ISPIN))
        TARDIS=IRPDPD(IRREPDO,ISYTYP(2,26+ISPIN))
        TARSIZ=TARDSZ*TARDIS
        I000=1
        I010=I000+IINTFP*TARSIZ
        CALL IZERO(ICORE,IINTFP*TARSIZ)
        IF(TERM1)THEN
        ALPHA=ONE
        IF(IFLAGS(93).EQ.2)RETURN
C
C FIRST CONTRACTION:
C
C           Z(amef) = T1(g,m) * W (efag)
C
C           Z(AMEF) =  T1(G,M * W (EFAG)    [ISPIN = 1]
C
C           Z(amef) = T1(g,m) * W (efag)    [ISPIN = 2]
C
        LISTW =230+ISPIN
        DISFUL=ABFULL(IRREPDO)
        WDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTW))
        WDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTW))
        I020=I010+IINTFP*DISFUL*WDSZ
        IF(I020.LE.MXCOR)THEN
C
C DO IN-CORE ALGORITHM
C
         CALL GETLST(ICORE(I010),1,WDIS,2,IRREPDO,LISTW)
C
C NOW EXPAND W(E<F,A<G) TO W(E<F,AG)
C
         CALL SYMEXP(IRREPDO,VRT(1,ISPIN),WDSZ,ICORE(I010))
C
C EVALUATE V CONTRIBUTION WITH THE MATRIX MULTIPLICATION
C
C             Z(E<F,AM) =  W(E<FA,G)*T(M,G)   [ISPIN=1]
C
C             Z(e<f,am) =  W(e<fa,g)*T(m,g)   [ISPIN=2]
C
C DO PREAMBLE AND DO THE CONTRACTION
C
         IOFFZ=I000
         IOFFW=I010
         DO 100 IRREPG=1,NIRREP
          IRREPA=DIRPRD(IRREPG,IRREPDO)
          IRREPM=IRREPG
          NROWZ=WDSZ*VRT(IRREPA,ISPIN)
          NCOLZ=POP(IRREPM,ISPIN)
          NROWW=NROWZ
          NCOLW=VRT(IRREPG,ISPIN)
          NROWT=VRT(IRREPG,ISPIN)
          NCOLT=POP(IRREPM,ISPIN)
          IOFFT=IOFFT1(IRREPG,ISPIN)
          IF(MIN(NROWZ,NCOLZ,NCOLW).GT.0)THEN
           CALL XGEMM('N','N',NROWZ,NCOLZ,NCOLW,ALPHA,ICORE(IOFFW),
     &                NROWW,ICORE(IOFFT),NROWT,BETA,ICORE(IOFFZ),
     &                NROWZ)
          ENDIF
          IOFFW=IOFFW+IINTFP*NROWW*NCOLW
          IOFFZ=IOFFZ+IINTFP*NROWZ*NCOLZ
100      CONTINUE
C
C WE HAVE TO DO IT OUT-OF-CORE
C
        ELSE
         CORLFT=MXCOR-I010
         CALL W5OAA(ICORE(I010),CORLFT,ICORE(IOFFT1(1,ISPIN)),
     &              ICORE(I000),IRREPDO,LISTW,ISPIN,WDSZ)
        ENDIF
        ENDIF
C
C FINALLY, INCREMENT THE LIST AND EXIT.
C
        I020=I010+IINTFP*TARDSZ*TARDIS
        IF(I020.LE.MAXCOR)THEN
         CALL GETLST(ICORE(I010),1,TARDIS,2,IRREPDO,LSTTAR)
         CALL SAXPY(TARSIZ,ONE,ICORE(I010),1,ICORE(I000),1)
        ELSE
         IOFF=I000
         DO 110 IDIS=1,TARDIS
          CALL GETLST(ICORE(I010),IDIS,1,2,IRREPDO,LSTTAR)
          CALL SAXPY (TARDSZ,ONE,ICORE(I010),1,ICORE(IOFF),1)
          IOFF=IOFF+TARDSZ*IINTFP
110      CONTINUE
        ENDIF
        CALL PUTLST(ICORE(I000),1,TARDIS,2,IRREPDO,LSTTAR)
C
20     CONTINUE
10    CONTINUE
      RETURN
      END
