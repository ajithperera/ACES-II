      SUBROUTINE W5AB1A(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE CALCULATES ONE OF THE TERMS WHICH CONTRIBUTES
C  TO THE W(eFaM) INTERMEDIATE WHICH IS NEEDED FOR CCSD GRADIENTS
C  AND CCSDT MODELS.  
C
C
C    Z(FeMa) = - SUM T1(n,e) * W (nFaM) - SUM T1(N,F) * W (eNaM)
C                 n                        N
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,ALPHA,BETA
      DIMENSION ICORE(MAXCOR),IOFFT1(8,2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      DATA ONE   /1.0/
      DATA ONEM  /-1.0/
      DATA ZILCH /0.0/
      CALL GETT1(ICORE,MAXCOR,MXCOR,IUHF,IOFFT1)
      LSTTAR=129
      ALPHA=ONEM
      BETA =ZILCH
C
C FIRST CONTRACTION: 
C                                
C    Z(FeMa) = - SUM T1(n,e) * W1(nFaM) 
C                 e                    
C
C HERE W1 IS SIMILAR TO THE W(MBEJ) INTERMEDIATE AND IS STORED
C  ON LIST 117, BUT UNFORTUNATELY AS an,FM
C
C FIRST PIECE 
C
C READ THE TRANSPOSE OF THE FIRST SET OF INTERMEDIATES INTO CORE FOR
C  ALL IRREPS I(an,FM)
C
      MAXTAR=0
      DO 5 IRREP=1,NIRREP
       MAXTAR=MAX(MAXTAR,IRPDPD(IRREP,ISYTYP(1,LSTTAR)))
5     CONTINUE
C
      LISTW =117
      ISIZE=ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
      I000=1
      I010=I000+ISIZE*IINTFP
      I020=I010+MAX(MAXTAR,ISIZE)*IINTFP
      CALL GETALL(ICORE(I010),ISIZE,1,LISTW)
C
C RESORT TO I(aM,Fn).
C
      CALL SSTRNG(ICORE(I010),ICORE(I000),ISIZE,ISIZE,ICORE(I020),
     &            'BBAA')
C
C NOW LOOP OVER IRREPS AND FORM PRODUCT
C
      IOFFIR=I000
      DO 10 IRREPDO=1,NIRREP
       WDSZ=IRPDPD(IRREPDO,ISYTYP(1,22))
       WDIS=IRPDPD(IRREPDO,ISYTYP(2,22))
       TARDSZ=IRPDPD(IRREPDO,ISYTYP(1,LSTTAR))
       TARDIS=IRPDPD(IRREPDO,ISYTYP(2,LSTTAR))
       TARSIZ=TARDIS*TARDSZ
C
C NOW FORM THE PRODUCT MATRIX [Z(aM,Fe)] WITH MATRIX MULTIPLICATION
C                                        +
C            Z(aM,Fe) = I(aMF,n) * T(e,n)
C
       I020=I010+IINTFP*TARDIS*TARDSZ
       IOFFZ=I010
       CALL IZERO(ICORE(IOFFZ),IINTFP*TARDIS*TARDSZ)
       DO 20 IRREPN=1,NIRREP
        IRREPF=DIRPRD(IRREPN,IRREPDO)
        IRREPE=IRREPN
        IOFFT=IOFFT1(IRREPE,2) 
        NROWI=WDSZ*VRT(IRREPF,1)
        NCOLI=POP(IRREPN,2)
        NROWZ=NROWI
        NCOLZ=VRT(IRREPE,2)
        NROWT=VRT(IRREPE,2)
        NCOLT=POP(IRREPN,2)
        NSUM =NCOLI
        IF(MIN(NROWZ,NCOLZ,NSUM).GT.0)THEN
         CALL XGEMM('N','T',NROWZ,NCOLZ,NSUM,ALPHA,ICORE(IOFFIR),
     &              NROWI,ICORE(IOFFT),NROWT,BETA,ICORE(IOFFZ),
     &              NROWZ)
        ENDIF
        IOFFZ=IOFFZ+IINTFP*NROWZ*NCOLZ
        IOFFIR=IOFFIR+IINTFP*NROWI*NCOLI
20     CONTINUE
C
C NOW CONVERT TARGET TO Z(Fe,Ma) AND AUGMENT STUFF ON DISK
C
       I030=I020+IINTFP*TARDSZ
       I040=I030+IINTFP*TARDSZ
       I050=I040+IINTFP*TARDSZ
       CALL SYMTR3(IRREPDO,VRT(1,2),POP(1,1),TARDIS,TARDSZ,
     &             ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040))
       I030=I020+IINTFP*TARSIZ
       IF(I030.LE.MAXCOR)THEN
        CALL TRANSP(ICORE(I010),ICORE(I020),TARDSZ,TARDIS)
c YAU : old
c       CALL ICOPY(IINTFP*TARDSZ*TARDIS,ICORE(I020),1,ICORE(I010),1)
c YAU : new
        CALL DCOPY(TARDSZ*TARDIS,ICORE(I020),1,ICORE(I010),1)
c YAU : end
        CALL GETLST(ICORE(I020),1,TARDIS,1,IRREPDO,LSTTAR)
        CALL SAXPY (TARSIZ,ONE,ICORE(I020),1,ICORE(I010),1)
        CALL PUTLST(ICORE(I010),1,TARDIS,1,IRREPDO,LSTTAR)
       ELSE
        I030=I020+IINTFP*TARDSZ
        IF(I030.GT.MAXCOR)THEN
         CALL INSMEM('W5AB2A',I030,MAXCOR)
         CALL ERREX
        ENDIF
        IOFF=I010
        DO 200 IDIS=1,TARDIS
         CALL SCOPY(TARDSZ,ICORE(IOFF),TARDIS,ICORE(I020),1)
         CALL GETLST(ICORE(I000),IDIS,1,1,IRREPDO,LSTTAR)
         CALL SAXPY (TARDSZ,ONE,ICORE(I020),1,ICORE(I000),1)
         CALL PUTLST(ICORE(I000),IDIS,1,1,IRREPDO,LSTTAR)
         IOFF=IOFF+IINTFP
200     CONTINUE
       ENDIF
10    CONTINUE 
C
C SECOND PIECE
C
C    Z(FeMa) =  - SUM T1(N,F) * W (eNaM)
C                  N
C
C READ IN THE I(aN,eM) INTERMEDIATES
C
      LISTW =126
      ISIZE=ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
      I000=1
      I010=I000+ISIZE*IINTFP
      I020=I010+ISIZE*IINTFP
      CALL GETALL(ICORE(I010),ISIZE,1,LISTW)
C
C RESORT TO I(aM,eN).
C
      CALL SSTRNG(ICORE(I010),ICORE(I000),ISIZE,ISIZE,ICORE(I020),
     &            'BABA')
C
C NOW LOOP OVER IRREPS AND FORM PRODUCT
C
      IOFFIR=I000
      DO 110 IRREPDO=1,NIRREP
       WDSZ=IRPDPD(IRREPDO,ISYTYP(1,26))
       WDIS=WDSZ
       TARDSZ=IRPDPD(IRREPDO,ISYTYP(1,LSTTAR))
       TARDIS=IRPDPD(IRREPDO,ISYTYP(2,LSTTAR))
       TARSIZ=TARDIS*TARDSZ
C
C NOW FORM THE PRODUCT MATRIX [Z(aM,eF)] WITH MATRIX MULTIPLICATION
C                                        +
C            Z(aM,eF) = I(aMe,N) * T(F,N)
C
       I020=I010+IINTFP*TARDIS*TARDSZ
       IOFFZ=I010
       CALL IZERO(ICORE(IOFFZ),IINTFP*TARDIS*TARDSZ)
       DO 120 IRREPN=1,NIRREP
        IRREPE=DIRPRD(IRREPN,IRREPDO)
        IRREPF=IRREPN
        IOFFT=IOFFT1(IRREPF,1) 
        NROWI=WDSZ*VRT(IRREPE,2)
        NCOLI=POP(IRREPN,1)
        NROWZ=NROWI
        NCOLZ=VRT(IRREPF,1)
        NROWT=VRT(IRREPF,1)
        NCOLT=POP(IRREPN,1)
        NSUM =NCOLI
        IF(MIN(NROWZ,NCOLZ,NSUM).GT.0)THEN
         CALL XGEMM('N','T',NROWZ,NCOLZ,NSUM,ALPHA,ICORE(IOFFIR),
     &              NROWI,ICORE(IOFFT),NROWT,BETA,ICORE(IOFFZ),
     &              NROWZ)
        ENDIF
        IOFFZ=IOFFZ+IINTFP*NROWZ*NCOLZ
        IOFFIR=IOFFIR+IINTFP*NROWI*NCOLI
120    CONTINUE
C
C NOW REORDER TARGET MATRIX TO Z(Fe,Ma)
C
C  Z(aM,eF) -> Z(Ma,eF) -> Z(Ma,Fe) -> Z(Fe,Ma)
C
       I030=I020+IINTFP*MAX(TARDIS,TARDSZ)
       I040=I030+IINTFP*MAX(TARDIS,TARDSZ)
       I050=I040+IINTFP*MAX(TARDIS,TARDSZ)
       CALL SYMTR1(IRREPDO,VRT(1,2),VRT(1,1),TARDIS,ICORE(I010),
     &             ICORE(I020),ICORE(I030),ICORE(I040))
       CALL SYMTR3(IRREPDO,VRT(1,2),POP(1,1),TARDIS,TARDSZ,
     &             ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040))
       I030=I020+IINTFP*TARDIS*TARDSZ
       IF(I030.LE.MAXCOR)THEN
        CALL TRANSP(ICORE(I010),ICORE(I020),TARDSZ,TARDIS)
        CALL GETLST(ICORE(I010),1,TARDIS,1,IRREPDO,LSTTAR)
        CALL SAXPY (TARSIZ,ONE,ICORE(I020),1,ICORE(I010),1)
        CALL PUTLST(ICORE(I010),1,TARDIS,1,IRREPDO,LSTTAR)
       ELSE
        I030=I020+IINTFP*TARDSZ
        IF(I030.GT.MAXCOR)THEN
         CALL INSMEM('W5AB2A',I030,MAXCOR)
         CALL ERREX
        ENDIF
        IOFF=I010
        DO 201 IDIS=1,TARDIS
         CALL SCOPY(TARDSZ,ICORE(IOFF),TARDIS,ICORE(I020),1)
         CALL GETLST(ICORE(I000),IDIS,1,1,IRREPDO,LSTTAR)
         CALL SAXPY (TARDSZ,ONE,ICORE(I020),1,ICORE(I000),1)
         CALL PUTLST(ICORE(I000),IDIS,1,1,IRREPDO,LSTTAR)
         IOFF=IOFF+IINTFP
201     CONTINUE
       ENDIF
110   CONTINUE 
      RETURN
      END
