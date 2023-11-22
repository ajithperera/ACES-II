      SUBROUTINE W4AB1A(ICORE,MAXCOR,IUHF,IOFFLIST)
C
C THIS ROUTINE CALCULATES ONE OF THE TERMS WHICH CONTRIBUTES
C  TO THE W(MnFi) INTERMEDIATE WHICH IS NEEDED FOR CCSD GRADIENTS
C  AND CCSDT MODELS.  
C
C    Z(MnFi) = SUM T1(n,e) * W (MeFi) + SUM T1(M,E) * W (EnFi)
C               e                        E
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
      LSTTAR=109+IOFFLIST
C
C ALPHA MUST BE ONE SINCE THE NEGATIVE OF THE W's ARE STORED
C
      ALPHA=ONE
      BETA =ZILCH
C
C FIRST CONTRACTION: 
C                                
C    Z(MnFi) =  SUM T1(n,e) * W1(MeFi) 
C                e                    
C
C HERE W1 AND W2 ARE SIMILAR TO THE W(MBEJ) INTERMEDIATE AND ARE STORED
C  ON LISTS 57 AND 58, BUT UNFORTUNATELY AS ei,FM AND Ei,Fn
C
C FIRST PIECE 
C
C READ THE TRANSPOSE OF THE FIRST SET OF INTERMEDIATES INTO CORE FOR
C  ALL IRREPS I(FM,ei)
C
      LISTW =117
      ISIZE=ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
      I000=1
      I010=I000+ISIZE*IINTFP
      I020=I010+ISIZE*IINTFP
      CALL GETALT(ICORE(I010),ISIZE,1,LISTW)
C
C RESORT TO I(Fi,eM).
C
      CALL SSTRNG(ICORE(I010),ICORE(I000),ISIZE,ISIZE,ICORE(I020),
     &            'AABB')
C
C NOW LOOP OVER IRREPS AND FORM PRODUCT
C
      IOFFIR=I000
      DO 10 IRREPDO=1,NIRREP
       WDSZ=IRPDPD(IRREPDO,ISYTYP(1,21))
       WDIS=IRPDPD(IRREPDO,ISYTYP(2,21))
       TARDSZ=IRPDPD(IRREPDO,ISYTYP(1,LSTTAR))
       TARDIS=IRPDPD(IRREPDO,ISYTYP(2,LSTTAR))
       TARSIZ=TARDIS*TARDSZ
C
C TRANSPOSE KET INDICES TO GIVE I(Fi,Me)
C
       I020=I010+IINTFP*WDSZ       
       I030=I020+IINTFP*WDSZ       
       I040=I030+IINTFP*WDSZ       
       CALL SYMTR1(IRREPDO,VRT(1,2),POP(1,1),WDSZ,ICORE(IOFFIR),
     &             ICORE(I010),ICORE(I020),ICORE(I030))
C
C NOW FORM THE PRODUCT MATRIX [Z(Fi,Mn)] WITH MATRIX MULTIPLICATION
C
C            Z(Fi,Mn) = I(Fi,Me) * T(e,n)
C
       I020=I010+IINTFP*TARDIS*TARDSZ
       IOFFZ=I010
       CALL IZERO(ICORE(IOFFZ),IINTFP*TARDIS*TARDSZ)
       DO 20 IRREPE=1,NIRREP
        IRREPM=DIRPRD(IRREPE,IRREPDO)
        IRREPN=IRREPE
        IOFFT=IOFFT1(IRREPE,2) 
        NROWI=WDSZ*POP(IRREPM,1)
        NCOLI=VRT(IRREPE,2)
        NROWZ=NROWI
        NCOLZ=POP(IRREPN,2)
        NROWT=VRT(IRREPE,2)
        NCOLT=POP(IRREPN,2)
        NSUM =NCOLI
        IF(MIN(NROWZ,NCOLZ,NSUM).GT.0)THEN
         CALL XGEMM('N','N',NROWZ,NCOLZ,NSUM,ALPHA,ICORE(IOFFIR),
     &              NROWI,ICORE(IOFFT),NROWT,BETA,ICORE(IOFFZ),
     &              NROWZ)
        ENDIF
        IOFFZ=IOFFZ+IINTFP*NROWZ*NCOLZ
        IOFFIR=IOFFIR+IINTFP*NROWI*NCOLI
20     CONTINUE
C
C NOW TRANSPOSE TARGET AND INCREMENT W(Mn,Fi) LIST ON DISK
C
       I030=I020+IINTFP*TARSIZ
       CALL TRANSP(ICORE(I010),ICORE(I020),TARDSZ,TARDIS)
       CALL GETLST(ICORE(I010),1,TARDIS,1,IRREPDO,LSTTAR)
       CALL SAXPY (TARSIZ,ONE,ICORE(I020),1,ICORE(I010),1)
       CALL PUTLST(ICORE(I010),1,TARDIS,1,IRREPDO,LSTTAR)
10    CONTINUE 
C
C SECOND PIECE
C    Z(MnFi) =  SUM T1(M,E) * W (EnFi)
C                E
C
C READ IN FIRST SET OF INTERMEDIATES INTO CORE FOR
C  ALL IRREPS I(Ei,Fn)
C
      LISTW =125
      ISIZE=ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
      I000=1
      I010=I000+ISIZE*IINTFP
      I020=I010+ISIZE*IINTFP
      CALL GETALL(ICORE(I010),ISIZE,1,LISTW)
C
C RESORT TO I(En,Fi).
C
      CALL SSTRNG(ICORE(I010),ICORE(I000),ISIZE,ISIZE,ICORE(I020),
     &            'ABAB')
C
C NOW LOOP OVER IRREPS AND FORM PRODUCT
C
      IOFFIR=I000
      DO 110 IRREPDO=1,NIRREP
       WDSZ=IRPDPD(IRREPDO,ISYTYP(1,25))
       WDIS=WDSZ
       TARDSZ=IRPDPD(IRREPDO,ISYTYP(1,LSTTAR))
       TARDIS=IRPDPD(IRREPDO,ISYTYP(2,LSTTAR))
       TARSIZ=TARDIS*TARDSZ
C
C DO IN-PLACE TRANSPOSITION OF I MATRIX AND THEN SWAP KET INDICES,
C  GIVING I(Fi,nE).
C
       CALL MTRAN2(ICORE(IOFFIR),WDIS)
       I020=I010+IINTFP*WDSZ       
       I030=I020+IINTFP*WDSZ       
       I040=I030+IINTFP*WDSZ       
       CALL SYMTR1(IRREPDO,VRT(1,1),POP(1,2),WDSZ,ICORE(IOFFIR),
     &             ICORE(I010),ICORE(I020),ICORE(I030))
C
C NOW FORM THE PRODUCT MATRIX [Z(Fi,nM)] WITH MATRIX MULTIPLICATION
C
C            Z(Fi,nM) = I(Fi,nE) * T(E,M)
C
       I020=I010+IINTFP*TARDIS*TARDSZ
       IOFFZ=I010
       CALL IZERO(ICORE(IOFFZ),IINTFP*TARDIS*TARDSZ)
       DO 120 IRREPE=1,NIRREP
        IRREPN=DIRPRD(IRREPE,IRREPDO)
        IRREPM=IRREPE
        IOFFT=IOFFT1(IRREPE,1) 
        NROWI=WDSZ*POP(IRREPN,2)
        NCOLI=VRT(IRREPE,1)
        NROWZ=NROWI
        NCOLZ=POP(IRREPM,1)
        NROWT=VRT(IRREPE,1)
        NCOLT=POP(IRREPM,1)
        NSUM =NCOLI
        IF(MIN(NROWZ,NCOLZ,NSUM).GT.0)THEN
         CALL XGEMM('N','N',NROWZ,NCOLZ,NSUM,ALPHA,ICORE(IOFFIR),
     &              NROWI,ICORE(IOFFT),NROWT,BETA,ICORE(IOFFZ),
     &              NROWZ)
        ENDIF
        IOFFZ=IOFFZ+IINTFP*NROWZ*NCOLZ
        IOFFIR=IOFFIR+IINTFP*NROWI*NCOLI
120    CONTINUE
C
C NOW SWAP KET INDICES OF TARGET AND THEN TRANSPOSE IT.
C
C  Z(Fi,nM) -> Z(Fi,Mn) -> Z(Mn,Fi)
C
       I030=I020+IINTFP*TARDIS
       I040=I030+IINTFP*TARDIS
       I050=I040+IINTFP*TARDIS
       CALL SYMTR1(IRREPDO,POP(1,2),POP(1,1),TARDIS,ICORE(I010),
     &             ICORE(I020),ICORE(I030),ICORE(I040))
       I030=I020+IINTFP*TARDIS*TARDSZ
       CALL TRANSP(ICORE(I010),ICORE(I020),TARDSZ,TARDIS)
C
C NOW INCREMENT LISTS ON DISK
C
       CALL GETLST(ICORE(I010),1,TARDIS,1,IRREPDO,LSTTAR)
       CALL SAXPY (TARSIZ,ONE,ICORE(I020),1,ICORE(I010),1)
       CALL PUTLST(ICORE(I010),1,TARDIS,1,IRREPDO,LSTTAR)
110   CONTINUE 
      RETURN
      END