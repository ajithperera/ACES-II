      SUBROUTINE W4AA1A(ICORE,MAXCOR,IUHF,IOFFLIST)
C
C THIS ROUTINE CALCULATES ONE FIVE CONTRACTIONS WHICH CONTRIBUTES
C  TO THE W(IFMN) INTERMEDIATE USED IN CCSD GRADIENTS AND CCSDT
C  MODELS FOR AAAA AND BBBB SPIN CASES. 
C
C           Z(ifmn) = T1(n,e) * W (mefi)
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,ALPHA,BETA
      CHARACTER*4 SPCASE(2)
      DIMENSION ICORE(MAXCOR),MNFULL(8),ABFULL(8),IOFFT1(8,2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      DATA ONE   /1.0/
      DATA ONEM  /-1.0/
      DATA ZILCH /0.0/
      DATA SPCASE /'AAAA','BBBB'/
      CALL GETT1(ICORE,MAXCOR,MXCOR,IUHF,IOFFT1)
      DO 10 ISPIN=1,1+IUHF
       DO 5 IRREP=1,NIRREP
        ABFULL(IRREP)=IRPDPD(IRREP,18+ISPIN)
        MNFULL(IRREP)=IRPDPD(IRREP,20+ISPIN)
5      CONTINUE
C
C                                
C           Z(ifmn) = - T1(n,e) * W (mefi)
C                                
C           Z(IFMN) = - T1(N,E) * W (MEFI)    [ISPIN = 1]
C
C           Z(ifmn) = - T1(n,e) * W (mefi)    [ISPIN = 2]
C
C HERE THE W IS SIMILAR TO THE W(MBEJ) INTERMEDIATE AND IS STORED AS
C  I(EI,FM) [I(ei,fm)] ON LISTS 53 AND 54 RESPECTIVELY.
C
C READ THE FIRST SET OF INTERMEDIATES INTO CORE FOR
C  ALL IRREPS I(EI,FM) [I(ei,fm)]
C
C
C ALPHA IS ONE SINCE WE STORE -W(MBEJ) ON DISK
C
       ALPHA=ONE
       BETA=ZILCH
       LISTW =122+ISPIN
       LSTTAR=106+ISPIN+IOFFLIST
       ISIZE=ISYMSZ(ISYTYP(1,LISTW),ISYTYP(2,LISTW))
       I000=1
       I010=I000+IINTFP*ISIZE
       I020=I010+IINTFP*ISIZE
       CALL GETALL(ICORE(I010),ISIZE,1,LISTW)
C
C RESORT TO I(EM,FI) [I(em,fi)]
C
       CALL SSTRNG(ICORE(I010),ICORE(I000),ISIZE,ISIZE,ICORE(I020),
     &              SPCASE(ISPIN))
C
C NOW LOOP OVER IRREPS AND FORM PRODUCT
C
       IOFFIR=I000
       DO 20 IRREPDO=1,NIRREP
        WDSZ=IRPDPD(IRREPDO,ISYTYP(1,18+ISPIN))
        WDIS=WDSZ
        TARDSZ=IRPDPD(IRREPDO,ISYTYP(1,LSTTAR))
        TARDIS=IRPDPD(IRREPDO,ISYTYP(2,LSTTAR))
        TARSIZ=TARDSZ*TARDIS
C
C DO IN-PLACE TRANSPOSITION TO FORM I(FI,EM) [I(fi,em)]
C
        CALL MTRAN2(ICORE(IOFFIR),WDIS)
C
C TRANSPOSE BRA AND KET INDICES TO GIVE I(IF,ME)
C
        I020=I010+IINTFP*WDIS
        I030=I020+IINTFP*WDIS
        I040=I030+IINTFP*WDIS
        CALL SYMTR1(IRREPDO,VRT(1,ISPIN),POP(1,ISPIN),WDIS,
     &              ICORE(IOFFIR),ICORE(I010),ICORE(I020),ICORE(I030))
        CALL SYMTR3(IRREPDO,VRT(1,ISPIN),POP(1,ISPIN),WDIS,WDIS,
     &              ICORE(IOFFIR),ICORE(I010),ICORE(I020),ICORE(I030))
C
C EVALUATE V CONTRIBUTION WITH THE MATRIX MULTIPLICATION
C                                          
C             Z(IF,MN) = W(IF,ME) * T(E,N)
C
C
        I020=I010+IINTFP*MNFULL(IRREPDO)*TARDIS
        IOFFZ=I010
        CALL IZERO(ICORE(IOFFZ),IINTFP*MNFULL(IRREPDO)*TARDIS)
        DO 130 IRREPE=1,NIRREP
         IRREPM=DIRPRD(IRREPE,IRREPDO)
         IRREPN=IRREPE
         NROWI=WDSZ*POP(IRREPM,ISPIN)
         NCOLI=VRT(IRREPE,ISPIN)
         NROWZ=NROWI
         NCOLZ=POP(IRREPN,ISPIN)
         NROWT=VRT(IRREPE,ISPIN)
         NCOLT=POP(IRREPN,ISPIN)
         IOFFT=IOFFT1(IRREPE,ISPIN)
         IF(MIN(NROWZ,NCOLZ,NCOLI).GT.0)THEN
          CALL XGEMM('N','N',NROWZ,NCOLZ,NCOLI,ALPHA,ICORE(IOFFIR),
     &               NROWI,ICORE(IOFFT),NROWT,BETA,ICORE(IOFFZ),
     &               NROWZ)
         ENDIF
         IOFFIR=IOFFIR+IINTFP*NROWI*NCOLI
         IOFFZ=IOFFZ+IINTFP*NROWZ*NCOLZ
130     CONTINUE
C
C THE PRODUCT MUST BE CONVERTED FROM Z(IF,MN) TO Z(M<N,IF).  DO 
C  THE REQUIRED OPERATIONS.
C
        I030=I020+IINTFP*TARDIS*TARDSZ
        CALL ASSYM(IRREPDO,POP(1,ISPIN),TARDIS,TARDIS,
     &             ICORE(I020),ICORE(I010))
c YAU : old
c       CALL ICOPY(IINTFP*TARSIZ,ICORE(I020),1,ICORE(I010),1)
c YAU : new
        CALL DCOPY(TARSIZ,ICORE(I020),1,ICORE(I010),1)
c YAU : end
        I030=I020+IINTFP*TARDIS*TARDSZ
        CALL TRANSP(ICORE(I010),ICORE(I020),TARDSZ,TARDIS)
C
C NOW INCREMENT THE STUFF ON THE DISK
C
        CALL GETLST(ICORE(I010),1,TARDIS,1,IRREPDO,LSTTAR)
        CALL SAXPY (TARSIZ,ONE,ICORE(I020),1,ICORE(I010),1)
        CALL PUTLST(ICORE(I010),1,TARDIS,1,IRREPDO,LSTTAR)
20     CONTINUE
10    CONTINUE
      RETURN
      END