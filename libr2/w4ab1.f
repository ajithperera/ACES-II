      SUBROUTINE W4AB1(ICORE,MAXCOR,IUHF,TERM1,TERM2,TERM3,TAU,
     &                 IOFFLIST)
C
C THIS ROUTINE CALCULATES THREE OF THE FIVE CONTRACTIONS WHICH 
C  CONTRIBUTE TO THE W(iFmN) INTERMEDIATE WHICH IS USED IN CCSD
C  GRADIENTS AND CCSDT MODELS FOR AAAA AND BBBB SPIN CASES.  
C
C  FIRST CONTRACTION:
C           Z(iFmN) = - SUM T1(F,O) * W (iOmN)
C  SECOND CONTRACTION:
C           Z(iFmN) = T2(eF,mN) * F(ie)
C  THIRD CONTRACTION
C           Z(iFmN) = TAU(eG,mN) <iF|eG>
C
C  THE REMAINING CONTRACTIONS ARE DONE BEFORE THIS ROUTINE IS
C   CALLED AND THE TARGET LIST CONTAINS THEM ON ENTRY. 
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,ALPHA,BETA
      LOGICAL INCORE,TERM1,TERM2,TERM3,TAU
      DIMENSION ICORE(MAXCOR),IOFFT1(8,2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      DATA ONE   /1.0D0/
      DATA ONEM  /-1.0D0/
      DATA ZILCH /0.0D0/
      CALL GETT1(ICORE,MAXCOR,MXCOR,IUHF,IOFFT1)
      IBOT=1
      BETA=ONE
      LSTTAR=109+IOFFLIST
      DO 20 IRREPDO=1,NIRREP
       TARDSZ=IRPDPD(IRREPDO,ISYTYP(1,9))
       TARDIS=IRPDPD(IRREPDO,ISYTYP(2,9))
       TARSIZ=TARDSZ*TARDIS
       LISTW =53
       WDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTW))
       WDIS=WDSZ
       ALPHA=ONEM
       CALL IZERO(ICORE,IINTFP*TARSIZ)
       I000=1
       I010=I000+IINTFP*TARSIZ
       IF(TERM1)THEN
C
C FIRST CONTRACTION: 
C                                
C           Z(iFmN) = - T1(F,O) * W (iOmN)
C
       I020=I010+IINTFP*WDSZ*WDIS
       I030=I020+IINTFP*WDSZ
       CALL GETLST(ICORE(I010),1,WDIS,2,IRREPDO,LISTW)
       CALL MTRAN2(ICORE(I010),WDIS)
C
C TRANSPOSE KET INDICES TO GIVE W(Nm,iO)
C
       I030=I020+IINTFP*WDSZ
       I040=I030+IINTFP*WDSZ
       I050=I040+IINTFP*WDSZ
       CALL SYMTR1(IRREPDO,POP(1,1),POP(1,2),WDSZ,ICORE(I010),
     &             ICORE(I020),ICORE(I030),ICORE(I040))
C
C EVALUATE V CONTRIBUTION WITH THE MATRIX MULTIPLICATION
C                                           +
C             Z(Nm,iF) = - W(Nm,iO) * T(F,O)
C
       IOFFZ=I000
       IOFFW=I010
       DO 330 IRREPO=1,NIRREP
        IRREPI=DIRPRD(IRREPO,IRREPDO)
        IRREPF=IRREPO
        NROWZ=WDSZ*POP(IRREPI,2)
        NCOLZ=VRT(IRREPF,1)
        NROWW=NROWZ
        NCOLW=POP(IRREPO,1)
        NROWT=VRT(IRREPF,1)
        NCOLT=NCOLW
        IOFFT=IOFFT1(IRREPO,1)
        IF(MIN(NROWZ,NCOLZ,NCOLW).GT.0)THEN
         CALL XGEMM('N','T',NROWZ,NCOLZ,NCOLW,ALPHA,ICORE(IOFFW),
     &              NROWW,ICORE(IOFFT),NROWT,BETA,ICORE(IOFFZ),
     &              NROWZ)
        ENDIF
        IOFFW=IOFFW+IINTFP*NROWW*NCOLW
        IOFFZ=IOFFZ+IINTFP*NROWZ*NCOLZ
330    CONTINUE
       ENDIF
C
C NOW REARRANGE THIS SO THAT IT IS ORDERED IN THE SAME WAY AS THE NEXT 
C  TWO TERMS (WHICH IS ALSO THE WAY THAT THE TARGET LIST IS ORDERED).
C
C                  Z(Nm,iF) -> Z(Nm,Fi)
C
       I020=I010+IINTFP*TARDSZ
       I030=I020+IINTFP*TARDSZ
       I040=I030+IINTFP*TARDSZ
       CALL SYMTR1(IRREPDO,POP(1,2),VRT(1,1),TARDSZ,ICORE(I000),
     &             ICORE(I010),ICORE(I020),ICORE(I030))
       IF(TERM2)THEN 
C
C SECOND CONTRACTION: 
C                                
C           Z(iFmN) = SUM T2(eF,mN) * F(ie)
C                      e
       ALPHA=ONE
       LISTT2=46
       NT2DSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTT2))
       NT2DIS=IRPDPD(IRREPDO,ISYTYP(2,LISTT2))
       I020=I010+IINTFP*NT2DSZ*NT2DIS
       I030=I020+IINTFP*NT2DSZ
       CALL GETTRN(ICORE(I010),ICORE(I020),NT2DSZ,NT2DIS,1,IRREPDO,
     &             LISTT2)
C
C REDEFINE DISTRIBUTION SIZE AND NUMBER OF DISTRIBUTIONS
C
       ITMP=NT2DIS
       NT2DIS=NT2DSZ
       NT2DSZ=ITMP
C
C WE NOW HAVE T(Nm,Fe) IN CORE.
C
C
C EVALUATE V CONTRIBUTION WITH THE MATRIX MULTIPLICATION
C                                          
C             Z(Nm,Fi) = T(Nm,Fe) * F(e,i)
C
       CALL GETLST(ICORE(I020),1,1,1,2,93)
       IOFFF =I020
       IOFFT2=I010
       IOFFZ =I000
       DO 430 IRREPE=1,NIRREP
        IRREPF=DIRPRD(IRREPE,IRREPDO)
        IRREPI=IRREPE
        NROWT2=NT2DSZ*VRT(IRREPF,1)
        NCOLT2=VRT(IRREPE,2)
        NROWZ =NROWT2
        NCOLZ =POP(IRREPI,2)
        NROWF =NCOLT2
        NCOLF =NCOLZ
        IF(MIN(NROWZ,NCOLZ,NCOLT2).GT.0)THEN
         CALL XGEMM('N','N',NROWZ,NCOLZ,NCOLT2,ALPHA,ICORE(IOFFT2),
     &               NROWT2,ICORE(IOFFF),NROWF,BETA,ICORE(IOFFZ),
     &               NROWZ)
        ENDIF
        IOFFF =IOFFF+IINTFP*NROWF*NCOLF
        IOFFT2=IOFFT2+IINTFP*NROWT2*NCOLT2
        IOFFZ =IOFFZ+IINTFP*NROWZ*NCOLZ
430    CONTINUE
       ENDIF
       IF(TERM3)THEN
C
C THIRD AND FINAL TERM (N**6 DEPENDENCE)
C                                
C           Z(iFmN) = TAU(eG,mN) <iF|eG>
C
        ALPHA=ONE
        LSTINT=30
        LISTT2=46
        INTDSZ=IRPDPD(IRREPDO,ISYTYP(1,LSTINT))
        INTDIS=IRPDPD(IRREPDO,ISYTYP(2,LSTINT))
        TAUDSZ=IRPDPD(IRREPDO,ISYTYP(1,LISTT2))
        TAUDIS=IRPDPD(IRREPDO,ISYTYP(2,LISTT2))
        I020=I010+IINTFP*TAUDSZ*TAUDIS
        CALL GETLST(ICORE(I010),1,TAUDIS,1,IRREPDO,LISTT2)
        IF(TAU)THEN
         CALL FTAU(ICORE(I010),ICORE(IOFFT1(1,1)),ICORE(IOFFT1(1,2)),
     &             TAUDSZ,TAUDIS,POP(1,1),POP(1,2),VRT(1,1),
     &             VRT(1,2),IRREPDO,3,ONE)
        ENDIF
C
C PERFORM CONTRACTION WITH GENERAL ALGORITHM
C
        MXCOR2=MXCOR-I020+1
        NINCOR=MXCOR2/(IINTFP*MAX(1,INTDSZ))
        nincor=1
        IFIRST=1
        NLEFT=INTDIS
        IOFFZ=I000
1       NREAD=MIN(NLEFT,NINCOR)
        CALL GETLST(ICORE(I020),IFIRST,NREAD,1,IRREPDO,LSTINT)
        NROW=TAUDIS
        NCOL=NREAD
        NSUM=TAUDSZ
        CALL XGEMM('T','N',NROW,NCOL,NSUM,ALPHA,ICORE(I010),
     &             NSUM,ICORE(I020),INTDSZ,BETA,ICORE(IOFFZ),
     &             NROW)
        IOFFZ=IOFFZ+IINTFP*NROW*NCOL
        NLEFT=NLEFT-NREAD
        IFIRST=IFIRST+NREAD
        IF(NLEFT.NE.0)GOTO 1      
       ENDIF
C
C FINALLY, INCREMENT THE TARGET LIST WITH Z(Nm,Fi) AND LEAVE.
C
       I020=I010+IINTFP*TARDSZ*TARDIS
       CALL GETLST(ICORE(I010),1,TARDIS,2,IRREPDO,LSTTAR)
       CALL SAXPY(TARSIZ,ONE,ICORE(I010),1,ICORE(I000),1)
       CALL PUTLST(ICORE(I000),1,TARDIS,2,IRREPDO,LSTTAR)
20    CONTINUE
      RETURN
      END
