      SUBROUTINE FORMQ2R(IRREPX,ICORE,MAXCOR)
C
C THIS ROUTINE FORMS THE FOLLOWING CONTRIBUTIONS TO THE THREE-BODY
C INTERMEDIATES FOR RHF REFERENCE FUNCTIONS:
C
C       Q(BE) = SUM C(MN,BF) * [2 <MN|EF> - <NM|EF>]
C
C       Q(MJ) = SUM C(EF,NJ) * [2 <EF|NM> - <FE|NM>]
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      DIMENSION ICORE(MAXCOR)
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMLOC/ISYMOFF(8,8,25)
C
      DATA ONE,ONEM,ZILCH /1.0D0,-1.0D0,0.0D0/
C
C ALLOCATE MEMORY FOR SINGLES (C) ARRAY AND TARGET ARRAYS
C
      LENQAB=IRPDPD(IRREPX,19)
      LENQIJ=IRPDPD(IRREPX,21)
      I0QAB=1
      I0QIJ=I0QAB+IINTFP*LENQAB
      I000 =I0QIJ+IINTFP*LENQIJ
      MXCOR=MAXCOR-I000+1
C
      LISTW=16
      LISTC=446
C
      CALL ZERO(ICORE(I0QAB),LENQAB)
      CALL ZERO(ICORE(I0QIJ),LENQIJ)
C
      DO 5 IRREPW=1,NIRREP
       DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW))
       NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW))
       MAXW=MAX(DISSYW,NUMDSW)
       I010=I000+IINTFP*DISSYW*NUMDSW
       CALL GETLST(ICORE(I000),1,NUMDSW,1,IRREPW,LISTW)
C
C SPIN-ADAPT INTEGRALS
C
       ITMP1=I010
       ITMP2=ITMP1+IINTFP*MAXW
       CALL SPINAD3(IRREPW,VRT(1,1),DISSYW,NUMDSW,ICORE(I000),
     &              ICORE(ITMP1),ICORE(ITMP2)) 
C
C FIRST EVALUATE Q(MJ) :
C
       IRREPCL=IRREPW
       IRREPCR=DIRPRD(IRREPCL,IRREPX)
       DISSYC=IRPDPD(IRREPCL,ISYTYP(1,LISTC))
       NUMDSC=IRPDPD(IRREPCR,ISYTYP(2,LISTC))
       I020=I010+IINTFP*DISSYC*NUMDSC
C
       CALL GETLST(ICORE(I010),1,NUMDSC,1,IRREPCR,LISTC)
C
C PERFORM CONTRACTION AS :
C
C                +
C        W(EFN,M) * C(EFN,J)
C
       IOFFQ=I0QIJ
       IOFFW0=I000
       IOFFC0=I010
       DO 10 IRREPJ=1,NIRREP
        IRREPM=DIRPRD(IRREPJ,IRREPX)
        IRREPN=DIRPRD(IRREPM,IRREPW)
        NUMM=POP(IRREPM,2)
        NUMJ=POP(IRREPJ,2)
        NUMN=POP(IRREPN,1)
        NROW=NUMM
        NCOL=NUMJ
        NSUM=DISSYW*NUMN
        IOFFW=IOFFW0+IINTFP*DISSYW*(ISYMOFF(IRREPM,IRREPW,21)-1)
        IOFFC=IOFFC0+IINTFP*DISSYC*(ISYMOFF(IRREPJ,IRREPCR,21)-1)
        CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFW),NSUM,
     &             ICORE(IOFFC),NSUM,ONE,ICORE(IOFFQ),NROW)
        IOFFQ=IOFFQ+NUMM*NUMJ*IINTFP
10     CONTINUE
C
C NOW EVALUATE Q(EB):
C
       CALL TRANSP(ICORE(I000),ICORE(I010),NUMDSW,DISSYW)
       CALL SCOPY (NUMDSW*DISSYW,ICORE(I010),1,ICORE(I000),1)
       IRREPCR=IRREPW
       IRREPCL=DIRPRD(IRREPCR,IRREPX)
       DISSYC=IRPDPD(IRREPCL,ISYTYP(1,LISTC))
       NUMDSC=IRPDPD(IRREPCR,ISYTYP(2,LISTC))
       I020=I010+IINTFP*DISSYC*NUMDSC
C
       CALL GETLST(ICORE(I010),1,NUMDSC,1,IRREPCR,LISTC)
       CALL TRANSP(ICORE(I010),ICORE(I020),NUMDSC,DISSYC)
       CALL SCOPY (NUMDSC*DISSYC,ICORE(I020),1,ICORE(I010),1)
C
C PERFORM CONTRACTION AS :
C
C                +
C        W(MN,FE) * C(MN,FB)
C
       IOFFQ=I0QAB
       IOFFW0=I000
       IOFFC0=I010
       DO 20 IRREPB=1,NIRREP
        IRREPE=DIRPRD(IRREPB,IRREPX)
        IRREPF=DIRPRD(IRREPE,IRREPW)
        NUME=VRT(IRREPE,2)
        NUMB=VRT(IRREPB,2)
        NUMF=VRT(IRREPF,1)
        NROW=NUME
        NCOL=NUMB
        NSUM=NUMDSW*NUMF
        IOFFW=IOFFW0+IINTFP*NUMDSW*(ISYMOFF(IRREPE,IRREPW,19)-1)
        IOFFC=IOFFC0+IINTFP*NUMDSC*(ISYMOFF(IRREPB,IRREPCL,19)-1)
        CALL XGEMM('T','N',NROW,NCOL,NSUM,ONEM,ICORE(IOFFW),NSUM,
     &             ICORE(IOFFC),NSUM,ONE,ICORE(IOFFQ),NROW)
        IOFFQ=IOFFQ+NUMB*NUME*IINTFP
20     CONTINUE
C
5     CONTINUE
C
      CALL GETLST(ICORE(I000),1,1,1,1,491)
      CALL SAXPY (LENQIJ,ONE,ICORE(I0QIJ),1,ICORE(I000),1)
      CALL PUTLST(ICORE(I000),1,1,1,1,491)
      CALL GETLST(ICORE(I000),1,1,1,1,492)
      CALL SAXPY (LENQAB,ONE,ICORE(I0QAB),1,ICORE(I000),1)
      CALL PUTLST(ICORE(I000),1,1,1,1,492)
C
      RETURN
      END