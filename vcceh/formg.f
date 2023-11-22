      SUBROUTINE FORMG(IRREPX,ICORE,MAXCOR)
C
C THIS ROUTINE FORMS THE FOLLOWING CONTRIBUTIONS TO THE THREE-BODY
C G INTERMEDIATES FOR RHF REFERENCE FUNCTIONS:
C
C       G(BE) = SUM C(MN,BF) * [2 T(MN,EF) - T(NM,EF)]
C
C       G(MJ) = SUM C(EF,NJ) * [2 T(EF,NM) - T(FE,NM)]
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,SNRM2
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
      LENGAB=IRPDPD(IRREPX,19)
      LENGIJ=IRPDPD(IRREPX,21)
      I0GAB=1
      I0GIJ=I0GAB+IINTFP*LENGAB
      I000 =I0GIJ+IINTFP*LENGIJ
      MXCOR=MAXCOR-I000+1
C
      LISTT=46
      LISTC=446
C
      CALL ZERO(ICORE(I0GAB),LENGAB)
      CALL ZERO(ICORE(I0GIJ),LENGIJ)
C
      DO 5 IRREPT=1,NIRREP
       DISSYT=IRPDPD(IRREPT,ISYTYP(1,LISTT))
       NUMDST=IRPDPD(IRREPT,ISYTYP(2,LISTT))
       MAXT=MAX(DISSYT,NUMDST)
       I010=I000+IINTFP*DISSYT*NUMDST
       CALL GETLST(ICORE(I000),1,NUMDST,1,IRREPT,LISTT)
C
C SPIN-ADAPT INTEGRALS
C
       ITMP1=I010
       ITMP2=ITMP1+IINTFP*MAXT
       CALL SPINAD3(IRREPT,VRT(1,1),DISSYT,NUMDST,ICORE(I000),
     &              ICORE(ITMP1),ICORE(ITMP2)) 
C
C FIRST EVALUATE G(MJ) :
C
       IRREPCL=IRREPT
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
C        T(EFN,M) * C(EFN,J)
C
       IOFFG=I0GIJ
       IOFFT0=I000
       IOFFC0=I010
       DO 10 IRREPJ=1,NIRREP
        IRREPM=DIRPRD(IRREPJ,IRREPX)
        IRREPN=DIRPRD(IRREPM,IRREPT)
        NUMM=POP(IRREPM,2)
        NUMJ=POP(IRREPJ,2)
        NUMN=POP(IRREPN,1)
        NROW=NUMM
        NCOL=NUMJ
        NSUM=DISSYT*NUMN
        IOFFT=IOFFT0+IINTFP*DISSYT*(ISYMOFF(IRREPM,IRREPT,21)-1)
        IOFFC=IOFFC0+IINTFP*DISSYC*(ISYMOFF(IRREPJ,IRREPCR,21)-1)
        CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFT),NSUM,
     &             ICORE(IOFFC),NSUM,ONE,ICORE(IOFFG),NROW)
        IOFFG=IOFFG+NUMM*NUMJ*IINTFP
10     CONTINUE
C
C NOW EVALUATE G(EB):
C
       CALL TRANSP(ICORE(I000),ICORE(I010),NUMDST,DISSYT)
       CALL SCOPY (NUMDST*DISSYT,ICORE(I010),1,ICORE(I000),1)
       IRREPCR=IRREPT
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
C        T(MN,FE) * C(MN,FB)
C
       IOFFG=I0GAB
       IOFFT0=I000
       IOFFC0=I010
       DO 20 IRREPB=1,NIRREP
        IRREPE=DIRPRD(IRREPB,IRREPX)
        IRREPF=DIRPRD(IRREPE,IRREPT)
        NUME=VRT(IRREPE,2)
        NUMB=VRT(IRREPB,2)
        NUMF=VRT(IRREPF,1)
        NROW=NUME
        NCOL=NUMB
        NSUM=NUMDST*NUMF
        IOFFT=IOFFT0+IINTFP*NUMDST*(ISYMOFF(IRREPE,IRREPT,19)-1)
        IOFFC=IOFFC0+IINTFP*NUMDSC*(ISYMOFF(IRREPB,IRREPCL,19)-1)
        CALL XGEMM('T','N',NROW,NCOL,NSUM,ONEM,ICORE(IOFFT),NSUM,
     &             ICORE(IOFFC),NSUM,ONE,ICORE(IOFFG),NROW)
        IOFFG=IOFFG+NUMB*NUME*IINTFP
20     CONTINUE
C
5     CONTINUE
C
      CALL PUTLST(ICORE(I0GIJ),1,1,1,1,491)
      CALL PUTLST(ICORE(I0GAB),1,1,1,1,492)
C
      RETURN
      END
