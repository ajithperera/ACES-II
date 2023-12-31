
      SUBROUTINE MEMSIZ(NAAAA1,NAAAA2,NAAAA3,NAAAA4,NAABB1,NAABB2,
     &                  NAABB3,NABAB1,NABAB2,NABAB3,NABCD1,NABCD2,
     &                  NABCD3,NBFTMP,IUHF,NAAAA11,NAABB11,NABAB11,
     &                  NABCD11)
      IMPLICIT INTEGER (A-Z)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      DIMENSION IDID(8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NF(2),NFMI(2),NFEA(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /AOPOPS/ AOPOP(8),MOPOP(8),NAO,NAOSPH,NMO
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
C
      NNP1O2(I)=(I*(I+1))/2
C
      IF(MBPT2) THEN
C
C COMPUTE MAXIMUM MEMORY REQUIRED FOR BUF1 ARRAY IN BXAAAA
C  
      IF(IUHF.EQ.0) THEN
       NAAAA1=0
       ITMP=IRPDPD(1,ISYTYP(1,18))*IRPDPD(1,ISYTYP(2,18))
       NAAAA1=ITMP      
       NAAAA4=ITMP
       DO 10 IRREP=1,NIRREP
        ITMP=NNP1O2(AOPOP(IRREP))*VRT(IRREP,1)*POP(IRREP,1)
        NAAAA1=MAX(NAAAA1,ITMP)
10     CONTINUE
C
C COMPUTE MAXIMUM MEMORY REQUIRED FOR BUF2 ARRAY IN BXAAAA
C
      NAAAA2=0
      DO 11 IRREP=1,NIRREP
       ITMP=VRT(IRREP,1)*POP(IRREP,1)*VRT(IRREP,1)*POP(IRREP,1)
       NAAAA2=NAAAA2+ITMP
11    CONTINUE
C
      NAAAA3=0
C
      ELSE
       NAAAA1=0
       ITMP1=IRPDPD(1,ISYTYP(1,18))*IRPDPD(1,ISYTYP(2,18))
       ITMP2=IRPDPD(1,ISYTYP(1,19))*IRPDPD(1,ISYTYP(2,19))
       ITMP3=IRPDPD(1,ISYTYP(1,20))*IRPDPD(1,ISYTYP(2,20))
       NAAAA1=MAX(NAAAA1,ITMP1,ITMP2,ITMP3)
       NAAAA4=NAAAA1
       DO 15 IRREP=1,NIRREP
        ITMP1=NNP1O2(AOPOP(IRREP))*VRT(IRREP,1)*POP(IRREP,1)
        ITMP2=NNP1O2(AOPOP(IRREP))*VRT(IRREP,2)*POP(IRREP,2)
        NAAAA1=MAX(NAAAA1,ITMP1,ITMP2)
15     CONTINUE
C
C COMPUTE MAXIMUM MEMORY REQUIRED FOR BUF2 AND BUF3 ARRAYS IN BXAAAA
C
      NAAAA2=0
      NAAAA3=0
      DO 16 IRREP=1,NIRREP
       ITMP1=VRT(IRREP,1)*POP(IRREP,1)*VRT(IRREP,1)*POP(IRREP,1)
       ITMP2=VRT(IRREP,2)*POP(IRREP,2)*VRT(IRREP,2)*POP(IRREP,2)
       ITMP3=VRT(IRREP,1)*POP(IRREP,1)*VRT(IRREP,2)*POP(IRREP,2)
       NAAAA2=NAAAA2+MAX(ITMP1,ITMP3)
       NAAAA3=NAAAA3+ITMP2
16    CONTINUE
      ENDIF
C
C*******************************************************************
C
C COMPUTE MAXIMUM MEMORY REQUIRED FOR BUF1 ARRAY IN BXAABB
C
      IF(IUHF.EQ.0) THEN
       NAABB1=0
       ITMP=IRPDPD(1,ISYTYP(1,18))*IRPDPD(1,ISYTYP(2,18))
       NAABB1=ITMP
       DO 20 IRREP1=2,NIRREP
        DO 21 IRREP2=1,IRREP1-1
         ITMP=NNP1O2(AOPOP(IRREP2))*VRT(IRREP1,1)*POP(IRREP1,1)
         ITMP2=2*AOPOP(IRREP2)*AOPOP(IRREP2)
         NAABB1=MAX(NAABB1,ITMP,ITMP2)
21      CONTINUE
20     CONTINUE
C 
C COMPUTE MAXIMUM MEMORY REQUIRED FOR BUF2 ARRAY IN BXAABB
C
       NAABB2=0
       NAABB3=0
       DO 22 IRREP1=2,NIRREP
        DO 23 IRREP2=1,IRREP1-1
         ITMP=VRT(IRREP2,1)*POP(IRREP2,1)*VRT(IRREP1,1)
     &         *POP(IRREP1,1)
         NAABB2=NAABB2+ITMP
23      CONTINUE
22     CONTINUE
C
      ELSE
       NAABB1=0
       ITMP1=IRPDPD(1,ISYTYP(1,18))*IRPDPD(1,ISYTYP(2,18))
       ITMP2=IRPDPD(1,ISYTYP(1,19))*IRPDPD(1,ISYTYP(2,19))
       ITMP3=IRPDPD(1,ISYTYP(1,20))*IRPDPD(1,ISYTYP(2,20))
       NAABB1=MAX(ITMP1,ITMP2,ITMP3)
       DO 24 IRREP1=2,NIRREP
        DO 25 IRREP2=1,IRREP1-1
         ITMP=NNP1O2(AOPOP(IRREP2))*MAX((VRT(IRREP1,1)
     &               *POP(IRREP1,1)),(VRT(IRREP1,2)*POP(IRREP1,2)))
         ITMP2=2*AOPOP(IRREP2)*AOPOP(IRREP2)
         NAABB1=MAX(NAABB1,ITMP,ITMP2)
25      CONTINUE
24     CONTINUE
C
       ITMP1=0
       ITMP2=0
       ITMP3=0
       ITMP4=0
       DO 26 IRREP1=2,NIRREP
        DO 27 IRREP2=1,IRREP1-1    
         ITMP1=ITMP1+VRT(IRREP2,1)*POP(IRREP2,1)*VRT(IRREP1,2)
     &             *POP(IRREP1,2)
         ITMP2=ITMP2+VRT(IRREP2,2)*POP(IRREP2,2)*VRT(IRREP1,1)
     &             *POP(IRREP1,1)
         ITMP3=ITMP3+VRT(IRREP2,1)*POP(IRREP2,1)*VRT(IRREP1,1)
     &             *POP(IRREP1,1)
         ITMP4=ITMP4+VRT(IRREP2,2)*POP(IRREP2,2)*VRT(IRREP1,2)
     &             *POP(IRREP1,2)
27      CONTINUE
26     CONTINUE
       NAABB2=MAX(ITMP1,ITMP2)
       NAABB3=MAX(ITMP3,ITMP4)
C
      ENDIF
C*******************************************************************
C
C COMPUTE MAXIMUM MEMORY REQUIRED FOR BUF1 ARRAY IN BXABAB
C
      IF(IUHF.EQ.0) THEN
      NABAB1=0
      DO 30 IRREP=1,NIRREP
       ITMP=IRPDPD(IRREP,ISYTYP(1,18))*IRPDPD(IRREP,ISYTYP(2,18))
       NABAB1=MAX(NABAB1,ITMP)
30    CONTINUE
      NABAB11=NABAB1
      DO 31 IRREP1=2,NIRREP
       DO 32 IRREP2=1,NIRREP-1
        ITMP=AOPOP(IRREP1)*AOPOP(IRREP2)
     &      *(VRT(IRREP1,1)*POP(IRREP2,1)
     &      +VRT(IRREP2,1)*POP(IRREP1,1))
        NABAB1=MAX(NABAB1,ITMP)
        ITMP=MAX(AOPOP(IRREP1),AOPOP(IRREP2))
     &      *(VRT(IRREP1,1)*POP(IRREP2,1)
     &      +VRT(IRREP2,1)*POP(IRREP1,1))
        NABAB11=MAX(NABAB11,ITMP)
32     CONTINUE
31    CONTINUE
C
C COMPUTE MAXIMUM MEMORY NEEDED FOR BUF2 ARRAY IN BXABAB
C
      NABAB2=0
      NABAB3=0
      DO 33 IRREP=2,NIRREP
       ITMP=0
       DO 34 IRREPA=1,NIRREP
        IRREPB=DIRPRD(IRREPA,IRREP)
        ITMP=(VRT(IRREPA,1)*POP(IRREPB,1)+VRT(IRREPB,1)
     &      *POP(IRREPA,1))
     &      *VRT(IRREPA,1)*POP(IRREPB,1)+ITMP
34     CONTINUE
       NABAB2=MAX(NABAB2,ITMP)
33    CONTINUE
      ELSE
       NABAB1=0
       DO 35 IRREP=1,NIRREP
        ITMP1=IRPDPD(IRREP,ISYTYP(1,18))*IRPDPD(IRREP,ISYTYP(2,18))
        ITMP2=IRPDPD(IRREP,ISYTYP(1,19))*IRPDPD(IRREP,ISYTYP(2,19))
        ITMP3=IRPDPD(IRREP,ISYTYP(1,20))*IRPDPD(IRREP,ISYTYP(2,20))
        NABAB1=MAX(NABAB1,ITMP1,ITMP2,ITMP3)
35     CONTINUE
       NABAB11=NABAB1
      DO 36 IRREP1=2,NIRREP
       DO 37 IRREP2=1,NIRREP-1
        ITMP1=AOPOP(IRREP1)*AOPOP(IRREP2)
     &      *(VRT(IRREP1,1)*POP(IRREP2,1)
     &      +VRT(IRREP2,1)*POP(IRREP1,1))
        ITMP2=AOPOP(IRREP1)*AOPOP(IRREP2)
     &      *(VRT(IRREP1,2)*POP(IRREP2,2)
     &      +VRT(IRREP2,2)*POP(IRREP1,2))
        NABAB1=MAX(NABAB1,ITMP1,ITMP2)
        ITMP1=MAX(AOPOP(IRREP1),AOPOP(IRREP2))
     &      *(VRT(IRREP1,1)*POP(IRREP2,1)
     &      +VRT(IRREP2,1)*POP(IRREP1,1))
        ITMP2=MAX(AOPOP(IRREP1),AOPOP(IRREP2))
     &      *(VRT(IRREP1,2)*POP(IRREP2,2)
     &      +VRT(IRREP2,2)*POP(IRREP1,2))
        NABAB11=MAX(NABAB11,ITMP1,ITMP2)
37     CONTINUE
36    CONTINUE
C
C COMPUTE MAXIMUM MEMORY NEEDED FOR BUF2 ARRAY IN BXABAB
C
      NABAB2=0
      NABAB3=0
      DO 38 IRREP=2,NIRREP
       ITMP1=0
       ITMP2=0
       ITMP3=0
       ITMP4=0
       DO 39 IRREPA=1,NIRREP
        IRREPB=DIRPRD(IRREPA,IRREP)
        ITMP1=(VRT(IRREPA,1)*POP(IRREPB,1)+VRT(IRREPB,1)
     &      *POP(IRREPA,1))
     &      *VRT(IRREPA,2)*POP(IRREPB,2)+ITMP1
        ITMP2=(VRT(IRREPA,2)*POP(IRREPB,2)+VRT(IRREPB,2)
     &      *POP(IRREPA,2))
     &      *VRT(IRREPA,1)*POP(IRREPB,1)+ITMP2
        ITMP3=(VRT(IRREPA,1)*POP(IRREPB,1)+VRT(IRREPB,1)
     &      *POP(IRREPA,1))
     &      *VRT(IRREPA,1)*POP(IRREPB,1)+ITMP3
        ITMP4=(VRT(IRREPA,2)*POP(IRREPB,2)+VRT(IRREPB,2)
     &      *POP(IRREPA,2))
     &      *VRT(IRREPA,2)*POP(IRREPB,2)+ITMP4
39     CONTINUE
       NABAB2=MAX(NABAB2,ITMP1,ITMP2)
       NABAB3=MAX(NABAB3,ITMP3,ITMP4)
38    CONTINUE
      ENDIF
C*******************************************************************
C
C COMPUTE MAXIMUM MEMORY REQUIRED FOR BUF2 ARRAY IN BXABCD
C
       IF(IUHF.EQ.0) THEN
       NABCD3=0
       NABCD1=0
       NABCD11=0
       NABCD2=0
       ITMP=0
       DO 100 IRREPDO=2,NIRREP
        ITMP=0
        DO 110 IRREPI=1,NIRREP
         IRREPXIA=DIRPRD(IRREPI,IRREPDO)
         IBOT=MAX(IRREPI,IRREPXIA)+1
         CALL IZERO(IDID,NIRREP)
         DO 111 IRREPTMP=IBOT,NIRREP
          IRREPXIC=DIRPRD(IRREPTMP,IRREPDO)
          IF(IRREPI.NE.IRREPTMP.AND.IRREPI.NE.IRREPXIC.AND.
     &      IDID(IRREPTMP).EQ.0)THEN
          IRREPXIB=MIN(IRREPTMP,IRREPXIC)
          IRREPXIC=MAX(IRREPTMP,IRREPXIC)
          IDID(IRREPXIB)=1
          IDID(IRREPXIC)=1
          NSIZVR=VRT(IRREPXIA,1)
          NSIZOR=POP(IRREPI,1)
          NSIZVL1=VRT(IRREPXIC,1)
          NSIZOL1=POP(IRREPXIB,1)
          NSIZVL2=VRT(IRREPXIB,1)
          NSIZOL2=POP(IRREPXIC,1)
          ITMP=(NSIZVL1*NSIZOL1+NSIZVL2*NSIZOL2)*NSIZVR*NSIZOR+ITMP
         ENDIF
111      CONTINUE
         NABCD2=MAX(NABCD2,ITMP)
110     CONTINUE
C
C COMPUTE MAXIMUM MEMORY NEEDED FOR BUF1 ARRAY IN BXABCD
C
        DO 120 IRREPI=1,NIRREP
         IRREPXIA=DIRPRD(IRREPI,IRREPDO)
         IF(IRREPXIA.LT.IRREPI)GOTO 117
         IBOT=MAX(IRREPI,IRREPXIA)+1
         CALL IZERO(IDID,8)
         DO 130 ITMP=IBOT,NIRREP
          IRREPXIC=DIRPRD(ITMP,IRREPDO)
          ITMP2=MIN(IRREPXIC,ITMP)
          IRREPXIC=MAX(ITMP,IRREPXIC)
          IRREPXIB=ITMP2
          IF(MAX(IDID(IRREPXIC),IDID(IRREPXIB)).NE.0)GOTO 130
          IDID(IRREPXIC)=1
          IDID(IRREPXIB)=1
          ITMPX=AOPOP(IRREPXIC)*AOPOP(IRREPXIB)*(VRT(IRREPXIA,1)
     &          *POP(IRREPI,1)+POP(IRREPXIA,1)*VRT(IRREPI,1))
          NABCD1=MAX(NABCD1,ITMPX)
          NABCD11=MAX(NABCD11,ITMPX)
130      CONTINUE
117      ITMPX=IRPDPD(IRREPI,ISYTYP(1,18))*IRPDPD(IRREPI,ISYTYP(2,18))
         NABCD1=MAX(NABCD1,ITMPX) 
         NABCD11=MAX(NABCD11,ITMPX) 
120     CONTINUE
100    CONTINUE
C
       ELSE
C
C COMPUTE MAXIMUM MEMORY NEEDED FOR BUF1 ARRAY IN BXABCD
C
       NABCD3=0
       NABCD1=0
       NABCD11=0
       NABCD2=0
       DO 200 IRREPDO=2,NIRREP
        ITMP1=0
        ITMP2=0
        ITMP3=0
        ITMP4=0
        DO 210 IRREPI=1,NIRREP
         IRREPXIA=DIRPRD(IRREPI,IRREPDO)
         IBOT=MAX(IRREPI,IRREPXIA)+1
         CALL IZERO(IDID,NIRREP)
         DO 211 IRREPTMP=IBOT,NIRREP
          IRREPXIC=DIRPRD(IRREPTMP,IRREPDO)
          IF(IRREPI.NE.IRREPTMP.AND.IRREPI.NE.IRREPXIC.AND.
     &      IDID(IRREPTMP).EQ.0)THEN
          IRREPXIB=MIN(IRREPTMP,IRREPXIC)
          IRREPXIC=MAX(IRREPTMP,IRREPXIC)
          IDID(IRREPXIB)=1
          IDID(IRREPXIC)=1
          NSIZVR=VRT(IRREPXIA,2)
          NSIZOR=POP(IRREPI,2)
          NSIZVL1=VRT(IRREPXIC,1)
          NSIZOL1=POP(IRREPXIB,1)
          NSIZVL2=VRT(IRREPXIB,1)
          NSIZOL2=POP(IRREPXIC,1)
          ITMP1=(NSIZVL1*NSIZOL1+NSIZVL2*NSIZOL2)*NSIZVR*NSIZOR+ITMP1
          NSIZVR=VRT(IRREPXIA,1)
          NSIZOR=POP(IRREPI,1)
          ITMP3=(NSIZVL1*NSIZOL1+NSIZVL2*NSIZOL2)*NSIZVR*NSIZOR+ITMP3
          NSIZVL1=VRT(IRREPXIC,2)
          NSIZOL1=POP(IRREPXIB,2)
          NSIZVL2=VRT(IRREPXIB,2)
          NSIZOL2=POP(IRREPXIC,2)
          ITMP2=(NSIZVL1*NSIZOL1+NSIZVL2*NSIZOL2)*NSIZVR*NSIZOR+ITMP2
          NSIZOR=POP(IRREPI,2)
          NSIZVR=VRT(IRREPXIA,2)
          ITMP4=(NSIZVL1*NSIZOL1+NSIZVL2*NSIZOL2)*NSIZVR*NSIZOR+ITMP4
         ENDIF
211      CONTINUE
         NABCD2=MAX(NABCD2,ITMP1,ITMP2)
         NABCD3=MAX(NABCD3,ITMP3,ITMP4)
210     CONTINUE
        DO 220 IRREPI=1,NIRREP
         IRREPXIA=DIRPRD(IRREPI,IRREPDO)
         IF(IRREPXIA.LT.IRREPI)GOTO 217
         IBOT=MAX(IRREPI,IRREPXIA)+1
         CALL IZERO(IDID,8)
         DO 230 ITMP=IBOT,NIRREP
          IRREPXIC=DIRPRD(ITMP,IRREPDO)
          ITMP2=MIN(IRREPXIC,ITMP)
          IRREPXIC=MAX(ITMP,IRREPXIC)
          IRREPXIB=ITMP2
          IF(MAX(IDID(IRREPXIC),IDID(IRREPXIB)).NE.0)GOTO 230
          IDID(IRREPXIC)=1
          IDID(IRREPXIB)=1
          ITMPX1=AOPOP(IRREPXIC)*AOPOP(IRREPXIB)*(VRT(IRREPXIA,1)
     &          *POP(IRREPI,1)+POP(IRREPXIA,1)*VRT(IRREPI,1))
          ITMPX2=AOPOP(IRREPXIC)*AOPOP(IRREPXIB)*(VRT(IRREPXIA,2)
     &          *POP(IRREPI,2)+POP(IRREPXIA,2)*VRT(IRREPI,2))
          NABCD1=MAX(NABCD1,ITMPX1,ITMPX2)
          ITMPX1=MAX(AOPOP(IRREPXIC),AOPOP(IRREPXIB))*
     &           (VRT(IRREPXIA,1)*POP(IRREPI,1)+POP(IRREPXIA,1)*
     &            VRT(IRREPI,1))
          ITMPX2=MAX(AOPOP(IRREPXIC),AOPOP(IRREPXIB))*
     &           (VRT(IRREPXIA,2)*POP(IRREPI,2)+POP(IRREPXIA,2)*
     &            VRT(IRREPI,2))
          NABCD11=MAX(NABCD11,ITMPX1,ITMPX2)
230      CONTINUE
217      ITMPX1=IRPDPD(IRREPI,ISYTYP(1,18))*IRPDPD(IRREPI,ISYTYP(2,18))
         ITMPX2=IRPDPD(IRREPI,ISYTYP(1,19))*IRPDPD(IRREPI,ISYTYP(2,19))
         ITMPX3=IRPDPD(IRREPI,ISYTYP(1,20))*IRPDPD(IRREPI,ISYTYP(2,20))
         NABCD1=MAX(NABCD1,ITMPX1,ITMPX2,ITMPX3)
         NABCD11=MAX(NABCD11,ITMPX1,ITMPX2,ITMPX3)
220     CONTINUE
200    CONTINUE 
       ENDIF
C
C COMPUTE MEMORY REQUIRED FOR BUFTMP ARRAY 
C
      NBFTMP=0
      DO 300 IRREP=1,NIRREP
       ITMP=4*AOPOP(IRREP)*AOPOP(IRREP)
       NBFTMP=MAX(ITMP,NBFTMP)
300   CONTINUE
C
C  MEMORY REQUIREMENTS FOR ALL OTHER METHODS
C
      ELSE
C
C  ARRAY BUF IN BXAAAA2
C
      NAAAA1=0
      DO 400 IRREP=1,NIRREP
       ITMP1=AOPOP(IRREP)*(AOPOP(IRREP)+1)/2
       ITMP2=ITMP1*ITMP1
       NAAAA1=MAX(NAAAA1,ITMP2)
400   CONTINUE
      NAAAA2=0
      IF(IUHF.NE.0) THEN
       DO 410 IRREP=1,NIRREP
        NAAAA2=MAX(NAAAA2,AOPOP(IRREP)*(AOPOP(IRREP)+1)/2)
410    CONTINUE
      ENDIF
      NAAAA3=0
      DO 420 IRREP=1,NIRREP
       NAAAA3=MAX(NAAAA3,3*AOPOP(IRREP)*AOPOP(IRREP))
420   CONTINUE
C
C BXAABB2
C
      NAABB1=0
      NAABB2=0
      NAABB3=0
      DO 500 IRREP1=2,NIRREP
       DO 510 IRREP2=1,IRREP1-1
        ITMP1=NNP1O2(AOPOP(IRREP1))
        ITMP2=NNP1O2(AOPOP(IRREP2))
        NAABB1=MAX(NAABB1,ITMP1*ITMP2) 
510    CONTINUE
500   CONTINUE
      IF(IUHF.NE.0)NAABB2=NAABB1
      NAABB3=0
      DO 520 IRREP=1,NIRREP
       NAABB3=MAX(NAABB3,3*AOPOP(IRREP)*AOPOP(IRREP))
520   CONTINUE
C
C BXABAB2
C
      NABAB1=0
      NABAB2=0
      NABAB3=0
      DO 600 IRREP1=2,NIRREP
       DO 610 IRREP2=1,IRREP1-1
        ITMP=AOPOP(IRREP1)*AOPOP(IRREP2)
        NABAB1=MAX(ITMP*ITMP,NABAB1)
610    CONTINUE
600   CONTINUE
      IF(IUHF.NE.0)NABAB2=NABAB1
      NABAB3=0
      DO 620 IRREP=1,NIRREP
       NABAB3=MAX(NABAB3,3*AOPOP(IRREP)*AOPOP(IRREP))
620   CONTINUE
C
C BXABCD2
C
      NABCD1=0
      NABCD2=0
      NABCD3=0
      DO 700 IRREPDO=2,NIRREP
      DO 710 IRREPI=1,NIRREP
       IRREPXIA=DIRPRD(IRREPI,IRREPDO)
       IBOT=MAX(IRREPI,IRREPXIA)+1
       CALL IZERO(IDID,NIRREP)
       DO 711 IRREPTMP=IBOT,NIRREP
        IRREPXIC=DIRPRD(IRREPTMP,IRREPDO)
        IF(IRREPI.NE.IRREPTMP.AND.IRREPI.NE.IRREPXIC.AND.
     &    IDID(IRREPTMP).EQ.0)THEN
        IRREPXIB=MIN(IRREPTMP,IRREPXIC)
        IRREPXIC=MAX(IRREPTMP,IRREPXIC)
        IDID(IRREPXIB)=1
        IDID(IRREPXIC)=1
        ITMP=AOPOP(IRREPXIC)*AOPOP(IRREPXIB)*AOPOP(IRREPXIA)
     &       *AOPOP(IRREPI)
        NABCD1=MAX(NABCD1,ITMP)
        ITMP1=AOPOP(IRREPXIC)*AOPOP(IRREPXIB)
        ITMP2=AOPOP(IRREPXIA)*AOPOP(IRREPI)
        NABCD3=MAX(NABCD3,3*ITMP1,3*ITMP2)
        ENDIF
711    CONTINUE
710   CONTINUE
      IF(IUHF.NE.0)NABCD2=NABCD1
700   CONTINUE
C
      ENDIF
C
C ALL DONE
C
      RETURN
      END
