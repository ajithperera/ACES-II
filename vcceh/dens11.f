      SUBROUTINE DENS11(C,DOO,DVV,IRREPC,ISPIN,IUHF)
C
C CALCULATES CONTRIBUTIONS TO THE TDA ONE-PARTICLE DENSITY MATRICES 
C
C                       +        
C    D(ij) = - SUM C(ai) C(aj)
C               a
C
C                           +
C    D(ab) = SUM C(am) C(bm)
C             m
C
C THE OCCUPIED-VIRTUAL DENSITY MATRIX ELEMENTS DO NOT INVOLVE 
C CONTRACTIONS AND ARE NOT CALCULATED HERE
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
      DIMENSION C(*),DOO(*),DVV(*)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
      COMMON/LISTDENS/LDENS
C
      DATA ONE,ONEM,ZILCH /1.0D0,-1.0D0,0.0D0/
C
      IRREPX=1
C
      IF(IUHF.EQ.0)THEN
       CALL SSCAL(NT(1),1.0D0/DSQRT(2.0D0),C,1)
      ENDIF
C 
C FORM OCCUPIED-OCCUPIED CONTRIBUTION
C
      DO 10 IRREPJ=1,NIRREP
       IRREPI=DIRPRD(IRREPJ,IRREPX)
       NUMJ=POP(IRREPJ,ISPIN)
       NUMI=POP(IRREPI,ISPIN)
       IRREPA=DIRPRD(IRREPI,IRREPC)
       NUMA=VRT(IRREPA,ISPIN)
       IOFFCL=ISYMOFF(IRREPI,IRREPC,8+ISPIN)
       IOFFCR=ISYMOFF(IRREPJ,IRREPC,8+ISPIN)
       IOFFD =ISYMOFF(IRREPJ,IRREPX,20+ISPIN)
       CALL XGEMM('T','N',NUMI,NUMJ,NUMA,ONEM,C(IOFFCL),NUMA,
     &             C(IOFFCR),NUMA,ZILCH,DOO(IOFFD),NUMI)
10    CONTINUE
      CALL PUTLST(DOO,1,1,1,ISPIN,ldens)
C
C FORM VIRTUAL-VIRTUAL CONTRIBUTION
C
      DO 20 IRREPB=1,NIRREP
       IRREPA=DIRPRD(IRREPB,IRREPX)
       NUMB=VRT(IRREPB,ISPIN)
       NUMA=VRT(IRREPA,ISPIN)
       IRREPM=DIRPRD(IRREPB,IRREPC)
       NUMM=POP(IRREPM,ISPIN)
       IOFFCL=ISYMOFF(IRREPM,IRREPC,8+ISPIN)
       IOFFCR=ISYMOFF(IRREPM,IRREPC,8+ISPIN)
       IOFFD =ISYMOFF(IRREPB,IRREPX,18+ISPIN)
       CALL XGEMM('N','T',NUMA,NUMB,NUMM,ONE,C(IOFFCR),NUMA,
     &            C(IOFFCL),NUMB,ZILCH,DVV(IOFFD),NUMA)
20    CONTINUE
      CALL PUTLST(DVV,1,1,1,2+ISPIN,ldens)
C
      CALL UPDMOI(1,NT(ISPIN),ISPIN,90,0,0)
      CALL UPDMOI(1,NT(ISPIN),ISPIN,190,0,0)
      CALL PUTLST(C,1,1,1,ISPIN,90)
      CALL PUTLST(C,1,1,1,ISPIN,190)
C
      RETURN
      END
