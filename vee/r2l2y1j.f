      SUBROUTINE R2L2Y1J(Y1,ICORE,MAXCOR,IUHF)
C
C Y2(ai) = -1/4 R(ef,mn)*L(ga,mn)*W(ef,gi)
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      DIMENSION ICORE(MAXCOR),Y1(*)
      COMMON/STATSYM/IRREPX
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA ONE,ONEM,ZILCH/1.0D0,-1.0D0,0.0D0/
C
      IF(IUHF.EQ.0)THEN
C
C SPIN ADAPTED RHF CODE.
C
C          
C Y2(AI)= - [2R(Mn,Ef)-R(Mn,Fe)]*L(Mn,Ga)*W(Gi,Ef)
C
C
       DO 10 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPX)
        LISTL=446
        LISTR=463
        LISTW=30
        DISSYR=IRPDPD(IRREPL,ISYTYP(1,LISTR))
        NUMDSR=IRPDPD(IRREPR,ISYTYP(2,LISTR))
        DISSYW=IRPDPD(IRREPL,ISYTYP(1,LISTW))
        NUMDSW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
        DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
        NUMDSL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
        DISSYQ=NUMDSR
        NUMDSQ=NUMDSW
        MAXT=MAX(NUMDSR,DISSYR,DISSYW,NUMDSW,NUMDSQ,DISSYQ,
     &           NUMDSL,DISSYL)
        I000=1
        I010=I000+IINTFP*MAX(NUMDSR*DISSYR,NUMDSL*DISSYL)
        I020=I010+IINTFP*NUMDSR*NUMDSW
        ITMP1=I020
        ITMP2=ITMP1+IINTFP*MAXT
        ITMP3=ITMP2+IINTFP*MAXT 
        CALL GETLST(ICORE(I000),1,NUMDSR,1,IRREPR,LISTR)
        CALL SPINAD3(IRREPL,VRT(1,1),DISSYR,NUMDSR,ICORE(I000),
     &               ICORE(ITMP1),ICORE(ITMP2))
C
        MAXDIS=(MAXCOR-I020)/(IINTFP*MAX(1,DISSYW)) 
        IF(MAXDIS.EQ.0) CALL INSMEM('R2L2Y1J',MAXDIS*DISSYW,DISSYW)
C
        NLEFT=NUMDSW
        IOFFLIST=1
C 
1       CONTINUE
C 
         NREAD=MIN(NLEFT,MAXDIS)
         NLEFT=NLEFT-NREAD
    
         CALL GETLST(ICORE(I020),IOFFLIST,NREAD,1,IRREPL,LISTW)
C                              _       +
C FORM INTERMEDIATE Q(Mn,Gi) = R(Ef,Mn) * W(Ef,Gi)
C
         IOFFR=I000
         IOFFW=I020
         IOFFQ=I010+IINTFP*NUMDSR*(IOFFLIST-1)
         CALL XGEMM('T','N',NUMDSR,NREAD,DISSYR,ONE,ICORE(IOFFR),
     &             DISSYR,ICORE(IOFFW),DISSYW,ZILCH,ICORE(IOFFQ),NUMDSR)
C
         IOFFLIST=IOFFLIST+NREAD
C
        IF(NLEFT.NE.0) GO TO 1
C
C NOW READ IN L(Ga,Mn) AND TRANSPOSE TO L(Mn,Ga)
C
        CALL GETTRN(ICORE(I000),ICORE(ITMP1),DISSYL,NUMDSL,
     &              1,IRREPR,LISTL)
C
C FORM PRODUCT Y(ai) = L(Mn,Ga) * Q(Mn,Gi)
C
        IOFFL=I000
        IOFFQ=I010
        IOFFY=1 
        DO 11 IRREPI=1,NIRREP
         IRREPA=IRREPI
         IRREPG=DIRPRD(IRREPA,IRREPL)
         NUMI=POP(IRREPI,1)
         NUMA=VRT(IRREPA,1)
         NUMG=VRT(IRREPG,1)
         NROW=NUMA
         NCOL=NUMI
         NSUM=NUMDSL*NUMG
         CALL XGEMM('T','N',NROW,NCOL,NSUM,ONEM,ICORE(IOFFL),NSUM,
     &              ICORE(IOFFQ),NSUM,ONE,Y1(IOFFY),NROW)
         IOFFL=IOFFL+IINTFP*NROW*NSUM
         IOFFQ=IOFFQ+IINTFP*NCOL*NSUM
         IOFFY=IOFFY+IINTFP*NROW*NCOL
11      CONTINUE
C
10     CONTINUE
C
      ELSE
C
C Y2(AI) <= -1/4 R(EF,MN)*W(EF,GI)*L(GA,MN) 
C
       IOFFY0=1 
       DO 20 ISPIN=1,2
        DO 21 IRREPR=1,NIRREP
         IRREPL=DIRPRD(IRREPR,IRREPX)
         LISTL=443+ISPIN
         LISTR=460+ISPIN
         LISTW=26+ISPIN
         DISSYR=IRPDPD(IRREPL,ISYTYP(1,LISTR))
         NUMDSR=IRPDPD(IRREPR,ISYTYP(2,LISTR))
         DISSYW=IRPDPD(IRREPL,ISYTYP(1,LISTW))
         NUMDSW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
         DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
         NUMDSL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
         DISSYLX=IRPDPD(IRREPL,18+ISPIN)
         DISSYQ=NUMDSR
         NUMDSQ=NUMDSW
         I000=1
         I010=I000+IINTFP*MAX(NUMDSR*DISSYR,NUMDSL*DISSYLX)
         I020=I010+IINTFP*NUMDSR*NUMDSW
         ITMP1=I020
         CALL GETLST(ICORE(I000),1,NUMDSR,1,IRREPR,LISTR)
C
         MAXDIS=(MAXCOR-I020)/(IINTFP*MAX(1,DISSYW))
         IF(MAXDIS.EQ.0) CALL INSMEM('R2L2Y1J',MAXDIS*DISSYW,DISSYW)
C
         NLEFT=NUMDSW
         IOFFLIST=1
C
2        CONTINUE
C
          NREAD=MIN(NLEFT,MAXDIS)
          NLEFT=NLEFT-NREAD
C
          CALL GETLST(ICORE(I020),IOFFLIST,NREAD,1,IRREPL,LISTW)
C                                         +
C FORM INTERMEDIATE Q(M<N,GI) = R(E<F,M<N) * W(E<F,GI)
C
          IOFFR=I000
          IOFFW=I020
          IOFFQ=I010+IINTFP*NUMDSR*(IOFFLIST-1)
          CALL XGEMM('T','N',NUMDSR,NREAD,DISSYR,ONE,ICORE(IOFFR),
     &               DISSYR,ICORE(IOFFW),DISSYW,ZILCH,ICORE(IOFFQ),
     &               NUMDSR)
C
          IOFFLIST=IOFFLIST+NREAD
C
         IF(NLEFT.NE.0) GO TO 2
C
C NOW READ IN L(G<A,M<N) AND TRANSPOSE TO L(M<N,G<A)
C
         CALL GETTRN(ICORE(I000),ICORE(ITMP1),DISSYL,NUMDSL,
     &               1,IRREPR,LISTL)
C
C EXPAND TO L(M<N,GA)
C
         CALL SYMEXP(IRREPL,VRT(1,ISPIN),NUMDSL,ICORE(I000))
C
C FORM PRODUCT Y(AI) <= L(M<N,GA) * Q(M<N,GI)
C
         IOFFL=I000
         IOFFQ=I010
         IOFFY=IOFFY0
         DO 22 IRREPI=1,NIRREP
          IRREPA=IRREPI
          IRREPG=DIRPRD(IRREPA,IRREPL)
          NUMI=POP(IRREPI,ISPIN)
          NUMA=VRT(IRREPA,ISPIN)
          NUMG=VRT(IRREPG,ISPIN)
          NROW=NUMA
          NCOL=NUMI
          NSUM=NUMDSL*NUMG
          CALL XGEMM('T','N',NROW,NCOL,NSUM,ONEM,ICORE(IOFFL),NSUM,
     &               ICORE(IOFFQ),NSUM,ONE,Y1(IOFFY),NROW)
          IOFFL=IOFFL+IINTFP*NROW*NSUM
          IOFFQ=IOFFQ+IINTFP*NCOL*NSUM
          IOFFY=IOFFY+IINTFP*NROW*NCOL
22       CONTINUE
C
C Y2(AI) <= - R(Ef,Mn)*W(Ef,Ig)*L(Ag,Mn) [ISPIN=1]
C Y2(AI) <= - R(Ef,Mn)*W(Ef,Gi)*L(Ga,Mn) [ISPIN=2]
C
         LISTL=446
         LISTR=463
         LISTW=28+ISPIN
         DISSYR=IRPDPD(IRREPL,ISYTYP(1,LISTR))
         NUMDSR=IRPDPD(IRREPR,ISYTYP(2,LISTR))
         DISSYW=IRPDPD(IRREPL,ISYTYP(1,LISTW))
         NUMDSW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
         DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
         NUMDSL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
         DISSYQ=NUMDSR
         NUMDSQ=NUMDSW
         MAXT=MAX(NUMDSR,DISSYR,DISSYW,NUMDSW,NUMDSQ,DISSYQ,
     &            NUMDSL,DISSYL)
         I000=1
         I010=I000+IINTFP*MAX(NUMDSR*DISSYR,NUMDSL*DISSYL)
         I020=I010+IINTFP*NUMDSR*NUMDSW
         ITMP1=I020
         ITMP2=ITMP1+IINTFP*MAXT
         ITMP3=ITMP2+IINTFP*MAXT 
         CALL GETLST(ICORE(I000),1,NUMDSR,1,IRREPR,LISTR)
C
         MAXDIS=(MAXCOR-I020)/(IINTFP*MAX(1,DISSYW)) 
         IF(MAXDIS.EQ.0) CALL INSMEM('R2L2Y1J',MAXDIS*DISSYW,DISSYW)
C
         NLEFT=NUMDSW
         IOFFLIST=1
C
3        CONTINUE
C
          NREAD=MIN(NLEFT,MAXDIS)
          NLEFT=NLEFT-NREAD
C
          CALL GETLST(ICORE(I020),IOFFLIST,NREAD,1,IRREPL,LISTW)
C                                       +
C FORM INTERMEDIATE Q(Mn,Ig)  = R(Ef,Mn) * W(Ef,Ig) [ISPIN=1] 
C                   Q(Mn,Gi)  = R(Ef,Mn) * W(Ef,Gi) [ISPIN=2] 
C
          IOFFR=I000
          IOFFW=I020
          IOFFQ=I010+IINTFP*NUMDSR*(IOFFLIST-1)
          CALL XGEMM('T','N',NUMDSR,NREAD,DISSYR,ONE,ICORE(IOFFR),
     &               DISSYR,ICORE(IOFFW),DISSYW,ZILCH,ICORE(IOFFQ),
     &               NUMDSR)
C
          IOFFLIST=IOFFLIST+NREAD
C
         IF(NLEFT.NE.0) GO TO 3
C
C NOW READ IN L(Ag,Mn) AND TRANSPOSE TO L(Mn,Ag) [ISPIN=1]
C NOW READ IN L(Ga,Mn) AND TRANSPOSE TO L(Mn,Ga) [ISPIN=2]
C
         CALL GETTRN(ICORE(I000),ICORE(ITMP1),DISSYL,NUMDSL,
     &               1,IRREPR,LISTL)
C
C Q(Mn,Ig) -> Q(Mn,gI) [ISPIN=1 ONLY]
C L(Mn,Ag) -> L(Mn,gA) [ISPIN=1 ONLY]
C
         IF(ISPIN.EQ.1)THEN
          CALL SYMTR1(IRREPL,VRT(1,1),VRT(1,2),NUMDSL,ICORE(I000),
     &                ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
          CALL SYMTR1(IRREPL,POP(1,1),VRT(1,2),DISSYQ,ICORE(I010),
     &                ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
         ENDIF
C
C FORM PRODUCT Y(AI) <= L(Mn,gA)*Q(Mn,gI) [ISPIN=1]
C FORM PRODUCT Y(AI) <= L(Mn,Ga)*Q(Mn,Gi) [ISPIN=2]
C
         IOFFL=I000
         IOFFQ=I010
         IOFFY=IOFFY0
         DO 23 IRREPI=1,NIRREP
          IRREPA=IRREPI
          IRREPG=DIRPRD(IRREPA,IRREPL)
          NUMI=POP(IRREPI,ISPIN)
          NUMA=VRT(IRREPA,ISPIN)
          NUMG=VRT(IRREPG,3-ISPIN)
          NROW=NUMA
          NCOL=NUMI
          NSUM=NUMDSL*NUMG
          CALL XGEMM('T','N',NROW,NCOL,NSUM,ONEM,ICORE(IOFFL),NSUM,
     &               ICORE(IOFFQ),NSUM,ONE,Y1(IOFFY),NROW)
          IOFFL=IOFFL+IINTFP*NROW*NSUM
          IOFFQ=IOFFQ+IINTFP*NCOL*NSUM
          IOFFY=IOFFY+IINTFP*NROW*NCOL
23       CONTINUE
C
21      CONTINUE
C
        IOFFY0=IOFFY0+IINTFP*NT(ISPIN)
C
20     CONTINUE
C
      ENDIF
C
      RETURN
      END