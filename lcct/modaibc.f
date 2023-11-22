
      SUBROUTINE MODAIBC(ICORE,MAXCOR,IUHF,FACT)
C
C FORMS THE AIBC HBAR ELEMENTS
C
C  HBAR(ai,bc) = <ai||bc> + T(am) <mi||bc>
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,FACT
      DIMENSION ICORE(MAXCOR)
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA ONE,ZILCH,ONEM/1.0D0,0.0D0,-1.0D0/
C
C DO SPIN CASE 4 FIRST
C
C  HBAR(Ai,Bc) = T(AM) * W(Bc,Mi)    as W(Bci,M) * T(AM)
C
      LISTW=16
      LISTZ=30
      LISTT=90
      I000=1
      I010=I000+IINTFP*NT(1)
      CALL GETLST(ICORE(I000),1,1,1,1,90)
      DO 10 IRREPDO=1,NIRREP
       DISSYW=IRPDPD(IRREPDO,ISYTYP(1,LISTW))
       NUMDSW=IRPDPD(IRREPDO,ISYTYP(2,LISTW))
       DISSYZ=IRPDPD(IRREPDO,ISYTYP(1,LISTZ))
       NUMDSZ=IRPDPD(IRREPDO,ISYTYP(2,LISTZ))
       MAXW=MAX(NUMDSW,DISSYW,NUMDSZ,DISSYZ)
       I020=I010+IINTFP*MAX(NUMDSW*DISSYW,3*DISSYZ,3*NUMDSZ)
       I030=I020+IINTFP*NUMDSZ*DISSYZ
       CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPDO,LISTW)
       CALL GETLST(ICORE(I020),1,NUMDSZ,1,IRREPDO,LISTZ)
       ITMP1=I030
       ITMP2=ITMP1+IINTFP*MAXW
       ITMP3=ITMP2+IINTFP*MAXW
       CALL SYMTR1(IRREPDO,POP(1,1),POP(1,2),DISSYW,ICORE(I010),
     &             ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
       CALL SYMTR1(IRREPDO,VRT(1,1),POP(1,2),DISSYW,ICORE(I020),
     &             ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
       IOFFT=I000
       IOFFW=I010
       IOFFZ=I020
       DO 20 IRREPM=1,NIRREP
        IRREPA=IRREPM
        IRREPI=DIRPRD(IRREPM,IRREPDO)
        NUMM=POP(IRREPM,1)
        NUMA=VRT(IRREPA,1)
        NUMI=POP(IRREPI,2)
        NROW=DISSYZ*NUMI
        NCOL=NUMA
        NSUM=NUMM
        CALL XGEMM('N','T',NROW,NCOL,NSUM,FACT*ONEM,ICORE(IOFFW),NROW,
     &             ICORE(IOFFT),NUMA,ONE,ICORE(IOFFZ),NROW)
        IOFFW=IOFFW+IINTFP*NROW*NSUM
        IOFFT=IOFFT+IINTFP*NUMM*NUMA
        IOFFZ=IOFFZ+IINTFP*NROW*NCOL
20     CONTINUE
       CALL SYMTR1(IRREPDO,POP(1,2),VRT(1,1),DISSYZ,ICORE(I020),
     &             ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
       CALL PUTLST(ICORE(I020),1,NUMDSZ,1,IRREPDO,LISTZ)
10    CONTINUE
C
C RETURN FOR RHF
C
      IF(IUHF.EQ.0)RETURN
C
C DO SPIN CASE 3
C
C  HBAR(aI,bC) = T(am) * W(bC,aI)    as W(CbI,a) * T(am)
C
      LISTW=16
      LISTZ=29
      LISTT=90
      I000=1
      I010=I000+IINTFP*NT(2)
      CALL GETLST(ICORE(I000),1,1,1,2,90)
      DO 210 IRREPDO=1,NIRREP
       DISSYW=IRPDPD(IRREPDO,ISYTYP(1,LISTW))
       NUMDSW=IRPDPD(IRREPDO,ISYTYP(2,LISTW))
       DISSYZ=IRPDPD(IRREPDO,ISYTYP(1,LISTZ))
       NUMDSZ=IRPDPD(IRREPDO,ISYTYP(2,LISTZ))
       I020=I010+IINTFP*NUMDSW*DISSYW
       I030=I020+IINTFP*NUMDSZ*DISSYZ
       CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPDO,LISTW)
       CALL GETLST(ICORE(I020),1,NUMDSZ,1,IRREPDO,LISTZ)
       IOFFT=I000
       IOFFW=I010
       IOFFZ=I020
       DO 220 IRREPM=1,NIRREP
        IRREPA=IRREPM
        IRREPI=DIRPRD(IRREPM,IRREPDO)
        NUMM=POP(IRREPM,2)
        NUMA=VRT(IRREPA,2)
        NUMI=POP(IRREPI,1)
        NROW=DISSYZ*NUMI
        NCOL=NUMA
        NSUM=NUMM
        CALL XGEMM('N','T',NROW,NCOL,NSUM,FACT*ONEM,ICORE(IOFFW),NROW,
     &             ICORE(IOFFT),NUMA,ONE,ICORE(IOFFZ),NROW)
        IOFFW=IOFFW+IINTFP*NROW*NSUM
        IOFFT=IOFFT+IINTFP*NUMM*NUMA
        IOFFZ=IOFFZ+IINTFP*NROW*NCOL
220    CONTINUE
       CALL PUTLST(ICORE(I020),1,NUMDSZ,1,IRREPDO,LISTZ)
210   CONTINUE
C
C DO SPIN CASES 1 AND 2
C
C  HBAR(ai,bc) = T(am) <mi||bc> as W(B<C I,M) * T(AM)
C
      DO 300 ISPIN=1,2
       LISTW=13+ISPIN
       LISTZ=26+ISPIN
       LISTT=90
       I000=1
       I010=I000+IINTFP*NT(ISPIN)
       CALL GETLST(ICORE(I000),1,1,1,ISPIN,90)
       DO 310 IRREPDO=1,NIRREP
        DISSYW=IRPDPD(IRREPDO,ISYTYP(1,LISTW))
        NUMDSW=IRPDPD(IRREPDO,ISYTYP(2,LISTW))
        NUMDSW0=IRPDPD(IRREPDO,20+ISPIN)
        DISSYZ=IRPDPD(IRREPDO,ISYTYP(1,LISTZ))
        NUMDSZ=IRPDPD(IRREPDO,ISYTYP(2,LISTZ))
        MAXZ=MAX(NUMDSZ,DISSYZ)
        I020=I010+IINTFP*NUMDSW0*DISSYW
        I030=I020+IINTFP*NUMDSZ*DISSYZ
        CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPDO,LISTW)
        CALL SYMEXP(IRREPDO,POP(1,ISPIN),DISSYW,ICORE(I010))
        CALL GETLST(ICORE(I020),1,NUMDSZ,1,IRREPDO,LISTZ)
        ITMP1=I030
        ITMP2=ITMP1+IINTFP*MAXZ
        ITMP3=ITMP2+IINTFP*MAXZ
        CALL SYMTR1(IRREPDO,VRT(1,ISPIN),POP(1,ISPIN),DISSYZ,
     &              ICORE(I020),ICORE(ITMP1),ICORE(ITMP2),
     &              ICORE(ITMP3))
        IOFFT=I000
        IOFFW=I010
        IOFFZ=I020
        DO 320 IRREPM=1,NIRREP
         IRREPA=IRREPM
         IRREPI=DIRPRD(IRREPM,IRREPDO)
         NUMM=POP(IRREPM,ISPIN)
         NUMA=VRT(IRREPA,ISPIN)
         NUMI=POP(IRREPI,ISPIN)
         NROW=DISSYZ*NUMI
         NCOL=NUMA
         NSUM=NUMM
         CALL XGEMM('N','T',NROW,NCOL,NSUM,FACT*ONE,ICORE(IOFFW),NROW,
     &              ICORE(IOFFT),NUMA,ONE,ICORE(IOFFZ),NROW)
         IOFFW=IOFFW+IINTFP*NROW*NSUM
         IOFFT=IOFFT+IINTFP*NUMM*NUMA
         IOFFZ=IOFFZ+IINTFP*NROW*NCOL
320     CONTINUE
        CALL SYMTR1(IRREPDO,POP(1,ISPIN),VRT(1,ISPIN),DISSYZ,
     &              ICORE(I020),ICORE(ITMP1),ICORE(ITMP2),
     &              ICORE(ITMP3))
        CALL PUTLST(ICORE(I020),1,NUMDSZ,1,IRREPDO,LISTZ)
310    CONTINUE
C
300   CONTINUE
C
      RETURN
      END
