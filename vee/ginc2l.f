      SUBROUTINE GINC2L(IRREPX,ICORE,MAXCOR,IUHF,LISTZ0)
C
C SOLVES FOR THE THREE-BODY CONTRIBUTIONS TO THE LEFT-HAND
C C2 VECTOR.
C
C LIKE DFT2INT2 EXCEPT THAT THE ONE-PARTICLE TERMS (G) ARE NOT
C SYMMETRIC, WHILE THE T AMPLITUDES ARE TOTALLY SYMMETRIC.
C
C   Z(ab,ij) = W(ae,ij) * G(be) - W(ab,im) * G(mj)
C
C
C THE CONTRACTIONS ARE CARRIED OUT USING RING-ORDERED LISTS :
C
C   Z(ai,bj) = - W(ai,bm) * G(mj) + W(ai,ej) * G(be)
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,HALF
C
      DIMENSION ICORE(MAXCOR)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMLOC/ ISYMOFF(8,8,25)
C
      DATA ONE  /1.0/
      DATA ONEM /-1.0/
      DATA HALF /0.5/
      DATA ZILCH/0.0/
C
C FIRST DO THE ALPHA-BETA SPIN CASE
C
C
C LOOP OVER IRREPS
C
      DO 100 IRREPZR=1,NIRREP
       IRREPZL=DIRPRD(IRREPZR,IRREPX)
       LISTW=18
       LISTZ=LISTZ0+2
       DISSYZ=IRPDPD(IRREPZL,ISYTYP(1,LISTZ))
       NUMDSZ=IRPDPD(IRREPZR,ISYTYP(2,LISTZ))
       I000=1
       I010=I000+IINTFP*NUMDSZ*DISSYZ
       IRREPWR=IRREPZL
       IRREPWL=IRREPZL
       DISSYW=IRPDPD(IRREPWL,ISYTYP(1,LISTW))
       NUMDSW=IRPDPD(IRREPWR,ISYTYP(2,LISTW))
       MAXW=MAX(DISSYW,NUMDSW,DISSYZ,NUMDSZ)
       I020=I010+IINTFP*NUMDSW*DISSYW
       I030=I020+IINTFP*IRPDPD(IRREPX,22)
       I040=I030+IINTFP*IRPDPD(IRREPX,20)
       I050=I040+IINTFP*MAXW
       I060=I050+IINTFP*MAXW
       I070=I060+IINTFP*MAXW
       CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPWR,LISTW)
       CALL GETLST(ICORE(I020),1,1,1,1+IUHF,491)
       CALL GETLST(ICORE(I030),1,1,1,1+IUHF,492)
C
C CARRY OUT  -W(AI,bm) * G(mj) CONTRACTION
C
       IOFFG0=I020
       IOFFW0=I010
       IOFFZ0=I000
       DO IRREPM=1,NIRREP
          IRREPJ=DIRPRD(IRREPM,IRREPX)
          IRREPB=DIRPRD(IRREPM,IRREPWR)
          NUMM=POP(IRREPM,2)
          NUMJ=POP(IRREPJ,2)
          NUMB=VRT(IRREPB,2)
          NROW=DISSYW*NUMB
          NCOL=NUMJ
          NSUM=NUMM
          IOFFG=IOFFG0+IINTFP*(ISYMOFF(IRREPJ,IRREPX,22)-1)
          IOFFW=IOFFW0+IINTFP*DISSYW*(ISYMOFF(IRREPM,IRREPWR,10)-1)
          IOFFZ=IOFFZ0+IINTFP*DISSYZ*(ISYMOFF(IRREPJ,IRREPZR,10)-1)
cYAU - save dgemm
          IF (NSUM.EQ.0) THEN
             CALL ZERO(ICORE(IOFFZ),NROW*NCOL)
          ELSE
             CALL XGEMM('N','N',NROW,NCOL,NSUM,
     &                  ONEM, ICORE(IOFFW),NROW,
     &                        ICORE(IOFFG),NSUM,
     &                  ZILCH,ICORE(IOFFZ),NROW)
          END IF
       END DO
C
C NOW TRANSPOSE KET INDICES
C
       CALL SYMTR1(IRREPZR,VRT(1,2),POP(1,2),DISSYW,ICORE(I000),
     &             ICORE(I040),ICORE(I050),ICORE(I060))
       CALL SYMTR1(IRREPWR,VRT(1,2),POP(1,2),DISSYW,ICORE(I010),
     &             ICORE(I040),ICORE(I050),ICORE(I060))
C
C NOW CARRY OUT W(AI,je) * G(eb) CONTRACTION
C
       IOFFG0=I030
       IOFFW0=I010
       IOFFZ0=I000
       DO 120 IRREPE=1,NIRREP
        IRREPB=DIRPRD(IRREPE,IRREPX)
        IRREPJ=DIRPRD(IRREPE,IRREPWR)
        NUMJ=POP(IRREPJ,2)
        NUME=VRT(IRREPE,2)
        NUMB=VRT(IRREPB,2)
        NROW=DISSYW*NUMJ
        NCOL=NUMB
        NSUM=NUME
        IOFFG=IOFFG0+IINTFP*(ISYMOFF(IRREPB,IRREPX,20)-1)
        IOFFW=IOFFW0+IINTFP*DISSYW*(ISYMOFF(IRREPE,IRREPWR,17)-1)
        IOFFZ=IOFFZ0+IINTFP*DISSYZ*(ISYMOFF(IRREPB,IRREPZR,17)-1)
        CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFW),NROW,
     &             ICORE(IOFFG),NSUM,ONE,ICORE(IOFFZ),NROW)
120    CONTINUE
C
C FOR RHF, TRANSPOSE KET INDICES OF TARGET AND AUGMENT RING LISTS
C
       IF(IUHF.EQ.0)THEN
        CALL SYMTR1(IRREPZR,POP(1,2),VRT(1,2),DISSYW,ICORE(I000),
     &              ICORE(I040),ICORE(I050),ICORE(I060))
        CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREPZR,LISTZ)
        CALL SAXPY (NUMDSZ*DISSYZ,ONE,ICORE(I000),1,ICORE(I010),1)
        CALL PUTLST(ICORE(I010),1,NUMDSZ,1,IRREPZR,LISTZ)
       ELSE
C
C DO TRANSPOSITION:
C
C      Z(AM,jb) -> Z(jb,AM)
C
        CALL TRANSP(ICORE(I000),ICORE(I020),NUMDSZ,DISSYZ)
c YAU : old
c       CALL ICOPY(NUMDSZ*DISSYZ*IINTFP,ICORE(I020),1,ICORE(I000),1)
c YAU : new
        CALL DCOPY(NUMDSZ*DISSYZ,ICORE(I020),1,ICORE(I000),1)
c YAU : end
C
        IRREPWR=IRREPZR
        IRREPWL=IRREPZR
        LISTW=17
        DISSYW=IRPDPD(IRREPWL,ISYTYP(1,LISTW))
        NUMDSW=IRPDPD(IRREPWR,ISYTYP(2,LISTW))
        MAXW=MAX(DISSYW,NUMDSW,DISSYZ,NUMDSZ)
        I020=I010+IINTFP*NUMDSW*DISSYW
        I030=I020+IINTFP*IRPDPD(IRREPX,21)
        I040=I030+IINTFP*IRPDPD(IRREPX,19)
        I050=I040+IINTFP*MAXW
        I060=I050+IINTFP*MAXW
        I070=I060+IINTFP*MAXW
        CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPWR,LISTW)
C
C      W(bj,AM) -> W(jb,AM)
C
        CALL SYMTR3(IRREPWL,VRT(1,2),POP(1,2),DISSYW,NUMDSW,
     &              ICORE(I010),ICORE(I040),ICORE(I050),ICORE(I060))
C
        CALL GETLST(ICORE(I020),1,1,1,1,491)
        CALL GETLST(ICORE(I030),1,1,1,1,492)
C
C CARRY OUT  CONTRACTION :
C
C    -W(jb,AM) * G(MI) 
C
        IOFFG0=I020
        IOFFW0=I010
        IOFFZ0=I000
        DO 210 IRREPM=1,NIRREP
         IRREPI=DIRPRD(IRREPM,IRREPX)
         IRREPA=DIRPRD(IRREPM,IRREPWL)
         NUMM=POP(IRREPM,1)
         NUMI=POP(IRREPI,1)
         NUMA=VRT(IRREPA,1)
         NROW=DISSYW*NUMA
         NCOL=NUMI
         NSUM=NUMM
         IOFFG=IOFFG0+IINTFP*(ISYMOFF(IRREPI,IRREPX,21)-1)
         IOFFW=IOFFW0+IINTFP*DISSYW*(ISYMOFF(IRREPM,IRREPWR,9)-1)
         IOFFZ=IOFFZ0+IINTFP*NUMDSZ*(ISYMOFF(IRREPI,IRREPZL,9)-1)
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONEM,ICORE(IOFFW),NROW,
     &              ICORE(IOFFG),NSUM,ONE,ICORE(IOFFZ),NROW)
210     CONTINUE
C
C NOW TRANSPOSE KET INDICES
C
C      Z(jb,EI) -> Z(jb,IE)
C      W(jb,EI) -> W(je,IE)
C
        CALL SYMTR1(IRREPZL,VRT(1,1),POP(1,1),NUMDSZ,ICORE(I000),
     &              ICORE(I040),ICORE(I050),ICORE(I060))
        CALL SYMTR1(IRREPWL,VRT(1,1),POP(1,1),DISSYW,ICORE(I010),
     &              ICORE(I040),ICORE(I050),ICORE(I060))
C
C
C NOW CARRY OUT CONTRACTION :
C
C        W(jb,IE) * G(EA)
C
        IOFFF0=I030
        IOFFW0=I010
        IOFFZ0=I000
        DO 220 IRREPE=1,NIRREP
         IRREPA=DIRPRD(IRREPE,IRREPX)
         IRREPI=DIRPRD(IRREPE,IRREPWL)
         NUMI=POP(IRREPI,1)
         NUME=VRT(IRREPE,1)
         NUMA=VRT(IRREPA,1)
         NROW=DISSYW*NUMI
         NCOL=NUMA
         NSUM=NUME
         IOFFF=IOFFF0+IINTFP*(ISYMOFF(IRREPA,IRREPX,19)-1)
         IOFFW=IOFFW0+IINTFP*DISSYW*(ISYMOFF(IRREPE,IRREPWL,16)-1)
         IOFFZ=IOFFZ0+IINTFP*NUMDSZ*(ISYMOFF(IRREPA,IRREPZL,16)-1)
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFW),NROW,
     &              ICORE(IOFFF),NSUM,ONE,ICORE(IOFFZ),NROW)
220     CONTINUE
C
C WE NOW HAVE WHAT WE WANT, BUT IN THE RATHER UNFORTUNATE ORDERING
C   Z(jb,IA).  LET US NOW REARRANGE THIS (TEDIOUSLY) TO Z(AI,bj)
C
        CALL SYMTR1(IRREPZL,POP(1,1),VRT(1,1),NUMDSZ,ICORE(I000),
     &              ICORE(I040),ICORE(I050),ICORE(I060))
        CALL TRANSP(ICORE(I000),ICORE(I010),DISSYZ,NUMDSZ)
c YAU : old
c       CALL ICOPY(NUMDSZ*DISSYZ*IINTFP,ICORE(I010),1,ICORE(I000),1)
c YAU : new
        CALL DCOPY(NUMDSZ*DISSYZ,ICORE(I010),1,ICORE(I000),1)
c YAU : end
        CALL SYMTR1(IRREPZR,POP(1,2),VRT(1,2),DISSYZ,ICORE(I000),
     &              ICORE(I040),ICORE(I050),ICORE(I060))
        CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREPZR,LISTZ)
        CALL SAXPY (NUMDSZ*DISSYZ,ONE,ICORE(I010),1,ICORE(I000),1)
        CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPZR,LISTZ)
C
       ENDIF
C
100   CONTINUE
C
      IF(IUHF.EQ.0)RETURN
C
C DO AAAA AND BBBB SPIN CASES
C
      DO 500 ISPIN=1,2
       DO 510 IRREPZR=1,NIRREP
        IRREPZL=DIRPRD(IRREPZR,IRREPX)
        IRREPTR=IRREPZL
        IRREPTL=IRREPZL
        LISTT=18+ISPIN
        LISTZ=LISTZ0-1+ISPIN
        DISSYZ=IRPDPD(IRREPZL,ISYTYP(1,LISTZ))
        NUMDSZ=IRPDPD(IRREPZR,ISYTYP(2,LISTZ))
        DISSYT=IRPDPD(IRREPTL,ISYTYP(1,LISTT))
        NUMDST=IRPDPD(IRREPTR,ISYTYP(2,LISTT))
        MAXT=MAX(DISSYT,NUMDST,DISSYZ,NUMDSZ)
        I000=1
        I010=I000+IINTFP*NUMDSZ*DISSYZ
        I020=I010+IINTFP*NUMDST*DISSYT
        I030=I020+IINTFP*IRPDPD(IRREPX,20+ISPIN)
        I040=I030+IINTFP*IRPDPD(IRREPX,18+ISPIN)
        I050=I040+IINTFP*MAXT
        I060=I050+IINTFP*MAXT
        I070=I060+IINTFP*MAXT
        CALL GETLST(ICORE(I010),1,NUMDST,1,IRREPTR,LISTT)
        CALL GETLST(ICORE(I020),1,1,1,ISPIN,491)
        CALL GETLST(ICORE(I030),1,1,1,ISPIN,492)
        CALL ZERO  (ICORE(I000),DISSYZ*NUMDSZ)
C
C CARRY OUT  CONTRACTION :
C
C    -T(AI,BM) * F(MJ) 
C
        IOFFF0=I020
        IOFFT0=I010
        IOFFZ0=I000
        DO IRREPM=1,NIRREP
           IRREPJ=DIRPRD(IRREPM,IRREPX)
           IRREPB=DIRPRD(IRREPM,IRREPTR)
           NUMM=POP(IRREPM,ISPIN)
           NUMJ=POP(IRREPJ,ISPIN)
           NUMB=VRT(IRREPB,ISPIN)
           NROW=DISSYT*NUMB
           NCOL=NUMJ
           NSUM=NUMM
           IOFFF=IOFFF0
     &          +IINTFP*(ISYMOFF(IRREPJ,IRREPX,20+ISPIN)-1)
           IOFFT=IOFFT0
     &          +IINTFP*DISSYT*(ISYMOFF(IRREPM,IRREPTR,8+ISPIN)-1)
           IOFFZ=IOFFZ0
     &          +IINTFP*DISSYZ*(ISYMOFF(IRREPJ,IRREPZR,8+ISPIN)-1)
cYAU - save dgemm
           IF (NSUM.EQ.0) THEN
              CALL ZERO(ICORE(IOFFZ),NROW*NCOL)
           ELSE
              CALL XGEMM('N','N',NROW,NCOL,NSUM,
     &                   ONEM, ICORE(IOFFT),NROW,
     &                         ICORE(IOFFF),NSUM,
     &                   ZILCH,ICORE(IOFFZ),NROW)
           END IF
        END DO
C
C NOW TRANSPOSE KET INDICES
C
        CALL SYMTR1(IRREPZR,VRT(1,ISPIN),POP(1,ISPIN),DISSYZ,
     &              ICORE(I000),ICORE(I040),ICORE(I050),ICORE(I060))
        CALL SYMTR1(IRREPTR,VRT(1,ISPIN),POP(1,ISPIN),DISSYT,
     &              ICORE(I010),ICORE(I040),ICORE(I050),ICORE(I060))
C
C NOW CARRY OUT CONTRACTION :
C
C        T(AI,JE) * F(EB)
C
C
        IOFFF0=I030
        IOFFT0=I010
        IOFFZ0=I000
        DO 530 IRREPE=1,NIRREP
         IRREPB=DIRPRD(IRREPE,IRREPX)
         IRREPJ=DIRPRD(IRREPE,IRREPTR)
         NUMJ=POP(IRREPJ,ISPIN)
         NUME=VRT(IRREPE,ISPIN)
         NUMB=VRT(IRREPB,ISPIN)
         NROW=DISSYT*NUMJ
         NCOL=NUMB
         NSUM=NUME
         IOFFF=IOFFF0+IINTFP*(ISYMOFF(IRREPB,IRREPX,18+ISPIN)-1)
         IOFFT=IOFFT0+IINTFP*DISSYT*(ISYMOFF(IRREPE,IRREPTR,15+ISPIN)-1)
         IOFFZ=IOFFZ0+IINTFP*DISSYZ*(ISYMOFF(IRREPB,IRREPZR,15+ISPIN)-1)
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFT),NROW,
     &              ICORE(IOFFF),NSUM,ONE,ICORE(IOFFZ),NROW)
530     CONTINUE
C
C TRANSPOSE KET INDICES OF TARGET AND WRITE TO RING LISTS
C
        CALL SYMTR1(IRREPZR,POP(1,ISPIN),VRT(1,ISPIN),DISSYZ,
     &              ICORE(I000),ICORE(I040),ICORE(I050),ICORE(I060))
        CALL SSCAL (NUMDSZ*DISSYZ,-HALF,ICORE(I000),1)
        CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREPZR,LISTZ)
        CALL SAXPY (NUMDSZ*DISSYZ,ONE,ICORE(I010),1,ICORE(I000),1)
        CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPZR,LISTZ)
510    CONTINUE
500   CONTINUE
C
      RETURN
      END
