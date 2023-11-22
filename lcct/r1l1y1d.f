      SUBROUTINE R1L1Y1D(R1,ICORE,MAXCOR,IUHF)
C
C Y1(ab,ij) = P(ij)P(ab) [R(me)*L(fb,mj)*W(ea,fi)]
CC
C
C
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,TWO,HALF
      DIMENSION ICORE(MAXCOR),R1(*)
      COMMON/STATSYM/IRREPX
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMLOC/ISYMOFF(8,8,25)
C
      DATA ONE,ONEM,TWO,ZILCH/1.0D0,-1.0D0,2.0D0,0.0D0/
      DATA HALF/0.5D0/
C
      IF(IUHF.EQ.0)THEN
C
C SPIN ADAPTED RHF CODE.
C
C   Y(Ab,Ij) = Q(Ab,Ij) + Q(Ba,Ji)
C
C WHERE
C
C   Q(ab,ij) = 1/2 R(me) * [2L(Fb,Mj) - F(Bf,Mj)] * [2W(Fi,Ea)-W(Fi,Ae)]
C
C            - 1/2 R(me)*L(Mj,Bf)*W(If,Ea) - R(me)*L(Mj,Af)*W(If,Eb)
C
       LISTW=30
       LISTSAV=434
       I000=1
       DO 10 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPX)
        DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
        NUMDSW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
        DISSYX=IRPDPD(IRREPR,ISYTYP(1,LISTSAV))
        NUMDSX=IRPDPD(IRREPL,ISYTYP(2,LISTSAV))
        MAXW=MAX(DISSYW,NUMDSW,DISSYX,NUMDSX)
        MAXSIZ=MAX(DISSYW*NUMDSW,DISSYX*NUMDSX)
        I010=I000+IINTFP*MAXSIZ
        I020=I010+IINTFP*MAXSIZ
        ITMP1=I020
        ITMP2=ITMP1+IINTFP*MAXW
        ITMP3=ITMP2+IINTFP*MAXW
C                                            _
C READ W(Ea,Fi), SPIN-ADAPT AND TRANSPOSE TO W(Fi,aE)
C
        CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPR,LISTW)
        CALL SPINAD3(IRREPR,VRT(1,1),DISSYW,NUMDSW,ICORE(I010),
     &               ICORE(ITMP1),ICORE(ITMP2))
        CALL TRANSP(ICORE(I010),ICORE(I000),NUMDSW,DISSYW)
        CALL SYMTR1(IRREPR,VRT(1,1),VRT(1,1),NUMDSW,ICORE(I000),
     &              ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
C
C FORM PRODUCT X(Fi,aM) = W(Fi,aE) * R(EM) AND PUT ON SCRATCH LIST
C
        DO 20 IRREPM=1,NIRREP
         IRREPE=DIRPRD(IRREPM,IRREPX)
         IRREPA=DIRPRD(IRREPE,IRREPR)
         NUMM=POP(IRREPM,1)
         NUMA=VRT(IRREPA,1)
         NUME=VRT(IRREPE,1)
         NROW=DISSYX*NUMA
         NCOL=NUMM
         NSUM=NUME
         IOFFX=I010+IINTFP*(ISYMOFF(IRREPM,IRREPL,9)-1)*DISSYX
         IOFFW=I000+IINTFP*(ISYMOFF(IRREPE,IRREPR,19)-1)*NUMDSW
         IOFFR=1+IINTFP*(ISYMOFF(IRREPM,IRREPX,9)-1)
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFW),NROW,
     &              R1(IOFFR),NSUM,ZILCH,ICORE(IOFFX),NROW)
20      CONTINUE
C
        CALL PUTLST(ICORE(I010),1,NUMDSX,1,IRREPL,LISTSAV)
C
10     CONTINUE
C
C NOW RESORT X(Fi,aM) -> X(FM,ai)
C
       ISIZE=IDSYMSZ(IRREPX,ISYTYP(1,LISTSAV),ISYTYP(2,LISTSAV))
       I010=I000+IINTFP*ISIZE
       I020=I010+IINTFP*ISIZE
       CALL GETALL(ICORE(I000),ISIZE,IRREPX,LISTSAV)
       CALL SSTGEN(ICORE(I000),ICORE(I010),ISIZE,VRT(1,1),POP(1,1),
     &             VRT(1,1),POP(1,1),ICORE(I020),IRREPX,'1432')
       CALL PUTALL(ICORE(I010),ISIZE,IRREPX,LISTSAV)
C
C NOW FORM SPIN ADAPTED L AMPLITUDES AND STORE IN L(FM,bj)
C
       DO 30 IRREPR=1,NIRREP
        LISTL=437
        LISTZ=42
        IRREPL=DIRPRD(IRREPR,IRREPX)
        DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
        NUMDSL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
        DISSYX=IRPDPD(IRREPL,ISYTYP(1,LISTL))
        NUMDSX=IRPDPD(IRREPR,ISYTYP(2,LISTL))
        DISSYZ=IRPDPD(IRREPR,ISYTYP(1,LISTZ))
        NUMDSZ=IRPDPD(IRREPR,ISYTYP(2,LISTZ))
        I010=I000+IINTFP*NUMDSL*DISSYL
        I020=I010+IINTFP*MAX(NUMDSL*DISSYL,NUMDSX*DISSYX)
        CALL GETLST(ICORE(I000),1,NUMDSL,1,IRREPR,LISTL)
        CALL GETLST(ICORE(I010),1,NUMDSL,1,IRREPR,LISTL+2)
        CALL SSCAL (NUMDSL*DISSYL,TWO,ICORE(I000),1)
        CALL SAXPY (NUMDSL*DISSYL,ONEM,ICORE(I010),1,ICORE(I000),1)
C
        CALL GETLST(ICORE(I010),1,NUMDSX,1,IRREPR,LISTSAV)
C
C FORM PRODUCT Y(ai,bj) = X(FM,ai) * L(FM,bj) AND STORE ON LIST 42
C
        CALL XGEMM('T','N',NUMDSX,NUMDSL,DISSYX,HALF,ICORE(I010),DISSYX,
     &             ICORE(I000),DISSYL,ZILCH,ICORE(I020),DISSYZ)
        CALL PUTLST(ICORE(I020),1,NUMDSZ,1,IRREPR,LISTZ)
30     CONTINUE
C
C NOW DEAL WITH EXCHANGE CONTRIBUTION
C
C    R(ME) * L(Af,Mj) * W(Be,Fi)
C
       LISTW=30
       LISTSAV=434
       I000=1
       DO 110 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPX)
        DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
        NUMDSW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
        DISSYX=IRPDPD(IRREPR,ISYTYP(1,LISTSAV))
        NUMDSX=IRPDPD(IRREPL,ISYTYP(2,LISTSAV))
        MAXW=MAX(DISSYW,NUMDSW,DISSYX,NUMDSX)
        MAXSIZ=MAX(DISSYW*NUMDSW,DISSYX*NUMDSX)
        I010=I000+IINTFP*MAXSIZ
        I020=I010+IINTFP*MAXSIZ
        ITMP1=I020
        ITMP2=ITMP1+IINTFP*MAXW
        ITMP3=ITMP2+IINTFP*MAXW
C
C READ W(Be,Fi) AND TRANSPOSE TO W(Fi,Be)
C
        CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPR,LISTW)
        CALL TRANSP(ICORE(I010),ICORE(I000),NUMDSW,DISSYW)
C
C FORM PRODUCT X(Fi,Bm) = W(Fi,Be) * R(em) AND PUT ON SCRATCH LIST
C
        DO 120 IRREPM=1,NIRREP
         IRREPE=DIRPRD(IRREPM,IRREPX)
         IRREPB=DIRPRD(IRREPE,IRREPR)
         NUMM=POP(IRREPM,1)
         NUMB=VRT(IRREPB,1)
         NUME=VRT(IRREPE,1)
         NROW=DISSYX*NUMB
         NCOL=NUMM
         NSUM=NUME
         IOFFX=I010+IINTFP*(ISYMOFF(IRREPM,IRREPL,9)-1)*DISSYX
         IOFFW=I000+IINTFP*(ISYMOFF(IRREPE,IRREPR,19)-1)*NUMDSW
         IOFFR=1+IINTFP*(ISYMOFF(IRREPM,IRREPX,9)-1)
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFW),NROW,
     &              R1(IOFFR),NSUM,ZILCH,ICORE(IOFFX),NROW)
120     CONTINUE
C
        CALL PUTLST(ICORE(I010),1,NUMDSX,1,IRREPL,LISTSAV)
C
110    CONTINUE
C
C NOW RESORT X(Fi,Bm) -> X(Fm,Bi)
C
       ISIZE=IDSYMSZ(IRREPX,ISYTYP(1,LISTSAV),ISYTYP(2,LISTSAV))
       I010=I000+IINTFP*ISIZE
       I020=I010+IINTFP*ISIZE
       CALL GETALL(ICORE(I000),ISIZE,IRREPX,LISTSAV)
       CALL SSTGEN(ICORE(I000),ICORE(I010),ISIZE,VRT(1,1),POP(1,1),
     &             VRT(1,1),POP(1,1),ICORE(I020),IRREPX,'1432')
       CALL PUTALL(ICORE(I010),ISIZE,IRREPX,LISTSAV)
C
C NOW READ IN L(aJ,Fm)
C
       DO 130 IRREPR=1,NIRREP
        LISTL=439
        LISTZ=43
        IRREPL=DIRPRD(IRREPR,IRREPX)
        DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
        NUMDSL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
        DISSYX=IRPDPD(IRREPR,ISYTYP(1,LISTL))
        NUMDSX=IRPDPD(IRREPL,ISYTYP(2,LISTL))
        DISSYZ=IRPDPD(IRREPL,ISYTYP(1,LISTZ))
        NUMDSZ=IRPDPD(IRREPL,ISYTYP(2,LISTZ))
        I010=I000+IINTFP*NUMDSL*DISSYL
        I020=I010+IINTFP*MAX(NUMDSL*DISSYL,NUMDSX*DISSYX)
        CALL GETLST(ICORE(I000),1,NUMDSL,1,IRREPR,LISTL)
        CALL GETLST(ICORE(I010),1,NUMDSX,1,IRREPL,LISTSAV)
C
C FORM PRODUCT Y(aJ,Bi) = L(aJ,Fm) * X(Fm,Bi)
C
        CALL XGEMM('N','N',DISSYL,NUMDSX,NUMDSL,ONE,ICORE(I000),DISSYL,
     &             ICORE(I010),DISSYX,ZILCH,ICORE(I020),DISSYZ)
        CALL PUTLST(ICORE(I020),1,NUMDSZ,1,IRREPL,LISTZ)
130    CONTINUE
C
C NOW PROCESS INCREMENTS
C
       ISIZE=ISYMSZ(ISYTYP(1,42),ISYTYP(2,42))
       I000=1
       I010=I000+IINTFP*ISIZE
       I020=I010+IINTFP*ISIZE
C
C READ Y(aJ,Bi) FROM LIST 43 AND FORM
C
C        X(aJ,Bi) = Y(aJ,Bi) - 1/2 Y(BJ,ai)
C
       CALL GETALL(ICORE(I000),ISIZE,1,43)
       CALL SSTGEN(ICORE(I000),ICORE(I010),ISIZE,VRT(1,1),POP(1,1),
     &             VRT(1,1),POP(1,1),ICORE(I020),1,'3214')
       CALL SAXPY(ISIZE,HALF,ICORE(I010),1,ICORE(I000),1)
       CALL SSTGEN(ICORE(I000),ICORE(I010),ISIZE,VRT(1,1),POP(1,1),
     &             VRT(1,1),POP(1,1),ICORE(I020),1,'1432')
       CALL GETALL(ICORE(I000),ISIZE,1,42)
       CALL SAXPY(ISIZE,ONEM,ICORE(I010),1,ICORE(I000),1)
       IOFF=I000
       DO 200 IRREP=1,NIRREP
        DISSYZ=IRPDPD(IRREP,ISYTYP(1,42))
        NUMDSZ=IRPDPD(IRREP,ISYTYP(2,42))
        CALL MPMT(ICORE(IOFF),NUMDSZ)
        IOFF=IOFF+IINTFP*NUMDSZ*NUMDSZ
200    CONTINUE
       CALL SSTGEN(ICORE(I000),ICORE(I010),ISIZE,VRT(1,1),POP(1,1),
     &             VRT(1,1),POP(1,1),ICORE(I020),1,'1324')
       CALL GETALL(ICORE(I000),ISIZE,1,63)
       CALL SAXPY(ISIZE,ONE,ICORE(I010),1,ICORE(I000),1)
       CALL PUTALL(ICORE(I000),ISIZE,1,63)
C
      ELSE
C
       write(6,*)' sorry, uhf is not coded '
C
      ENDIF
C
      RETURN
      END
