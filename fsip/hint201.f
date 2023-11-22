      SUBROUTINE HINT201(HEFF,ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE CALCULATES THE FOLLOWING CONTRIBUTIONS TO THE
C DOUBLES AMPLITUDES FOR THE (0,1) FOCK SPACE SECTOR
C
C
C    Z(ij,ka) = - T(ij,ma) * Heff(km)
C      NN AN        NN AN         AA
C
CEND
      IMPLICIT INTEGER (A-H,O-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      LOGICAL RHF
      DIMENSION ICORE(MAXCOR),HEFF(1)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /FSSYMPOP/ FSDPDAN(8,22),FSDPDNA(8,22),FSDPDAA(8,22),
     &                  FSDPDIN(8,22),FSDPDNI(8,22),FSDPDII(8,22),
     &                  FSDPDAI(8,22),FSDPDIA(8,22)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /FSSYM/ POPA(8,2),POPI(8,2),VRTA(8,2),VRTI(8,2),
     &               NTAN(2),NTNA(2),NTAA(2),
     &               NTIN(2),NTNI(2),NTII(2),
     &               NTAI(2),NTIA(2),
     &               NF1AN(2),NF1AA(2),NF1IN(2),NF1II(2),NF1AI(2),
     &               NF2AN(2),NF2AA(2),NF2IN(2),NF2II(2),NF2AI(2)
C
      DATA ONE/1.0/, ONEM/-1.0/, ZILCH /0.0/
C
      RHF=IUHF.EQ.0
C
      DO 10 IRREP=1,NIRREP
C
C SPIN CASE ABAB -> T(Ij,Ma)*Heff(MK)
C                     NN AN       AA
C
       LISTT=99
       LISTZ=199
       NUMDST=FSDPDAN(IRREP,ISYTYP(2,LISTT))
       DISSYT=IRPDPD (IRREP,ISYTYP(1,LISTT))
       I000=1
       I010=I000+IINTFP*NUMDST*DISSYT
       I020=I010+IINTFP*NUMDST*DISSYT
       MAXT=NUMDST*DISSYT
       I030=I020+IINTFP*MAXT
       I040=I030+IINTFP*MAXT
       I050=I040+IINTFP*MAXT
       CALL FSGET (ICORE(I000),1,NUMDST,1,IRREP,LISTT,'NNAN')
       CALL FSGETT1(HEFF,3,91,'AA',21)
       CALL SYMTR1(IRREP,POPA(1,1),VRT(1,2),DISSYT,ICORE(I000),
     &             ICORE(I020),ICORE(I030),ICORE(I040))
C
C DO MULTIPLICATION T(Ija,M)*HEFF(KM)
C                     NNN A       AA
C
       IOFFT=I000
       IOFFZ=I010
       IOFFH=1
       DO 20 IRREPM=1,NIRREP
        IRREPK=IRREPM
        IRREPA=DIRPRD(IRREPM,IRREP)
        NROW=DISSYT*VRT(IRREPA,2)
        NCOL=POPA(IRREPK,1)
        NSUM=POPA(IRREPM,1)
        CALL XGEMM('N','T',NROW,NCOL,NSUM,ONEM,ICORE(IOFFT),NROW,
     &             HEFF(IOFFH),NSUM,ZILCH,ICORE(IOFFZ),NROW)
        IOFFT=IOFFT+NROW*NSUM*IINTFP
        IOFFZ=IOFFZ+NROW*NCOL*IINTFP
        IOFFH=IOFFH+NSUM*NCOL*IINTFP
20     CONTINUE
       CALL SYMTR1(IRREP,VRT(1,2),POPA(1,1),DISSYT,ICORE(I010),
     &             ICORE(I020),ICORE(I030),ICORE(I040))
       CALL FSGET(ICORE(I000),1,NUMDST,1,IRREP,LISTZ,'NNAN')
       CALL SAXPY(NUMDST*DISSYT,ONE,ICORE(I010),1,ICORE(I000),1)
       CALL FSPUT(ICORE(I000),1,NUMDST,1,IRREP,LISTZ,'NNAN')
C
       IF(RHF)GOTO 10
C
C SPIN CASE BABA -> T(Ij,Am)*Heff(km)
C                     NN NA       AA
C
       LISTT=98
       LISTZ=198
       NUMDST=FSDPDNA(IRREP,ISYTYP(2,LISTT))
       DISSYT=IRPDPD (IRREP,ISYTYP(1,LISTT))
       I000=1
       I010=I000+IINTFP*NUMDST*DISSYT
       I020=I010+IINTFP*NUMDST*DISSYT
       MAXT=NUMDST*DISSYT
       I030=I020+IINTFP*MAXT
       I040=I030+IINTFP*MAXT
       I050=I040+IINTFP*MAXT
       CALL FSGET (ICORE(I000),1,NUMDST,1,IRREP,LISTT,'NNNA')
       CALL FSGETT1(HEFF,4,91,'AA',22)
C
C DO MULTIPLICATION T(IjA,m)*HEFF(mk)
C                     NNN A       AA
C
       IOFFT=I000
       IOFFZ=I010
       IOFFH=1
       DO 30 IRREPM=1,NIRREP
        IRREPK=IRREPM
        IRREPA=DIRPRD(IRREPM,IRREP)
        NROW=DISSYT*VRT(IRREPA,1)
        NCOL=POPA(IRREPK,2)
        NSUM=POPA(IRREPM,2)
        CALL XGEMM('N','T',NROW,NCOL,NSUM,ONEM,ICORE(IOFFT),NROW,
     &             HEFF(IOFFH),NSUM,ZILCH,ICORE(IOFFZ),NROW)
        IOFFT=IOFFT+NROW*NSUM*IINTFP
        IOFFZ=IOFFZ+NROW*NCOL*IINTFP
        IOFFH=IOFFH+NSUM*NCOL*IINTFP
30     CONTINUE
       CALL FSGET(ICORE(I000),1,NUMDST,1,IRREP,LISTZ,'NNNA')
       CALL SAXPY(NUMDST*DISSYT,ONE,ICORE(I010),1,ICORE(I000),1)
       CALL FSPUT(ICORE(I000),1,NUMDST,1,IRREP,LISTZ,'NNNA')
C
C SPIN CASES AAAA AND BBBB -> T(I<J,MA)*Heff(KM)
C                               N N AN       AA
C
       DO 100 ISPIN=1,2
        LISTT=95+ISPIN
        LISTZ=195+ISPIN
        NUMDST=FSDPDAN(IRREP,ISYTYP(2,LISTT))
        DISSYT=IRPDPD (IRREP,ISYTYP(1,LISTT))
        I000=1
        I010=I000+IINTFP*NUMDST*DISSYT
        I020=I010+IINTFP*NUMDST*DISSYT
        MAXT=MAX(NUMDST,DISSYT)
        I030=I020+IINTFP*MAXT
        I040=I030+IINTFP*MAXT
        I050=I040+IINTFP*MAXT
        CALL FSGET (ICORE(I000),1,NUMDST,1,IRREP,LISTT,'NNAN')
        CALL SYMTR1(IRREP,POPA(1,ISPIN),VRT(1,ISPIN),DISSYT,ICORE(I000),
     &             ICORE(I020),ICORE(I030),ICORE(I040))
        CALL FSGETT1(HEFF,ISPIN+2,91,'AA',20+ISPIN)
C
C DO MULTIPLICATION T(I<JA,M)*HEFF(KM)
C                     N NN A       AA
C
        IOFFT=I000
        IOFFZ=I010
        IOFFH=1
        DO 40 IRREPM=1,NIRREP
         IRREPK=IRREPM
         IRREPA=DIRPRD(IRREPM,IRREP)
         NROW=DISSYT*VRT(IRREPA,ISPIN)
         NCOL=POPA(IRREPK,ISPIN)
         NSUM=POPA(IRREPM,ISPIN)
         CALL XGEMM('N','T',NROW,NCOL,NSUM,ONEM,ICORE(IOFFT),NROW,
     &              HEFF(IOFFH),NSUM,ZILCH,ICORE(IOFFZ),NROW)
         IOFFT=IOFFT+NROW*NSUM*IINTFP
         IOFFZ=IOFFZ+NROW*NCOL*IINTFP
         IOFFH=IOFFH+NSUM*NCOL*IINTFP
40      CONTINUE
        CALL SYMTR1(IRREP,VRT(1,ISPIN),POPA(1,ISPIN),DISSYT,ICORE(I010),
     &              ICORE(I020),ICORE(I030),ICORE(I040))
        CALL FSGET(ICORE(I000),1,NUMDST,1,IRREP,LISTZ,'NNAN')
        CALL SAXPY(NUMDST*DISSYT,ONE,ICORE(I010),1,ICORE(I000),1)
        CALL FSPUT(ICORE(I000),1,NUMDST,1,IRREP,LISTZ,'NNAN')
100    CONTINUE
C
10    CONTINUE
C
      RETURN
      END
