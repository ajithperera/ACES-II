      SUBROUTINE HHL01(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE CALCULATES THE FOLLOWING CONTRIBUTION TO THE
C DOUBLES AMPLITUDES FOR THE (0,1) FOCK SPACE SECTOR
C
C
C    Z(ij,ka) = T(mn,ka) * W(mn,ij)
C      NN AN      NN AN      NN NN      
C
C THE WORLD'S EASIEST CODE  
C
CEND
      IMPLICIT INTEGER (A-H,O-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      LOGICAL RHF
      DIMENSION ICORE(MAXCOR)
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
C SPIN CASE ABAB -> T(Mn,Ka)*W(Mn,Ij)
C
       LISTT=99
       LISTW=53
       LISTZ=199
       NUMDSW=IRPDPD (IRREP,ISYTYP(2,LISTW))
       DISSYW=IRPDPD (IRREP,ISYTYP(1,LISTW))
       NUMDST=FSDPDAN(IRREP,ISYTYP(2,LISTT))
       DISSYT=IRPDPD (IRREP,ISYTYP(1,LISTT))
       I000=1
       I010=I000+IINTFP*NUMDSW*DISSYW
       I020=I010+IINTFP*NUMDST*DISSYT
       I030=I020+IINTFP*NUMDST*DISSYT
       CALL GETLST(ICORE(I000),1,NUMDSW,1,IRREP,LISTW)
       CALL FSGET (ICORE(I010),1,NUMDST,1,IRREP,LISTT,'NNAN')
       NROW=DISSYW
       NCOL=NUMDST
       NSUM=DISSYW
       CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(I000),NSUM,
     &            ICORE(I010),NSUM,ZILCH,ICORE(I020),NROW)
       CALL FSGET(ICORE(I010),1,NUMDST,1,IRREP,LISTZ,'NNAN')
       CALL SAXPY(NUMDST*DISSYT,ONE,ICORE(I020),1,ICORE(I010),1)
       CALL FSPUT(ICORE(I010),1,NUMDST,1,IRREP,LISTZ,'NNAN')
C
C SKIP REMAINING SPIN CASES FOR RHF
C
       IF(RHF)GOTO 10
C
C SPIN CASE BABA -> T(Mn,Ak)*W(Mn,Ij)
C
       LISTT=98
       LISTW=53
       LISTZ=198
       NUMDSW=IRPDPD (IRREP,ISYTYP(2,LISTW))
       DISSYW=IRPDPD (IRREP,ISYTYP(1,LISTW))
       NUMDST=FSDPDNA(IRREP,ISYTYP(2,LISTT))
       DISSYT=IRPDPD (IRREP,ISYTYP(1,LISTT))
       I000=1
       I010=I000+IINTFP*NUMDSW*DISSYW
       I020=I010+IINTFP*NUMDST*DISSYT
       I030=I020+IINTFP*NUMDST*DISSYT
       CALL GETLST(ICORE(I000),1,NUMDSW,1,IRREP,LISTW)
       CALL FSGET (ICORE(I010),1,NUMDST,1,IRREP,LISTT,'NNNA')
       NROW=DISSYW
       NCOL=NUMDST
       NSUM=DISSYW
       CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(I000),NSUM,
     &            ICORE(I010),NSUM,ZILCH,ICORE(I020),NROW)
       CALL FSGET(ICORE(I010),1,NUMDST,1,IRREP,LISTZ,'NNNA')
       CALL SAXPY(NUMDST*DISSYT,ONE,ICORE(I020),1,ICORE(I010),1)
       CALL FSPUT(ICORE(I010),1,NUMDST,1,IRREP,LISTZ,'NNNA')
C
C SPIN CASES AAAA AND BBBB -> T(M<N,KA)*W(M<N,I<J)
C
       DO 100 ISPIN=1,2
        LISTT=95+ISPIN
        LISTW=50+ISPIN
        LISTZ=195+ISPIN
        NUMDSW=IRPDPD (IRREP,ISYTYP(2,LISTW))
        DISSYW=IRPDPD (IRREP,ISYTYP(1,LISTW))
        NUMDST=FSDPDAN(IRREP,ISYTYP(2,LISTT))
        DISSYT=IRPDPD (IRREP,ISYTYP(1,LISTT))
        I000=1
        I010=I000+IINTFP*NUMDSW*DISSYW
        I020=I010+IINTFP*NUMDST*DISSYT
        I030=I020+IINTFP*NUMDST*DISSYT
        CALL GETLST(ICORE(I000),1,NUMDSW,1,IRREP,LISTW)
        CALL FSGET (ICORE(I010),1,NUMDST,1,IRREP,LISTT,'NNAN')
        NROW=DISSYW
        NCOL=NUMDST
        NSUM=DISSYW
        CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(I000),NSUM,
     &             ICORE(I010),NSUM,ZILCH,ICORE(I020),NROW)
        CALL FSGET(ICORE(I010),1,NUMDST,1,IRREP,LISTZ,'NNAN')
        CALL SAXPY(NUMDST*DISSYT,ONE,ICORE(I020),1,ICORE(I010),1)
        CALL FSPUT(ICORE(I010),1,NUMDST,1,IRREP,LISTZ,'NNAN')
100    CONTINUE
C
10    CONTINUE
C
      RETURN
      END