
      SUBROUTINE PCCD_VVINOO(IOO,DVV,ICORE,MXCOR,IUHF)
C
C THIS ROUTINE COMPUTES THE CONTRIBUTIONS:
C
C      IOO(I,J) = - [SUM D(E,F) * <EI||FJ> + D(e,f) * <eI|fJ>] (ISPIN=1)
C
C      IOO(i,j) = - [SUM D(e,f) * <ei||fj> + D(E,F) * <Ei|Fj>] (ISPIN=2)
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION IOO,DVV,ONE
      DIMENSION DVV(1),IOO(1),ICORE(MXCOR)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD

      DATA ONE /1.0/
C
      IF(IUHF.EQ.0)THEN
C
C SPIN ADAPTED RHF CODE
C
       LISTW1=25 
       LISTW2=23 
       DO 120 IRREP=1,NIRREP
        DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW1))
        NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW2))
        I000=1
        I010=I000+IINTFP*NUMDIS*DISSIZ
        I020=I010+IINTFP*NUMDIS*DISSIZ
        CALL GETLST(ICORE(I000),1,NUMDIS,1,IRREP,LISTW1)
        CALL GETLST(ICORE(I010),1,NUMDIS,1,IRREP,LISTW2)
        CALL SAXPY(NUMDIS*DISSIZ,ONE,ICORE(I010),1,ICORE(I000),1)
        CALL VMINUS(ICORE(I000),NUMDIS*DISSIZ)
        CALL DOT24(IRREP,IOO,DVV,ICORE(I000),ICORE(I010),DISSIZ,
     &             POP(1,1),POP(1,1),VRT(1,1),POP(1,1),VRT(1,1),
     &             POP(1,1),'STST')
120    CONTINUE
      ELSEIF(IUHF.NE.0)THEN
       DO 10 ISPIN=1,1+IUHF
        LISTW1=22+ISPIN 
        LISTW2=27-ISPIN 
        DO 20 IRREP=1,NIRREP
         NDSZ1=IRPDPD(IRREP,ISYTYP(1,LISTW1))
         NDIS1=IRPDPD(IRREP,ISYTYP(2,LISTW1))
         NDSZ2=IRPDPD(IRREP,ISYTYP(1,LISTW2))
         NDIS2=IRPDPD(IRREP,ISYTYP(2,LISTW2))
         IOFFI=1+(ISPIN-1)*NFMI(1)
C
C DO FIRST PART - D(E,F) * <EI||FJ> [D(e,f) * <ei||fj>]
C
         I000=1
         I010=I000+IINTFP*NDSZ1*NDIS1
         I020=I010+IINTFP*NFEA(ISPIN)
         IOFFD=1+(ISPIN-1)*NFEA(1)
         CALL GETLST(ICORE(I000),1,NDIS1,2,IRREP,LISTW1)
         CALL VMINUS(ICORE(I000),NDIS1*NDSZ1)
         CALL DOT24(IRREP,IOO(IOFFI),DVV(IOFFD),ICORE(I000),
     &              ICORE(I010),NDSZ1,POP(1,ISPIN),POP(1,ISPIN),
     &              VRT(1,ISPIN),POP(1,ISPIN),VRT(1,ISPIN),POP(1,ISPIN),
     &              'STST')
C
C NOW DO SECOND PART 
C
C          D(e,f) * <eI|fJ>   [ISPIN=1]
C
C          D(E,F) * <Ei|Fj>   [ISPIN=2]
C
         I000=1
         I010=I000+IINTFP*NDSZ2*NDIS2
         I020=I010+NFEA(3-ISPIN)
         IOFFD=1+(2-ISPIN)*NFEA(1)
         CALL GETLST(ICORE(I000),1,NDIS2,2,IRREP,LISTW2)
         CALL VMINUS(ICORE(I000),NDIS2*NDSZ2)
         IF(ISPIN.EQ.1)THEN
          CALL DOT24(IRREP,IOO(IOFFI),DVV(IOFFD),ICORE(I000),
     &               ICORE(I010),NDSZ2,POP(1,ISPIN),POP(1,ISPIN),
     &               VRT(1,2),POP(1,1),VRT(1,2),POP(1,1),'STST')
         ELSE
          CALL DOT24(IRREP,IOO(IOFFI),DVV(IOFFD),ICORE(I000),
     &               ICORE(I010),NDSZ2,POP(1,ISPIN),POP(1,ISPIN),
     &               VRT(1,1),POP(1,2),VRT(1,1),POP(1,2),'STST')
         ENDIF
20      CONTINUE
10     CONTINUE
      ENDIF
      RETURN
      END
