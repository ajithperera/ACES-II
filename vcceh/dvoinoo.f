      SUBROUTINE DVOINOO(IRREPX,IOO,DVO,ICORE,MXCOR,IUHF,LISTW0)
C
C THIS ROUTINE COMPUTES THE CONTRIBUTIONS:
C
C IOO(I,J) =  SUM [D(M,E) * <IM||JE> + D(M,E) * <IM|JE>] (ISPIN=1)
C                                             
C
C IOO(I,J) =  SUM [D(M,E) * <IM||JE> + D(M,E) * <MI|EJ>] (ISPIN=1)
C                                             
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION IOO,DVO,ONE
      DIMENSION DVO(1),IOO(1),ICORE(MXCOR)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      DATA ONE /1.0/
C
      IF(IUHF.EQ.0)THEN
C
C SPIN ADAPTED RHF CODE
C      
       IOFFI=1
       LISTW=LISTW0-1+4
       LENTAR=IRPDPD(IRREPX,21)
       DO 120 IRREP=1,NIRREP
        DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
        NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
        I000=1
        I010=I000+IINTFP*LENTAR
        I020=I010+IINTFP*NUMDIS*DISSIZ
        I030=I020+IINTFP*NUMDIS
        I040=I030+IINTFP*NUMDIS
        CALL IZERO(ICORE(I000),IINTFP*LENTAR)
        CALL GETLST(ICORE(I010),1,NUMDIS,1,IRREP,LISTW)
        CALL SPINAD3(IRREP,POP(1,1),DISSIZ,NUMDIS,ICORE(I010),
     &               ICORE(I020),ICORE(I030))
        CALL DDOT24(IRREPX,IRREP,ICORE(I000),DVO,ICORE(I010),
     &             ICORE(I020),DISSIZ,POP(1,1),POP(1,1),
     &             POP(1,1),POP(1,1),POP(1,1),VRT(1,1),'TSTS')
        CALL SAXPY(LENTAR,ONE,ICORE(I000),1,IOO(IOFFI),1)
120    CONTINUE
      ELSEIF(IUHF.NE.0)THEN
       DO 10 ISPIN=1,1+IUHF
        LISTW1=LISTW0-1+ISPIN
        LISTW2=LISTW0-1+5-ISPIN
        LENTAR=IRPDPD(IRREPX,20+ISPIN)
        LEND  =IRPDPD(IRREPX,8+ISPIN)
        DO 20 IRREP=1,NIRREP
         NDSZ1=IRPDPD(IRREP,ISYTYP(1,LISTW1))
         NDIS1=IRPDPD(IRREP,ISYTYP(2,LISTW1))
         NDSZ1F=IRPDPD(IRREP,20+ISPIN)
         NDSZ2=IRPDPD(IRREP,ISYTYP(1,LISTW2))
         NDIS2=IRPDPD(IRREP,ISYTYP(2,LISTW2))
         IOFFI=1+(ISPIN-1)*IRPDPD(IRREPX,21)
C
C DO FIRST PART   D(M,E) * <IM||JE>  (ISPIN=1)
C                 D(M,E) * <IM||JE>  (ISPIN=2)
C
C AND KEEP RESULT AT BOTTOM OF CORE
C
         I000=1
         I010=I000+IINTFP*LENTAR
         I020=I010+IINTFP*NDSZ1F*NDIS1
         I030=I020+IINTFP*LEND
         CALL IZERO(ICORE(I000),LENTAR*IINTFP)
         IOFFD=1+(ISPIN-1)*IRPDPD(IRREPX,9)
         CALL GETLST(ICORE(I010),1,NDIS1,2,IRREP,LISTW1)
         CALL SYMEXP2(IRREP,POP(1,ISPIN),NDSZ1F,NDSZ1,NDIS1,ICORE(I010),
     &                ICORE(I010))
         CALL DDOT24(IRREPX,IRREP,ICORE(I000),DVO(IOFFD),ICORE(I010),
     &              ICORE(I020),NDSZ1F,POP(1,ISPIN),POP(1,ISPIN),
     &              POP(1,ISPIN),POP(1,ISPIN),POP(1,ISPIN),VRT(1,ISPIN),
     &              'TSTS')
C
C NOW DO SECOND PART 
C
C         D(M,E) * <IM|JE> (ISPIN = 1)
C         D(M,E) * <MI|EJ> (ISPIN = 2)
C
         I020=I010+IINTFP*NDSZ2*NDIS2
         I030=I020+IRPDPD(IRREPX,11-ISPIN)
         IOFFD=1+(2-ISPIN)*IRPDPD(IRREPX,9)
         CALL GETLST(ICORE(I010),1,NDIS2,2,IRREP,LISTW2)
         IF(ISPIN.EQ.1)THEN
          CALL DDOT24(IRREPX,IRREP,ICORE(I000),DVO(IOFFD),ICORE(I010),
     &               ICORE(I020),NDSZ2,POP(1,1),POP(1,1),
     &               POP(1,1),POP(1,2),POP(1,1),VRT(1,2),'TSTS')
         ELSE
          CALL DDOT24(IRREPX,IRREP,ICORE(I000),DVO(IOFFD),ICORE(I010),
     &               ICORE(I020),NDSZ2,POP(1,2),POP(1,2),
     &               POP(1,1),POP(1,2),VRT(1,1),POP(1,2),'STST')
         ENDIF
C
C NOW FORM TARGET(I,J)+TARGET(J,I) [I,J;J,I] AND DECREMENT IOO.
C
         CALL SAXPY(LENTAR,ONE,ICORE(I000),1,IOO(IOFFI),1)
20      CONTINUE
10     CONTINUE
      ENDIF
      RETURN
      END
