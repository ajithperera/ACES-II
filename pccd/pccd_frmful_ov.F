      SUBROUTINE PCCD_FRMFUL_OV(DOV,DVO,DFULL,ISCR,MAXCOR,NBAS,
     +                          FUNC,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POP,VRT
      CHARACTER*4 FUNC

      DIMENSION DVO(1),DOV(1),DFULL(NBAS,NBAS),ISCR(MAXCOR)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM2/ POP(8,2),VRT(8,2),NDOVA,NDOVB,NDOOA,NDOOB,NDVVA,
     &              NDVVB
#include "sympop.com"

      DATA IONE /1/
      IOFFOV=1
      DO 100 ISPIN=1,1+IUHF
       LENVO=NVRTO(ISPIN)*NOCCO(ISPIN)
       I000=1
       I010=I000+LENVO*IINTFP
       I020=I010+LENVO
       IF(I020.GT.MAXCOR)CALL INSMEM('PCCD_FRMFUL_OV',I020,MAXCOR)
       IF (FUNC .EQ. "COPY") THEN
       CALL PUTBLK(DFULL,DVO(IOFFOV),'VO',ISCR(I000),ISCR(I010),NBAS,
     &             IRPDPD(1,8+ISPIN),LENVO,ISPIN)
       CALL PUTBLK(DFULL,DOV(IOFFOV),'OV',ISCR(I000),ISCR(I010),NBAS,
     &             IRPDPD(1,15+ISPIN),LENVO,ISPIN)
       ELSEIF (FUNC .EQ. "ADD ") THEN
       CALL ADDBLK(DFULL,DVO(IOFFOV),'VO',ISCR(I000),ISCR(I010),NBAS,
     &             IRPDPD(1,8+ISPIN),LENVO,ISPIN)
       CALL ADDBLK(DFULL,DOV(IOFFOV),'OV',ISCR(I000),ISCR(I010),NBAS,
     &             IRPDPD(1,15+ISPIN),LENVO,ISPIN)
       ENDIF 
       IOFFOV=IOFFOV+NDOVA
100   CONTINUE
      RETURN
      END 
