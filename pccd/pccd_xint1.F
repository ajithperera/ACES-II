      SUBROUTINE PCCD_XINT1(XIA,DOO,ICORE,MAXCOR,IUHF)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,DISSYW,POP,VRT
      DIMENSION XIA(1),DOO(1),ICORE(MAXCOR)
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      MXCOR=MAXCOR 
      ISHIFT=0
C
      DO 1000 ISPIN=1,IUHF+1
C
       IF(ISPIN.EQ.1) THEN
        IOFFX=1
        IOFFDA=1
        IOFFDB=1+NF1AA*IUHF
       ELSE
        IOFFX=1+NTAA
        IOFFDA=1+NF1AA
        IOFFDB=1
       ENDIF
C
C  FIRST DO AAAA PART (UHF ONLY)
C
       IF(IUHF.NE.0) THEN
C
        LISTW=6+ISPIN + ISHIFT 
C
        DO 100 IRREP=1,NIRREP
C
         NOCCSQ=0
         DO 101 IRREPJ=1,NIRREP
          NOCCSQ=NOCCSQ+POP(IRREPJ,ISPIN)*
     &                  POP(DIRPRD(IRREP,IRREPJ),ISPIN)
101      CONTINUE
C
         NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
         DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
C
         I001=1
         I002=I001+IINTFP*NOCCSQ*NUMSYW
         IF(I002.LT.MXCOR) THEN
C
          CALL XIA1AA(ICORE(I001),DOO(IOFFDA),XIA(IOFFX),ISPIN,
     &                POP(1,ISPIN),VRT(1,ISPIN),DISSYW,NUMSYW,
     &                NOCCSQ,LISTW,IRREP)
         ELSE   
          CALL INSMEM('pccd_xint1',I002,MXCOR)
         ENDIF
100     CONTINUE
C
       ENDIF
C 
       LISTW=8+ISPIN + ISHIFT 
       IF(IUHF.EQ.0) LISTW=10 + ISHIFT 
C
       DO 200 IRREP=1,NIRREP
C
        NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
        DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
C
        I001=1
        I002=I001+IINTFP*NUMSYW*DISSYW
        I003=I002
        IF(IUHF.EQ.0) I003=I002+2*IINTFP*NUMSYW
        IF(I003.LT.MXCOR) THEN
C
         CALL XIA1AB(ICORE(I001),DOO(IOFFDB),XIA(IOFFX),ISPIN,
     &               POP(1,ISPIN),VRT(1,ISPIN),POP(1,3-ISPIN),
     &               VRT(1,3-ISPIN),DISSYW,NUMSYW,LISTW,IRREP,
     &               IUHF,ICORE(I002))
        ELSE
         CALL INSMEM('pccd_xint1',I003,MXCOR)
        ENDIF
200    CONTINUE
1000  CONTINUE
      RETURN
      END
