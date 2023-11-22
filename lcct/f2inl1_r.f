      SUBROUTINE F2INL1_R(ICORE,MAXCOR,IUHF,IOFFLIST,LISTL1,LISTL1OFF)
C
C THIS ROUTINE CALCULATES THE TERM
C
C
C - SUM M L(M,A) F(IM)
C
C - SUM m L(m,a) F(im)   (UHF only)
C
CEND
C
C CODED JG AUGUST/90
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,DISSYZ,POP,VRT
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWO
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,
     &             NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C
C
      DATA ONE,ONEM,HALF/1.0D0,-1.0D0,0.5D0/
C
C   CALCULATE SIZE OF F(IM) ARRAY
C
      NFAA=NF1AA
      NFBB=NF1BB
      I0AA=MAXCOR+1-NFAA*IINTFP
      MXCOR=MAXCOR-NFAA*IINTFP
      CALL GETLST(ICORE(I0AA),1,1,1,1,91+IOFFLIST)
      IF(IUHF.EQ.0) THEN
       I0BB=I0AA
      ELSE
       I0BB=I0AA-NFBB*IINTFP
       MXCOR=MXCOR-NFBB*IINTFP
       CALL GETLST(ICORE(I0BB),1,1,1,2,91+IOFFLIST)
      ENDIF
C
C ADDITIONAL CODE FOR QCISD AND CCSD METHODS
C
        I0TA=I0BB-NTAA*IINTFP
        I0ZA=I0TA-NTAA*IINTFP
        MXCOR=MXCOR-2*NTAA*IINTFP
        CALL GETLST(ICORE(I0TA),1,1,2,1+LISTL1OFF,LISTL1)
        CALL GETLST(ICORE(I0ZA),1,1,1,3,90)
        IF(IUHF.EQ.0) THEN
         I0TB=I0TA
         I0ZB=I0ZA
        ELSE
         I0TB=I0ZA-NTBB*IINTFP
         I0ZB=I0TB-NTBB*IINTFP
         MXCOR=MXCOR-2*NTBB*IINTFP
         CALL GETLST(ICORE(I0TB),1,1,2,2+LISTL1OFF,LISTL1)
         CALL GETLST(ICORE(I0ZB),1,1,1,4,90)
        ENDIF


       DO 300 ISPIN=1,IUHF+1
C
        IF(ISPIN.EQ.1) THEN
         IOFFF=I0AA
         IOFFT=I0TA
         IOFFZ=I0ZA
         I0Z=I0ZA
        ELSE
         IOFFF=I0BB
         IOFFT=I0TB
         IOFFZ=I0ZB
         I0Z=I0ZB
        ENDIF
        DO 250 IRREP=1,NIRREP
C
         NOCC=POP(IRREP,ISPIN)
         NVRT=VRT(IRREP,ISPIN)
         IF(MIN(NVRT,NOCC).GT.0) THEN
         CALL XGEMM('N','T',NVRT,NOCC,NOCC,ONEM,ICORE(IOFFT),NVRT,
     &              ICORE(IOFFF),NOCC,ONE,ICORE(IOFFZ),NVRT)
         ENDIF
         IOFFF=IOFFF+NOCC*NOCC*IINTFP
         IOFFT=IOFFT+NOCC*NVRT*IINTFP
         IOFFZ=IOFFZ+NOCC*NVRT*IINTFP
250     CONTINUE
        CALL PUTLST(ICORE(I0Z),1,1,1,ISPIN+2,90)
300    CONTINUE

      RETURN
      END
