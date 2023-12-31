      SUBROUTINE RSL1INL2(ICORE,MAXCOR,IUHF,IRREPX,LISTT2RS)
C
C This subroutine calculates the terms 
C        L(AI,BJ) = L(AI) * F(BJ)
C
C These terms are usually calculated in DRHFRNG and DUHFRNG, but
C   this routine is for the times when those routines are not called.
C
C      IMPLICIT NONE
      INTEGER ICORE, MAXCOR, IUHF, IRREPX, LISTT2RS
      DIMENSION ICORE(MAXCOR)
C
      INTEGER IINTLN, IFLTLN, IINTFP, IALONE, IBITWD
      INTEGER POP, VRT, NT, NFMI, NFEA, NSTART, NIRREP, IRREPS, DIRPRD
      INTEGER IRPDPD, ISYTYP, ID
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      DOUBLE PRECISION ONE, ONEM, ZILCH
      PARAMETER (ONE = 1.0D0, ONEM = -1.0D0, ZILCH = 0.0D0)
      INTEGER I000, I010, I020, ISPIN, NCOL, NROW, IRREPL, IRREPR
      INTEGER LSTTMP, DISSYZ, NUMDSZ
C
      IRREPR = 1
      IRREPL=DIRPRD(IRREPR,IRREPX)
      IF (IUHF .EQ. 0) THEN
C
C  do the contribution for RHF
C        L(AI,BJ) = L(AI) * F(BJ)
C
        LSTTMP=LISTT2RS + 2
        DISSYZ=IRPDPD(IRREPL,ISYTYP(1,LSTTMP))
        NUMDSZ=IRPDPD(IRREPR,ISYTYP(2,LSTTMP))
        I000=1
        I010=I000+IINTFP*DISSYZ*NUMDSZ
        I020=I010+IINTFP*NT(1)
        CALL GETLST(ICORE(I010),1,1,1,1,93)
        CALL GETLST(ICORE(I020),1,1,1,1,490)
        NROW=IRPDPD(IRREPX,9)
        NCOL=NT(1)
        CALL XGEMM('N','N',NROW,NCOL,1,ONEM,ICORE(I020),NROW,
     &     ICORE(I010),1,ZILCH,ICORE(I000),NROW)
C
C symmetrize the list
C
C This is a special case of the symmetrization from DRHFRNG.  In
C   the present case, we know that the only non-zero contribution
C   is for irrepr=1.  Therefore, unless irrepx=1, the irrepl list
C   is zero, and we do not read and add it.
C
        IF(IRREPX.EQ.1)THEN
         CALL MPMT  (ICORE(I000),DISSYZ)
         CALL PUTLST(ICORE(I000),1,DISSYZ,1,IRREPR,LSTTMP)
        ELSE
         CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPR,LSTTMP)
         CALL TRANSP(ICORE(I000),ICORE(I010),NUMDSZ,DISSYZ)
         CALL PUTLST(ICORE(I010),1,DISSYZ,1,IRREPL,LSTTMP)
        ENDIF
      ELSE
C
C do the contribution for UHF
C        L(AI,BJ) = L(AI) * F(BJ)
C
        DO 10 ISPIN = 1,2
          LSTTMP=LISTT2RS-1+ISPIN
          DISSYZ=IRPDPD(IRREPL,ISYTYP(1,LSTTMP))
          NUMDSZ=IRPDPD(IRREPR,ISYTYP(2,LSTTMP))
          I000=1
          I010=I000+IINTFP*DISSYZ*NUMDSZ
          I020=I010+IINTFP*NT(ISPIN)
          CALL GETLST(ICORE(I010),1,1,1,ISPIN,93)
          CALL GETLST(ICORE(I020),1,1,1,ISPIN,490)
          NROW=IRPDPD(IRREPX,8+ISPIN)
          NCOL=NT(ISPIN)
          CALL XGEMM('N','N',NROW,NCOL,1,ONE,ICORE(I020),NROW,
     &       ICORE(I010),1,ZILCH,ICORE(I000),NROW)
          CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPR,LSTTMP)
 10     CONTINUE
C
C        L(AI,bj) = L(AI) * F(bj)
C
        LSTTMP = LISTT2RS + 2
        DISSYZ=IRPDPD(IRREPL,ISYTYP(1,LSTTMP))
        NUMDSZ=IRPDPD(IRREPR,ISYTYP(2,LSTTMP))
        I000=1
        I010=I000+IINTFP*DISSYZ*NUMDSZ
        I020=I010+IINTFP*NT(2)
        CALL GETLST(ICORE(I010),1,1,1,2,93)
        CALL GETLST(ICORE(I020),1,1,1,1,490)
        NROW=IRPDPD(IRREPX,9)
        NCOL=NT(2)
        CALL XGEMM('N','N',NROW,NCOL,1,ONEM,ICORE(I020),NROW,
     &             ICORE(I010),1,ZILCH,ICORE(I000),NROW)
        IF (IRREPL .NE. 1) CALL PUTLST(ICORE(I000),1,NUMDSZ,1,
     &     IRREPR,LSTTMP)
C
C        L(AI,bj) = F(AI) * L(bj)
C
        IRREPL = 1
        IRREPR=DIRPRD(IRREPL,IRREPX)
        DISSYZ=IRPDPD(IRREPL,ISYTYP(1,LSTTMP))
        NUMDSZ=IRPDPD(IRREPR,ISYTYP(2,LSTTMP))
        I000=1
        I010=I000+IINTFP*DISSYZ*NUMDSZ
        I020=I010+IINTFP*IRPDPD(IRREPX,10)
        CALL GETLST(ICORE(I010),1,1,1,2,490)
        CALL GETLST(ICORE(I020),1,1,1,1,93)
        NROW=NT(1)
        NCOL=IRPDPD(IRREPX,10)
        IF (IRREPR .EQ. 1) THEN
          CALL XGEMM('N','N',NROW,NCOL,1,ONEM,ICORE(I020),NROW,
     &       ICORE(I010),1,ONE,ICORE(I000),NROW)
        ELSE
          CALL XGEMM('N','N',NROW,NCOL,1,ONEM,ICORE(I020),NROW,
     &       ICORE(I010),1,ZILCH,ICORE(I000),NROW)
        ENDIF
        CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPR,LSTTMP)
      ENDIF
      RETURN
      END
