      SUBROUTINE TPDIJAB4(ICORE,MAXCOR,IUHF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL  MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      INTEGER DIRPRD,DISSYR,DISSYT,DISSYL,DISSYG,POP,VRT
      DIMENSION ICORE(MAXCOR)
      DIMENSION I0T(2),I0R(2)
C
C CALCULATION OF THE FOURTH GROUP OF IJAB CONTRIBUTIONS TO
C THE EOM-CCSD TWO-PARTICLE DENSITY MATRIX
C
C     1/16 R(IJ,EF) L(EF,MN) (T(MN,AB) + T(M,A) T(N,B) - T(M,B) T(N,A))
C
C +   1/16 R(MN,AB) L(EF,MN) (T(IJ,EF) + T(I,E) T(J,F) - T(I,F) T(J,E))
C
C +   1/16 (R(I,E) T(J,F) + T(I,E) R(J,F) - R(I,F) T(J,E) - R(J,E) T(I,F)
C
C            L(EF,MN) (T(MN,AB) + T(M,A) T(N,B) - T(,M,B) T(N,A))
C  
C +   1/16 (R(M,A) T(N,B) + T(M,A) R(N,B) - R(M,B) T(N,A) - R(N,A) T(M,B)
C
C            L(EF,MN) (T(IJ,EF) + T(I,E) T(J,F) - T(I,F) T(J,E))
C
C  WHICH GIVES IN THE MOST COMPACT NOTATION (SUITABLE FOR CODING):
C
C   = 1/16 RTAU(IJ,EF) L(EF,MN) TAU(MN,AB)
C
C     + 1/16 TAU(IJ,EF) L(EF,MN) RTAU(MN,AB)
C
C   WITH TAU(IJ,EF) = T(IJ,EF) + T(I,E) T(J,F) - T(J,E) T(I,F)
C
C   AND  RTA(IJ,EF) = R(IJ,EF) + R(I,E) T(J,F) + T(I,E) R(J,F) 
C
C                              - R(I,F) T(J,E) - R(I,E) T(J,F)
C
C
C NOTE THAT ACTUALLY 4*G(IJ,AB) IS CALCULATED WHICH MUST BE
C DIVIDED BY TWO AFTER ALL TERMS ARE SUMMED UP.
C
CEND 
C
C CODED JG OCTOBER/93
C
      COMMON/STATSYM/IRREPX
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/REFTYPE/MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
C
      DATA AZERO,ONE,TWO /0.D0,1.0D0,2.D0/
C
C READ IN T1 AND R1 AMPLITUDES
C
      I0T(1)=1
      I0T(2)=I0T(1)+IRPDPD(1,9)*IINTFP*IUHF
      I0R(1)=I0T(2)+IRPDPD(1,10)*IINTFP
      I0R(2)=I0R(1)+IRPDPD(IRREPX,9)*IINTFP*IUHF
      ISTART=I0R(2)+IRPDPD(IRREPX,10)*IINTFP
C
      CALL GETLST(ICORE(I0R(1)),1,1,1,3,490)
      CALL GETLST(ICORE(I0T(1)),1,1,1,1,90)
      IF(IUHF.EQ.1) THEN
       CALL GETLST(ICORE(I0R(2)),1,1,1,4,490)
       CALL GETLST(ICORE(I0T(2)),1,1,1,2,90)
      ENDIF
C
C LOOP OVER SPIN CASES
C
      DO 1000 ISPIN=3-2*IUHF,3
C
       LISTT=43+ISPIN
       LISTR=460+ISPIN
       LISTL=443+ISPIN
       LISTG=113+ISPIN
C
       I000=ISTART
C
       DO 100 IRREPR=1,NIRREP
C
        IRREPL=DIRPRD(IRREPR,IRREPX)
C
        NUMSYG=IRPDPD(IRREPR,ISYTYP(2,LISTG))
        DISSYG=IRPDPD(IRREPR,ISYTYP(1,LISTG))
C
        CALL ZERO(ICORE(I000),NUMSYG*DISSYG)
c        CALL GETLST(ICORE(I000),1,NUMSYG,1,IRREPR,LISTG)
C
C FIRST TERM: 4* 1/16 (RTAU(IJ,EF)) L(MN,EF) TAU(MN,AB)
C
        NUMSYT=IRPDPD(IRREPR,ISYTYP(2,LISTT))
        DISSYT=IRPDPD(IRREPR,ISYTYP(1,LISTT))
        NUMSYR=IRPDPD(IRREPR,ISYTYP(2,LISTR))
        DISSYR=IRPDPD(IRREPL,ISYTYP(1,LISTR))
        NUMSYL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
        DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
C
C ALLOCATE MEMORY
C
        I010=I000+IINTFP*NUMSYG*DISSYG
        I020=I010+IINTFP*MAX(NUMSYR*DISSYR,NUMSYT*DISSYT)
        I030=I020+IINTFP*NUMSYL*DISSYL
        IEND=I030+IINTFP*NUMSYL*NUMSYR
C
        IF(IEND.GT.MAXCOR) CALL INSMEM('TPDIJAB4',IEND,MAXCOR)
C
C FIRST CONTRACTION
C
C   GET R AND L AMPLITUDES
C
        IF (CC2) THEN
        CALL ZERO(ICORE(I010),NUMSYR*DISSYR)
        ELSE
        CALL GETLST(ICORE(I010),1,NUMSYR,1,IRREPR,LISTR)
        ENDIF 
C
        IF(ISPIN.LE.2) THEN
         CALL DTAU(IRREPL,IRREPR,1,IRREPX,ICORE(I010),
     &             ICORE(I0T(ISPIN)),ICORE(I0T(ISPIN)),
     &             ICORE(I0R(ISPIN)),ICORE(I0R(ISPIN)),
     &             ISPIN,ONE)
        ELSE
         CALL DTAU(IRREPL,IRREPR,1,IRREPX,ICORE(I010),
     &             ICORE(I0T(1)),ICORE(I0T(2)),
     &             ICORE(I0R(1)),ICORE(I0R(2)),
     &             3,ONE)
        ENDIF
C

        CALL GETLST(ICORE(I020),1,NUMSYL,1,IRREPR,LISTL)
C
c        CALL XGEMM('T','N',NUMSYL,NUMSYR,DISSYR,ONE,ICORE(I010),DISSYL,
c     &             ICORE(I020),DISSYR,AZERO,ICORE(I030),NUMSYL)
        CALL XGEMM('T','N',NUMSYL,NUMSYR,DISSYR,ONE,ICORE(I020),DISSYL,
     &             ICORE(I010),DISSYR,AZERO,ICORE(I030),NUMSYL)
C
C   GET T AMPLITUDES
C
        IF (CC2) THEN
        CALL ZERO(ICORE(I010),NUMSYT*DISSYT)
        ELSE
        CALL GETLST(ICORE(I010),1,NUMSYT,1,IRREPR,LISTT)
        ENDIF
C
C FORM TAU AMPLITUDES
C
        IF(ISPIN.LT.3) THEN
C
          CALL FTAU(ICORE(I010),ICORE(I0T(ISPIN)),ICORE(I0T(ISPIN)),
     &              DISSYT,NUMSYT,POP(1,ISPIN),POP(1,ISPIN),
     &              VRT(1,ISPIN),VRT(1,ISPIN),IRREPR,ISPIN,ONE)
C
        ELSE
C
          CALL FTAU(ICORE(I010),ICORE(I0T(1)),ICORE(I0T(2)),
     &              DISSYT,NUMSYT,POP(1,1),POP(1,2),
     &              VRT(1,1),VRT(1,2),IRREPR,ISPIN,ONE)
C
        ENDIF
C
        CALL XGEMM('N','N',DISSYT,NUMSYR,NUMSYL,ONE,ICORE(I010),DISSYT,
     &             ICORE(I030),NUMSYL,ONE,ICORE(I000),DISSYT)
C
C SECOND TERM: 4* 1/16 (RTAU(MN,AB) L(MN,EF) TAU(IJ,EF)
C
C
        NUMSYT=IRPDPD(IRREPR,ISYTYP(2,LISTT))
        DISSYT=IRPDPD(IRREPR,ISYTYP(1,LISTT))
        NUMSYR=IRPDPD(IRREPL,ISYTYP(2,LISTR))
        DISSYR=IRPDPD(IRREPR,ISYTYP(1,LISTR))
        NUMSYL=IRPDPD(IRREPL,ISYTYP(2,LISTL))
        DISSYL=IRPDPD(IRREPR,ISYTYP(1,LISTL))
C
C ALLOCATE MEMORY
C
        I020=I010+IINTFP*MAX(NUMSYR*DISSYR,NUMSYT*DISSYT)
        I030=I020+IINTFP*NUMSYL*DISSYL
        IEND=I030+IINTFP*NUMSYL*NUMSYT
C
        IF(IEND.GT.MAXCOR) CALL INSMEM('TPDIJAB4',IEND,MAXCOR)
C
C FIRST CONTRACTION
C
C   GET T AND L AMPLITUDES
C
        IF (CC2) THEN
        CALL ZERO(ICORE(I010),NUMSYT*DISSYT)
        ELSE
        CALL GETLST(ICORE(I010),1,NUMSYT,1,IRREPR,LISTT)
        ENDIF 
C
C FORM TAU AMPLITUDES
C
        IF(ISPIN.LT.3) THEN
C
          CALL FTAU(ICORE(I010),ICORE(I0T(ISPIN)),ICORE(I0T(ISPIN)),
     &              DISSYT,NUMSYT,POP(1,ISPIN),POP(1,ISPIN),
     &              VRT(1,ISPIN),VRT(1,ISPIN),IRREPR,ISPIN,ONE)
C
        ELSE
C
          CALL FTAU(ICORE(I010),ICORE(I0T(1)),ICORE(I0T(2)),
     &              DISSYT,NUMSYT,POP(1,1),POP(1,2),
     &              VRT(1,1),VRT(1,2),IRREPR,ISPIN,ONE)
C
        ENDIF
C
        CALL GETLST(ICORE(I020),1,NUMSYL,1,IRREPL,LISTL)
C
        CALL XGEMM('T','N',NUMSYL,NUMSYT,DISSYT,ONE,ICORE(I020),DISSYL,
     &             ICORE(I010),DISSYT,AZERO,ICORE(I030),NUMSYL)
C
C   GET R AMPLITUDES
C
        IF (CC2) THEN
        CALL ZERO(ICORE(I010),NUMSYR*DISSYR)
        ELSE
        CALL GETLST(ICORE(I010),1,NUMSYR,1,IRREPL,LISTR)
        ENDIF 
C
C FORM DTAU AMPLITUDES
C
        IF(ISPIN.LE.2) THEN
         CALL DTAU(IRREPR,IRREPL,1,IRREPX,ICORE(I010),
     &             ICORE(I0T(ISPIN)),ICORE(I0T(ISPIN)),
     &             ICORE(I0R(ISPIN)),ICORE(I0R(ISPIN)),
     &             ISPIN,ONE)
        ELSE
         CALL DTAU(IRREPR,IRREPL,1,IRREPX,ICORE(I010),
     &             ICORE(I0T(1)),ICORE(I0T(2)),
     &             ICORE(I0R(1)),ICORE(I0R(2)),
     &             3,ONE)
        ENDIF
C
        CALL XGEMM('N','N',DISSYR,NUMSYT,NUMSYL,ONE,ICORE(I010),DISSYR,
     &             ICORE(I030),NUMSYL,ONE,ICORE(I000),DISSYR)
C
        I010=I000+IINTFP*NUMSYG*DISSYG
        IEND=I010+IINTFP*NUMSYG*DISSYG
        IF(IEND.GE.MAXCOR) CALL INSMEM('TPDIJAB4',IEND,MAXCOR)
        CALL GETLST(ICORE(I010),1,NUMSYG,1,IRREPR,LISTG)
        CALL SAXPY(NUMSYG*DISSYG,ONE,ICORE(I010),1,ICORE(I000),1)
        CALL PUTLST(ICORE(I000),1,NUMSYG,1,IRREPR,LISTG)
C
100    CONTINUE
C
1000  CONTINUE

      IF(.NOT.CC2) THEN
       TWO=0.D0
      ELSE
       TWO=2.D0
       DO ISPIN=3-2*IUHF,3
       LISTG=113+ISPIN
       NSIZE=ISYMSZ(ISYTYP(1,LISTG),ISYTYP(2,LISTG))
       CALL GETALL(ICORE,NSIZE,1,LISTG)
       CALL SSCAL(NSIZE,0.5D0,ICORE,1)
       CALL PUTALL(ICORE,NSIZE,1,LISTG)
       ENDDO 
      ENDIF 
C
C ALL DONE, RETURN
C
      TWO=0.D0
      IF(IUHF.EQ.0)
     &   CALL CHECKGAM1(ICORE,16,116,TWO,IUHF,2,VRT)
      IF(IUHF.EQ.1) THEN
       CALL CHECKGAM(ICORE,16,116,TWO)
       CALL CHECKGAM(ICORE,14,114,TWO)
       CALL CHECKGAM(ICORE,15,115,TWO)
      ENDIF
C
      RETURN
      END