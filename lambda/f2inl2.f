      SUBROUTINE F2INL2(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE CALCULATES THE TERM
C
C  - P(IJ) SUM M L(IM,AB) F(JM)
C
C IN RHF :
C
C - SUM m L(Im,Ab) F(jm) + SUM M L(Mj,Ab) F(JM)
C
C IN UHF
C
C - SUM M L(IM,AB) F(JM) + SUM M L(MJ,AB) F(IM)
C
C - SUM m L(Im,Ab) F(jm) + SUM M L(Mj,Ab) F(IM)
C
C - SUM m L(im,ab) F(jm) + SUM m L(mj,ab) F(im)
C
C 
C FOR QCISD AND CCSD IN ADDITION THE TERMS
C
C - SUM M L(M,A) F(IM) 
C
C - SUM m L(m,a) F(im)   (UHF only)
C
C ARE CALCULATED
C
C FOR CCSD METHODS THE FOLLOWING TERM
C
C 1/2 SUM E F(E,J) T(E,M)
C
C HAS TO BE ADDED TO THE F(IM) BEFORE THE
C CONTRACTION WITH THE L2 AMPLITUDES IS PERFORMED
C
CEND
C
C CODED JG AUGUST/90
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,DISSYZ,POP,VRT
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,CC2
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWO  
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,
     &             NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON /FLAGS/IFLAGS(100)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /METH/ MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &              CC2
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
      CALL GETLST(ICORE(I0AA),1,1,1,1,91)
      IF(IUHF.EQ.0) THEN
       I0BB=I0AA
      ELSE
       I0BB=I0AA-NFBB*IINTFP
       MXCOR=MXCOR-NFBB*IINTFP
       CALL GETLST(ICORE(I0BB),1,1,1,2,91)
      ENDIF
C
C ADDITIONAL CODE FOR QCISD AND CCSD METHODS
C
       IF(QCISD.OR.CCSD) THEN
        I0TA=I0BB-NTAA*IINTFP
        I0ZA=I0TA-NTAA*IINTFP
        MXCOR=MXCOR-2*NTAA*IINTFP
        CALL GETLST(ICORE(I0TA),1,1,2,1,190)
        CALL GETLST(ICORE(I0ZA),1,1,1,3,90)
        IF(IUHF.EQ.0) THEN
         I0TB=I0TA
         I0ZB=I0ZA
        ELSE
         I0TB=I0ZA-NTBB*IINTFP
         I0ZB=I0TB-NTBB*IINTFP
         MXCOR=MXCOR-2*NTBB*IINTFP
         CALL GETLST(ICORE(I0TB),1,1,2,2,190)
         CALL GETLST(ICORE(I0ZB),1,1,1,4,90)
        ENDIF
CC
C  IF CCSD METHODS ADD TO FMI THE FOLLOWING TERM
C
C  1/2 SUM E  F(E,J) T(E,M)
C
       IF(.NOT.QCISD) THEN
C
C  ALLOCATE MEMORY BUT DO NOT UPDATE MXCOR
C
       I01=I0ZB-IINTFP*MAX(NTAA,NTBB)
       I02=I01-IINTFP*MAX(NTAA,NTBB)
C
C  GET T1 AMPLITUDES AND F(ME) INTERMEDIATES
C 
        CALL GETLST(ICORE(I01),1,1,1,1,93)
        CALL GETLST(ICORE(I02),1,1,1,1,90)
C
        IOFFMI=I0AA
        IOFFEM=I01
        IOFT=I02
        DO 20 IRREP=1,NIRREP
         NOCC=POP(IRREP,1)
         NVRT=VRT(IRREP,1)
         CALL XGEMM('T','N',NOCC,NOCC,NVRT,HALF,ICORE(IOFFEM),NVRT,
     &              ICORE(IOFT),NVRT,ONE,ICORE(IOFFMI),NOCC)
         IOFFMI=IOFFMI+NOCC*NOCC*IINTFP
         IOFFEM=IOFFEM+NOCC*NVRT*IINTFP
         IOFT=IOFT+NOCC*NVRT*IINTFP
20      CONTINUE
C
        IF(IUHF.EQ.1) THEN
C
C  GET T1 AMPLITUDES AND F(ME) INTERMEDIATES
C
         CALL GETLST(ICORE(I01),1,1,1,2,93)
         CALL GETLST(ICORE(I02),1,1,1,2,90)
C
         IOFFMI=I0BB
         IOFFEM=I01
         IOFT=I02
         DO 21 IRREP=1,NIRREP
          NOCC=POP(IRREP,2)
          NVRT=VRT(IRREP,2)
          CALL XGEMM('T','N',NOCC,NOCC,NVRT,HALF,ICORE(IOFFEM),NVRT,
     &               ICORE(IOFT),NVRT,ONE,ICORE(IOFFMI),NOCC)
          IOFFMI=IOFFMI+NOCC*NOCC*IINTFP
          IOFFEM=IOFFEM+NOCC*NVRT*IINTFP
          IOFT=IOFT+NOCC*NVRT*IINTFP
21       CONTINUE
        ENDIF
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
      ENDIF
C
C CC2 and f(m,e) != 0, we need to add +1/2f(m,e)*t(i,e) contribution to
C F(m,i) intermediate.
C
      IF (CC2) THEN

         IF (IFLAGS(38) .EQ. 0) RETURN

            CALL GETLST(ICORE(I01),1,1,1,3,93)
            CALL GETLST(ICORE(I02),1,1,1,1,90)
            CALL GETLST(ICORE(I0AA),1,1,1,3,91)
            
            IOFFMI=I0AA
            IOFFEM=I01
            IOFT=I02

            DO IRREP=1,NIRREP
               NOCC=POP(IRREP,1)
               NVRT=VRT(IRREP,1)
               CALL XGEMM('T','N',NOCC,NOCC,NVRT,ONE,ICORE(IOFFEM),
     &                     NVRT,ICORE(IOFT),NVRT,ONE,ICORE(IOFFMI),
     &                     NOCC)
               IOFFMI=IOFFMI+NOCC*NOCC*IINTFP
               IOFFEM=IOFFEM+NOCC*NVRT*IINTFP
               IOFT=IOFT+NOCC*NVRT*IINTFP 
            ENDDO

            IF (IUHF.EQ.1) THEN
               CALL GETLST(ICORE(I01),1,1,1,4,93)
               CALL GETLST(ICORE(I02),1,1,1,2,90)
               CALL GETLST(ICORE(I0BB),1,1,1,4,91)

               IOFFMI=I0BB
               IOFFEM=I01
               IOFT=I02

               DO IRREP=1,NIRREP
                  NOCC=POP(IRREP,2)
                  NVRT=VRT(IRREP,2)
                  CALL XGEMM('T','N',NOCC,NOCC,NVRT,ONE,ICORE(IOFFEM),
     &                        NVRT,ICORE(IOFT),NVRT,ONE,ICORE(IOFFMI),
     &                        NOCC)
                  IOFFMI=IOFFMI+NOCC*NOCC*IINTFP
                  IOFFEM=IOFFEM+NOCC*NVRT*IINTFP
                  IOFT=IOFT+NOCC*NVRT*IINTFP
               ENDDO

           ENDIF

      ENDIF
C
C     AA AND BB SPIN CASES
C
      IF(IUHF.EQ.1) THEN
C
C
       DO 100 ISPIN=1,2
C
        IF(ISPIN.EQ.1) THEN
         I000=I0AA
        ELSE
         I000=I0BB
        ENDIF
        LISTT=ISPIN+143
        LISTZ=ISPIN+60
C
        DO 50 IRREP=1,NIRREP 
C
        NOCCSQ=0
        DO 45 IRREPJ=1,NIRREP
         IRREPI=DIRPRD(IRREP,IRREPJ)
         NOCCSQ=NOCCSQ+POP(IRREPJ,ISPIN)*POP(IRREPI,ISPIN)
45      CONTINUE
C
         DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT)) 
         DISSYZ=IRPDPD(IRREP,ISYTYP(1,LISTZ))
         NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
         NUMSYZ=IRPDPD(IRREP,ISYTYP(2,LISTZ))
         I001=1
         I002=I001+IINTFP*NOCCSQ*DISSYT
         I003=I002+IINTFP*NOCCSQ*DISSYZ
         IF(MIN(NUMSYT,NUMSYZ,DISSYT,DISSYZ).NE.0) THEN
          I004=I003+IINTFP*MAX(DISSYT,DISSYZ)
          IF(I004.LT.MXCOR) THEN
C
C
C    IN CORE VERSION
C
          CALL F2L2AA(ICORE(I001),ICORE(I002),ICORE(I000),
     &                POP(1,ISPIN),NOCCSQ,DISSYT,DISSYZ,NUMSYT,
     &                NUMSYZ,NFAA,LISTT,LISTZ,IRREP,ICORE(I003))
          ELSE
           STOP 'F2L2AA'
          ENDIF
         ENDIF
50      CONTINUE
100    CONTINUE
      ENDIF
  
C
C      AB SPIN CASE
C
      LISTT=146
      LISTZ=63
C
C   LOOP OVER IRREPS
C
      DO 200 IRREP=1,NIRREP
C
C   RETRIEVE T2 AMPLITUDES AND CALCULATE Z-AMPLITUDES
C
       DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
       DISSYZ=IRPDPD(IRREP,ISYTYP(1,LISTZ))
       NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
       NUMSYZ=IRPDPD(IRREP,ISYTYP(2,LISTZ))
       I001=1
       I002=I001+IINTFP*NUMSYT*DISSYT
       I003=I002+IINTFP*NUMSYZ*DISSYZ
       IF(MIN(NUMSYT,NUMSYZ,DISSYT,DISSYZ).NE.0) THEN
        I004=I003+IINTFP*MAX(DISSYT,DISSYZ,NUMSYT,NUMSYZ)*3
        IF(I004.LT.MXCOR) THEN
C
C       IN CORE ALGORITHM
C
         CALL F2L2AB(ICORE(I001),ICORE(I002),ICORE(I0AA),ICORE(I0BB),
     &               POP(1,1),POP(1,2),VRT(1,1),VRT(1,2),DISSYT,
     &               DISSYZ,NUMSYT,NUMSYZ,NFAA,NFBB,LISTT,LISTZ,
     &               IRREP,IUHF,ICORE(I003))
        ELSE
C
C       OUT CORE ALGORITHM
C
         STOP 'F2L2AB'
        ENDIF
       ENDIF
200   CONTINUE
      RETURN
      END
