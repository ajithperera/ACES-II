      SUBROUTINE G5AA(T2,G,T1,T1A,GCA,FACT,LAMBDA,CCSD,TRIP,
     &                ISPIN,POP,VRT,NT,DISSYT,NUMSYT,DISSYG,
     &                NUMSYG,LISTT,LISTL,LISTG,IRREP,TMP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL LAMBDA,CCSD,TRIP
      INTEGER DIRPRD,DISSYT,DISSYG,POP,VRT
      DIMENSION T2(DISSYT,1),G(DISSYG,1),T1(NT,2),TMP(1) 
      DIMENSION POP(8),VRT(8),T1A(1),GCA(1)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
      DATA AZERO,ONE /0.D0,1.D0/
C
      IF(CCSD.OR.TRIP) THEN
       CALL GETLST(G,1,NUMSYG,1,IRREP,LISTG)
       IF(MIN(NUMSYT,DISSYT).NE.0) THEN
       CALL SYMTR1(IRREP,VRT,POP,DISSYG,G,TMP,TMP(1+DISSYG),
     &             TMP(1+2*DISSYG))
       ENDIF
      ELSE
       CALL ZERO(G,NUMSYG*DISSYG)
      ENDIF
C
      IF(MIN(NUMSYT,DISSYT).NE.0) THEN
      IPASS=1
C
C     IPASS=1 CORRESPONDS TO LAMBDA(M,C) T(MI,AB)
C     IPASS=2 CORRESPONDS TO T(M,C) LAMBDA(MI,AB)
C
1     CONTINUE
C
C  SET PARAMETERS FOR LOOP
C
      IF(IPASS.EQ.1) THEN
       LIST1=LISTT
      ELSE
       LIST1=LISTL
      ENDIF
C
C    GET FIRST T2 AMPLITUDES FROM DISK
C    AND EXPAND THE OCCUPIED-OCCUPIED BLOCK
C
      CALL GETLST(T2,1,NUMSYT,1,IRREP,LIST1)
C
C  FOR CCSD FORM TAU AMPLITUDES IN FIRST PASS
C
      IF(IPASS.EQ.1.AND.CCSD) THEN
       CALL FTAU(T2,T1A,T1A,DISSYT,NUMSYT,POP,POP,VRT,VRT,
     &           IRREP,ISPIN,ONE)
      ENDIF
C
      CALL SYMEXP(IRREP,POP,DISSYT,T2)
C
C     NOW PERFORM MULTIPLICATION
C
C  IOFF OFFSET IN T1
C  JOFFT2 OFFSET IN OCC. OCC. BLOCK OF T2
C  JOFFG OFFSET IN OCC VIRT BLOCK OF G
C
      IOFF=1
      JOFFT2=1
      JOFFG=1
C
      DO 100 IRREPJ=1,NIRREP
C
C  OCCUPIED AND VIRTAUl ORBITALS FOR MULTIPLICATION
C
      NOCCJ=POP(IRREPJ)
      NVRTJ=VRT(IRREPJ) 
      IRREPI=DIRPRD(IRREP,IRREPJ)
      NOCCI=POP(IRREPI)
C
C  IF ANY OF THE POPULATION IS ZERO, SKIP LOOP
C
      IF(MIN(NOCCI,NVRTJ,NOCCJ).NE.0) THEN
C
       CALL XGEMM('N','T',DISSYT*NOCCI,NVRTJ,NOCCJ,FACT,T2(1,JOFFT2),
     &            DISSYT*NOCCI,T1(IOFF,IPASS),NVRTJ,ONE,G(1,JOFFG),
     &            DISSYT*NOCCI)
      ENDIF
C
C   UPDATE POINTERS
C
      IOFF=IOFF+NVRTJ*NOCCJ
      JOFFT2=JOFFT2+NOCCI*NOCCJ
      JOFFG=JOFFG+NOCCI*NVRTJ
C
100   CONTINUE
C
      IPASS=IPASS+1
C
C  ONLY FOR QCISD AND CCSD TWO PASSES ARE REQUIRED
C
      IF(LAMBDA.AND.IPASS.EQ.2) GO TO 1
C
C  TRANSPOSE INDICES SINCE WE HAVE CALCULATES -G(AB,IC) = G(AB,CI)
C  WITH AN ORDERING OF  A,B,I,C
C
      CALL SYMTR1(IRREP,POP,VRT,DISSYG,G,TMP,TMP(1+DISSYG),
     &            TMP(1+2*DISSYG))
C
      ENDIF
C
C  FOR CCSD ADD THE -1/2 P(AB) G(CA)*T(IB) CONTRIBUTION HERE
C
      IF(CCSD) THEN
C
C  GET GCA FROM DISK
C
       CALL GETLST(GCA,1,1,1,ISPIN,192)
C
C   CALCULATE THE CORRESPONDING CONTRIBUTION IN G5TAU
C
       CALL G5TAU(G,GCA,T1A,DISSYG,NUMSYG,VRT,POP,VRT,
     &            IRREP,ISPIN)
      ENDIF
C
C  SAVE GAMMA ON LIST
C
      CALL PUTLST(G,1,NUMSYG,2,IRREP,LISTG)
C
      RETURN
      END
