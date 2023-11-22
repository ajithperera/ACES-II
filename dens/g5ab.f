      SUBROUTINE G5AB(T2,G,T1,T1A,T1B,GCA,FACT,LAMBDA,CCSD,TRIP,
     &                ISPIN,POP1,POP2,VRT1,VRT2,NT,DISSYT,NUMSYT,
     &                DISSYG,NUMSYG,LISTT,LISTL,LISTG,IRREP,TMP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL LAMBDA,CCSD,TRIP
      INTEGER DIRPRD,DISSYT,DISSYG,POP1,POP2,VRT1,VRT2
      DIMENSION T2(DISSYT,1),G(DISSYG,1),T1(NT,2),T1A(1),T1B(1) 
      DIMENSION POP1(8),POP2(8),VRT1(8),VRT2(8),TMP(1),GCA(1)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
      DATA AZERO,ONE /0.D0,1.D0/
C
      IF(CCSD.OR.TRIP) THEN
       CALL GETLST(G,1,NUMSYG,1,IRREP,LISTG)
       IF(ISPIN.EQ.1.AND.(MIN(NUMSYT,DISSYT).NE.0)) THEN
        CALL SYMTR1(IRREP,VRT1,POP2,DISSYG,G,TMP,TMP(1+DISSYG),
     &              TMP(1+2*DISSYG))
       ENDIF
      ELSE
       CALL ZERO(G,NUMSYG*DISSYG)
      ENDIF
C
      IF(MIN(NUMSYT,DISSYT).NE.0) THEN

      IPASS=1
C
C   IPASS =1 CORRESPONDS TO LAMBDA(M,C) T(MI,AB)
C   IPASS =2 CORRESPONDS TO T(M,C) LAMBDA(MI,AB)
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
C
C    IF ISPIN=1 TRANSPOSE FROM A,b ; I,j to A,b ; j,I
C    IF ISPIN=2 DO NOTHING WITH T2
C
      CALL GETLST(T2,1,NUMSYT,1,IRREP,LIST1)
C
C    FOR CCSD FORM TAU AMPLITUDES IN THE FIRST PASS
C
      IF(IPASS.EQ.1.AND.CCSD) THEN
       IF(ISPIN.EQ.1) THEN
        CALL FTAU(T2,T1A,T1B,DISSYT,NUMSYT,POP1,POP2,VRT1,VRT2,
     &            IRREP,3,ONE)
       ELSE
        CALL FTAU(T2,T1A,T1B,DISSYT,NUMSYT,POP2,POP1,VRT2,VRT1,
     &            IRREP,3,ONE)
       ENDIF
      ENDIF
C
      IF(ISPIN.EQ.1) THEN
       CALL SYMTR1(IRREP,POP1,POP2,DISSYT,T2,TMP,TMP(1+DISSYT),
     &             TMP(1+2*DISSYT))
      ENDIF
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
      NOCCJ=POP1(IRREPJ)
      NVRTJ=VRT1(IRREPJ) 
      IRREPI=DIRPRD(IRREP,IRREPJ)
      NOCCI=POP2(IRREPI)
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
      IF(LAMBDA.AND.IPASS.EQ.2) GO TO 1
C
C  IF ISPIN=1 THEN
C  TRANSPOSE INDICES SINCE WE HAVE CALCULATES -G(AB,IC) = G(AB,CI)
C  WITH AN ORDERING OF  A,B,I,C
C
C  FOR ISPIN=2 DO NOTHING
C
      IF(ISPIN.EQ.1) THEN
       CALL SYMTR1(IRREP,POP2,VRT1,DISSYG,G,TMP,TMP(1+DISSYG),
     &             TMP(1+2*DISSYG))
      ENDIF
      ENDIF
C
C  FOR CCSD ADD HERE THE -1/2 P(AB) G(CA)*T(IB) CONTRIBUTION
C
      IF(CCSD) THEN
C
C  GET GCA FROM DISK
C
       CALL GETLST(GCA,1,1,1,ISPIN,192)
C
C  CALL G5TAU IN ORDER TO CALCULATE THE CONTRIBUTION
C
       IF(ISPIN.EQ.1) THEN
        CALL G5TAU(G,GCA,T1B,DISSYG,NUMSYG,VRT1,POP2,
     &             VRT2,IRREP,3)
       ELSE
        CALL G5TAU(G,GCA,T1A,DISSYG,NUMSYG,VRT1,POP2,
     &             VRT2,IRREP,4)
       ENDIF
      ENDIF
C
C  SAVE GAMMA ON LIST
C
      CALL PUTLST(G,1,NUMSYG,2,IRREP,LISTG)
C
      RETURN
      END
