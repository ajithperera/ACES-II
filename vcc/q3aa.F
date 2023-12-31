      SUBROUTINE Q3AA(W,T,T1,Z,T1A,TAU,ISPIN,POP,VRT,NVRTSQ,
     &                NVRT2SQ,DISSYW,DISSYT,NUMSYW,NUMSYT,
     &                IRREP,TMP,FACT,UCC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TAU,UCC,ADC2
      INTEGER DISSYT, DISSYW, DIRPRD,POP,VRT
      DIMENSION W(NUMSYW,NVRT2SQ),T(NUMSYT,NVRT2SQ)
      DIMENSION T1(DISSYT,NUMSYT),Z(NVRTSQ)
      DIMENSION TMP(1)
      DIMENSION POP(8),VRT(8),T1A(1)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
C
      DATA ONEM,ONE,HALF,TWO /-1.0D0,1.D0,0.5D0,2.D0/
C
      IND(I,J)=((I-2)*(I-1))/2+J
C
C ADC2 check; CALC=MBPT(2),EXCITE=EOMEE,EOMREF=ADC2
C
      ADC2 = .FALSE.
      ADC2 =  (IFLAGS(2) .EQ. 1 .AND. IFLAGS(87) .EQ. 3 .AND.
     &         IFLAGS2(117) .EQ. 10)

       FACTM=-FACT
C      PICK UP FIRST THE RELEVANT T2 AND W PIECES
C
      CALL GETLST(T1,1,NUMSYW,2,IRREP,ISPIN+13)
      CALL TRANSP(T1,W,NUMSYW,DISSYW)
      CALL GETLST(T1,1,NUMSYT,1,IRREP,ISPIN+43)   

      IF(TAU) THEN
       CALL FTAU(T1,T1A,T1A,DISSYT,NUMSYT,POP,POP,VRT,VRT,
     &           IRREP,ISPIN,HALF)
      ENDIF
      CALL TRANSP(T1,T,NUMSYT,DISSYT)
C
      CALL SYMEXP(IRREP,VRT,NUMSYW,W)
      CALL SYMEXP(IRREP,VRT,NUMSYT,T)
C 
         IOFF=1
         JOFF=1
         DO 90 IRREPJ=1,NIRREP
C          
C        GET NUMBER OF VIRTUAL ORBITALS FOR IRREPJ     
C
           NVRTJ=VRT(IRREPJ)
C
C        IF ZERO, NOTHING TO COMPUTE
C
           IF(NVRTJ.EQ.0) GO TO 90
C
C        DETERMINE IRREPI WHOSE DIRECT PRODUCT WITH IRREPJ GIVES IRREP
C
           IRREPI=DIRPRD(IRREP,IRREPJ)
C
C        GET NUMBER  OF VIRTUAL ORBITALS FOR IRREPI
C
           NVRTI=VRT(IRREPI)
C
C        IF ZERO, NOTHING TO COMPUTE
C
           IF(NVRTI.EQ.0) GO TO 80
C
            CALL XGEMM('T','N',NVRTJ,NVRTJ,NVRTI*NUMSYW,FACTM,
     &                  W(1,JOFF),NVRTI*NUMSYW,T(1,JOFF),
     &                  NVRTI*NUMSYW,ONE,Z(IOFF),NVRTJ)

C
C The ADC(2) variants needs symmetrization of this contribution. So,
C check whether we are doing EOM-MBPT(2) with the EOM-REF=ADC(2).
C
C            IF (ADC2) CALL SYMMET2(Z(IOFF), NVRTJ)

            IF(UCC) CALL DT(ONE,Z(IOFF),NVRTJ,NVRTJ,'ADD',Z(IOFF))
C
            JOFF=JOFF+NVRTI*NVRTJ
80         CONTINUE
           IOFF=IOFF+NVRTJ*NVRTJ
C
90        CONTINUE
C
      RETURN
      END
