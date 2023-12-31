      SUBROUTINE RCL_Q3AB(W,T,T1,Z,T1A,T1B,ISPIN,POP1,POP2,
     &                    VRT1,VRT2,NVRTSQ,DISSYW,DISSYT,NUMSYW,
     &                    NUMSYT,IRREP,TMP,IUHF,FACT,IMULT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      INTEGER DISSYT,DISSYW,DIRPRD,POP1,POP2,VRT1,VRT2
      DIMENSION W(NUMSYW,DISSYW),T(NUMSYT,DISSYT),Z(NVRTSQ)
      DIMENSION T1(DISSYW,NUMSYW),T1A(1),T1B(1)
      DIMENSION TMP(1)
      DIMENSION POP1(8),POP2(8),VRT1(8),VRT2(8) 
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
C
      DATA ONE,ONEM,TWO,HALF /1.0D0,-1.0D0,2.0D0,0.5D0/
      FACTM=-FACT
C
C  SPIN ADAPTED CODE FOR RHF
C
       IF(IUHF.EQ.0) THEN
          
        IF (RCCD) THEN
           IF (IMULT .EQ. 1) THEN 
               CALL GETLST(T1,1,NUMSYW,2,IRREP,197)
               CALL TRANSP(T1,W,NUMSYW,DISSYW)
           ELSE IF (IMULT .EQ. 2) THEN
              CALL GETLST(T1,1,NUMSYW,2,IRREP,198)
              CALL TRANSP(T1,W,NUMSYW,DISSYW)
           ENDIF 
        ELSE IF (DRCCD) THEN

C If antisymmetric  <ab||ij> is used 

            CALL GETLST(T1,1,NUMSYW,2,IRREP,16)
            CALL TRANSP(T1,W,NUMSYW,DISSYW)
            CALL SPINAD1(IRREP,VRT1,NUMSYW,W,TMP,TMP(1+NUMSYW))

C If Coulomb only integrals are used

CSSS       CALL DSCAL(NUMSYW*DISSYW,TWO,W,1)
        ENDIF 
       ELSE
C
C      PICK UP FIRST THE RELEVANT T2 AND W PIECES
C
          CALL GETLST(T1,1,NUMSYW,2,IRREP,16)
          CALL TRANSP(T1,W,NUMSYW,DISSYW)
       ENDIF 

C
       IF  (IMULT .EQ. 2 .AND. IUHF .EQ. 0) THEN
          CALL GETLST(T1,1,NUMSYT,1,IRREP,44)
          CALL TRANSP(T1,T,NUMSYT,DISSYT)
       ELSE 
          CALL GETLST(T1,1,NUMSYT,1,IRREP,46)
          CALL TRANSP(T1,T,NUMSYT,DISSYT)
       ENDIF 

#ifdef _DEBUG_LVL0
      call checksum("@q3b-W:",W,NUMSYW*DISSYW)
      call checksum("@q3b-T:",T,NUMSYT*DISSYT)
#endif
C
      IF(IUHF.EQ.1.AND.ISPIN.EQ.1) THEN
        CALL SYMTR1(IRREP,VRT1,VRT2,NUMSYW,W,TMP,TMP(1+NUMSYW),
     &              TMP(1+2*NUMSYW))
        CALL SYMTR1(IRREP,VRT1,VRT2,NUMSYT,T,TMP,TMP(1+NUMSYT),
     &              TMP(1+2*NUMSYT))
      ENDIF
C
         IOFF=1
         JOFF=1
         DO 90 IRREPI=1,NIRREP
C          
C        GET NUMBER OF VIRTUAL ORBITALS FOR IRREPJ
C
           NVRTI=VRT1(IRREPI)
C
C        IF ZERO, NOTHING TO COMPUTE
C
           IF(NVRTI.EQ.0) GO TO 90
C
C        DETERMINE IRREPI WHOSE DIRECT PRODUCT WITH IRREPJ GIVES IRREP
C
           IRREPJ=DIRPRD(IRREP,IRREPI)
C
C        GET NUMBER OF VIRTUAL ORBITALS FOR KRREP
C
           NVRTJ=VRT2(IRREPJ)
C
C        IF ZERO, NOTHING TO COMPUTE
C
           IF(NVRTJ.EQ.0) GO TO 80
C
            CALL XGEMM('T','N',NVRTI,NVRTI,NUMSYW*NVRTJ,FACTM,
     &                  W(1,JOFF),NVRTJ*NUMSYW,T(1,JOFF),
     &                  NVRTJ*NUMSYT,ONE,Z(IOFF),NVRTI)

            JOFF=JOFF+NVRTJ*NVRTI
C
80         CONTINUE
C
C          UPDATE THE OFFSET FOR THE INTERMEDIATE Z
C
           IOFF=IOFF+NVRTI**2
90        CONTINUE
C
      RETURN
      END
