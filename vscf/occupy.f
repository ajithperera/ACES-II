      SUBROUTINE OCCUPY(NIRREP,NBASIR,NBASTOT,EVAL,ISCR,IUHF)
C
C DETERMINES THE OCCUPANCY FOR AN SCF DETERMINANT BASED ON
C  AN EXISTING SET OF EIGENVALUES
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EVAL(NBASTOT,1+IUHF),ISCR(NBASTOT,2),NSPIN(2)
      DIMENSION NBASIR(NIRREP),NOCOLD(8,2)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /POPUL/ NOCC(8,2)
C
      DATA I16 /16/
C
      CALL ICOPY(I16,NOCC,1,NOCOLD,1)
      IONE=1
      CALL GETREC(20,'JOBARC','NMPROTON',IONE,NPROTON)
      ICHARG=IFLAGS(28)
      IMULT=IFLAGS(29)
C
C DETERMINE ALPHA AND BETA OCCUPANCIES
C
      NUMEL=NPROTON-ICHARG
      IALPEX=IMULT-1
      NRHS=NUMEL-IALPEX
      IF(MOD(NRHS,2).NE.0)THEN
       WRITE(6,100)
100    FORMAT(T3,'@OCCUPY-F, Specified charge and multiplicity are ',
     &           'impossible.   Try again.')
       CALL ERREX
      ENDIF
      NSPIN(2)=NRHS/2
      NSPIN(1)=NSPIN(2)+IALPEX
C
      DO 30 ISPIN=1,1+IUHF
       IBASIS=0
       DO 10 I=1,NIRREP
        DO 20 J=1,NBASIR(I)
         IBASIS=IBASIS+1
         ISCR(IBASIS,1)=IBASIS
         ISCR(IBASIS,2)=I
20      CONTINUE
10     CONTINUE
       CALL IZERO(NOCC(1,ISPIN),8)
       CALL PIKSR2(NBASTOT,EVAL(1,ISPIN),ISCR)
       DO 40 IOCC=1,NSPIN(ISPIN)
        IPOS=ISCR(IOCC,1)
        IRREP=ISCR(IPOS,2)
        NOCC(IRREP,ISPIN)=NOCC(IRREP,ISPIN)+1
40     CONTINUE
30    CONTINUE
C
      IF(IUHF.EQ.0) THEN
        DO 110 I=1,8
          NOCC(I,2)=NOCC(I,1)
  110   CONTINUE
      ENDIF
C
      RETURN
      END