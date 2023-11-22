
      SUBROUTINE FORMXIJ(AIOO,DOO,XOO,EVAL,IUHF,QRHF)
C
C   FORM FIRST X(IJ) WITH I > J, THEN CONSTRUCT THE REAL
C   I(IJ) INTERMEDIATES, GET THE Z(IJ) AS
C
C   Z(IJ) = X(IJ)/(f(II)-f(JJ)
C
C   D(IJ) = - Z(IJ)
C
C   AND AUGMENT I(IJ) BY f(I,I) D(I,J)
C
C NOTE, SO FAR ONLY IMPLEMENTED FOR HF BASED CORRELATION METHODS
C
CEND
C
C CODED 3/91 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER POP,VRT,DIRPRD
      CHARACTER*1 SPCASE
CJDW KKB stuff
      LOGICAL QRHF
CJDW END
C
      DIMENSION AIOO(1),DOO(1),XOO(1),EVAL(1)
      DIMENSION SPCASE(2)
C
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPP(255,2),DIRPRD(8,8)
CJDW KKB stuff
C     COMMON/SYM /POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYM2/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
CJDW END
C
      DATA TRESH /1.D-12/,ONE /1.D0/
      DATA SPCASE /'A','B'/
C
C  LOOP OVER SPIN CASES
C
      DO 1000 ISPIN=1,IUHF+1
C
C DETERMINE OFFSETS
C
       IOO=1+(ISPIN-1)*NFMI(1)
       IOV=1+(ISPIN-1)*NT(1)
C
       CALL ZERO(XOO(IOO),NFMI(ISPIN))
C
C ALLOCATE CORE FOR EIGENVALUES
C
CJDW KKB stuff
C       IF (QRHF) THEN
C        CALL GETREC(20,'JOBARC','RHFEVAL ',
C     &             IINTFP*(NOCCO(ISPIN)+NVRTO(ISPIN)),
C     &             EVAL)
C        write (6,9091)
C 9091   format (///4x,'==== analytic gradient with REF=QRHF,',
C     &             ' PERT_ORB=CANONICAL is incorrect right now.'///)
C        call errex 
c-------------------
c--     canonical pert_orb in qrhf gradient is not correct at present.
c--     I'll be back to this point later    ----- Mar. 95 -- KB  -----
c-------------------
C       ELSE
CJDW END
       CALL GETREC(20,'JOBARC','SCFEVAL'//SPCASE(ISPIN),
     &             IINTFP*(NOCCO(ISPIN)+NVRTO(ISPIN)),
     &             EVAL)
CJDW KKB stuff
C       ENDIF
CJDW END
C
C  X(I,J) = I(J,I) - I(I,J)  I > J
C
C THEN SYMMETRIZE I(I,J) BY I(I,J) = I(J,I) FOR ALL I < J
C
       IOFF=0
       NOFF=0
       DO 100 IRREP=1,NIRREP
C
        NOCC=POP(IRREP,ISPIN)
C
        DO 20 I=1,NOCC
         DO 10 J=1,I-1
C
         IND1=IOO+IOFF+I+(J-1)*NOCC-1
         IND2=IOO+IOFF+J+(I-1)*NOCC-1
C
         EVDIFF=EVAL(J+NOFF)-EVAL(I+NOFF)
C
         XOOI=AIOO(IND2)-AIOO(IND1)
C
C Marcel observed that the DROPMO gradient of s-triazine fails as
C a result of the instability introduced by small EVDIFF. It is 
C normally assumed that when EVDIFF small XOOI is also small and
C as a result the only cheked performed was for the size of XOII.
C But in this case that assumption was not accuarte and the size
C of EVDIFF must also be checked. 01/2011, Ajith Perera.
C 

         IF (ABS(EVDIFF) .GT. 1.0D-9) THEN
            IF(ABS(XOOI).GT.TRESH) XOO(IND1)=XOOI/EVDIFF
         ENDIF
C
10       CONTINUE
C
         DO 12 J=1,I-1
C
         IND1=IOO+IOFF+I+(J-1)*NOCC-1
         IND2=IOO+IOFF+J+(I-1)*NOCC-1
C
         AIOO(IND2)=AIOO(IND2)-EVAL(J+NOFF)*XOO(IND1)
C
12       CONTINUE 
C
         DO 15 J=1,I-1
C
         IND1=IOO+IOFF+I+(J-1)*NOCC-1
         IND2=IOO+IOFF+J+(I-1)*NOCC-1
C
         XOO(IND2)=XOO(IND1)
         AIOO(IND1)=AIOO(IND2)
C
15       CONTINUE
C
20      CONTINUE
C
        NOFF=NOFF+NOCC
        IOFF=IOFF+NOCC*NOCC
C
100    CONTINUE
C
       CALL SAXPY(NFMI(ISPIN),ONE,XOO(IOO),1,DOO(IOO),1)
C
1000  CONTINUE
C
      RETURN
      END
