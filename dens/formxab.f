      SUBROUTINE FORMXAB(AIVV,DVV,XVV,EVAL,IUHF,QRHF)
C
C   FORM FIRST X(AB) WITH A > B, THEN CONSTRUCT THE REAL
C   I(AB) INTERMEDIATES, GET THE Z(AB) AS
C
C   Z(AB) = X(AB)/(f(AA)-f(BB)
C
C   D(AB) = - Z(AB)
C
C   AND AUGMENT I(A,B) IN THE USUAL WAY :   f(AA) D(A,B)
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
      DIMENSION AIVV(1),DVV(1),XVV(1),EVAL(1)
      DIMENSION SPCASE(2)
C
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPP(255,2),DIRPRD(8,8)
CJDW KKB stuff
C     COMMON /SYM/  POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYM2/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
CJDW END
C
      DATA TRESH /1.D-12/,ONE/1.D0/
      DATA SPCASE /'A','B'/
C
C  LOOP OVER SPIN CASES
C
      DO 1000 ISPIN=1,IUHF+1
C
C DETERMINE OFFSETS
C
       IVV=1+(ISPIN-1)*NFEA(1)
       IOV=1+(ISPIN-1)*NT(1)
C
       CALL ZERO(XVV(IVV),NFEA(ISPIN))
C
C ALLOCATE CORE FOR EIGENVALUES
C
CJDW KKB stuff. Beginning of decision block is new.
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
CJDW END
C      ELSE
       CALL GETREC(20,'JOBARC','SCFEVAL'//SPCASE(ISPIN),
     &             IINTFP*(NOCCO(ISPIN)+NVRTO(ISPIN)),
     &             EVAL)
CJDW KKB stuff
C       ENDIF
CJDW END
C
C  X(A,B) = I(B,A) - I(A,B)  A > B
C
C THEN SYMMETRIZE I(A,B) BY I(A,B) = I(B,A) FOR ALL A,B
C
       IOFF=0
       NOFF=NOCCO(ISPIN)
       DO 100 IRREP=1,NIRREP
C
        NVRT=VRT(IRREP,ISPIN)
C
        DO 20 IA=1,NVRT
         DO 10 IB=1,IA-1
C
         IND1=IVV+IOFF+IA+(IB-1)*NVRT-1
         IND2=IVV+IOFF+IB+(IA-1)*NVRT-1
C
         EVDIFF=EVAL(IB+NOFF)-EVAL(IA+NOFF)
C
         XVVI=AIVV(IND2)-AIVV(IND1)     
C
C Marcel observed that the DROPMO gradient of s-triazine fails as
C a result of the instability introduced by small EVDIFF. It is
C normally assumed that when EVDIFF small XOOI is also small and
C as a result the only cheked performed was for the size of XOII.               C But in this case that assumption was not accuarte and the size                C of EVDIFF must also be checked. 01/2011, Ajith Perera.
C
         IF (ABS(EVDIFF) .GT. 1.0D-9) THEN
   
            IF(ABS(XVVI).GT.TRESH) XVV(IND1)=XVVI/EVDIFF
         
         ENDIF 
C
10       CONTINUE
C
         DO 12 IB=1,IA-1
C
         IND1=IVV+IOFF+IA+(IB-1)*NVRT-1
         IND2=IVV+IOFF+IB+(IA-1)*NVRT-1
C
         AIVV(IND2)=AIVV(IND2)-EVAL(IB+NOFF)*XVV(IND1)
12       CONTINUE
C
         DO 15 IB=1,IA-1
C
         IND1=IVV+IOFF+IA+(IB-1)*NVRT-1
         IND2=IVV+IOFF+IB+(IA-1)*NVRT-1
C
         XVV(IND2)=XVV(IND1)
         AIVV(IND1)=AIVV(IND2)
C
15       CONTINUE
C
20      CONTINUE
C
        NOFF=NOFF+NVRT
        IOFF=IOFF+NVRT*NVRT
C
100    CONTINUE
C
       CALL SAXPY(NFEA(ISPIN),ONE,XVV(IVV),1,DVV(IVV),1)
C
1000  CONTINUE
C
      RETURN
      END
