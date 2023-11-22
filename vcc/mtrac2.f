      SUBROUTINE MTRAC2(F,Z,DISSYZ,POP1,POP2,IRREP,IUHF)
C
C  THIS ROUTINE CALCULATES THE CONTRACTION
C
C   F(M,I) = F(M,I) + SUM N Z(MN,IN) (UHF)
C
C   F(M,I) = F(M,I) + SUM N (TWO Z(MN,IN) - Z(NM,IN)) (RHF)
C
C   WITH N1 : DIMENSION Of M AND I
C        N2 : DIMENSION OF N
C       
C  THIS CONTRIBUTION IS REQUIRED IN CCSD FOR THE
C  CONSTRUCTION OF THE F(MI) INTERMEDIATES
C
CEND
C
C  CODED JULY/90  JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYZ,DIRPRD,POP1,POP2 
      DIMENSION F(1),Z(DISSYZ,1),POP1(8),POP2(8),IP(8),IPF(8)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      DATA ONE,TWO,ONEM /1.D0,2.0D0,-1.D0/
C
      IP(1)=0
      IPF(1)=0
      DO 1 IRREPJ=1,NIRREP-1
       IRREPI=DIRPRD(IRREP,IRREPJ)
       IP(IRREPJ+1)=IP(IRREPJ)+POP1(IRREPI)*POP2(IRREPJ)
       IPF(IRREPJ+1)=IPF(IRREPJ)+POP1(IRREPJ)*POP1(IRREPJ)
1     CONTINUE
C
      IF(IUHF.EQ.1) THEN
      IOFFZ=0
      DO 100 IRREPJ=1,NIRREP
       NOCCJ=POP2(IRREPJ)
       IRREPI=DIRPRD(IRREP,IRREPJ)
       NOCCI=POP1(IRREPI)
       INDF=IPF(IRREPI)
       DO 20 I=1,NOCCI
       DO 20 M=1,NOCCI
       INDF=INDF+1
       CALL PSADD(Z(IOFFZ+M,IOFFZ+I),ONE,(DISSYZ+1)*NOCCI,F(INDF),NOCCJ)
 20    CONTINUE
      IOFFZ=IOFFZ+NOCCI*NOCCJ
100   CONTINUE
      ELSE
      IOFFZ=0
      DO 200 IRREPJ=1,NIRREP
       NOCCJ=POP2(IRREPJ)
       IRREPI=DIRPRD(IRREP,IRREPJ)
       NOCCI=POP1(IRREPI)
       IOFFZ2=IP(IRREPI)
       INDF=IPF(IRREPI)
       DO 120 I=1,NOCCI
       DO 120 M=1,NOCCI
       INDF=INDF+1
       CALL PSADD(Z(IOFFZ+M,IOFFZ+I),TWO,(DISSYZ+1)*NOCCI,F(INDF),NOCCJ)
       CALL PSADD(Z(IOFFZ2+(M-1)*NOCCJ+1,IOFFZ+I),ONEM,DISSYZ*NOCCI+1,
     $      F(INDF),NOCCJ)
 120  CONTINUE
       IOFFZ=IOFFZ+NOCCJ*NOCCI
200    CONTINUE
       ENDIF
      RETURN
      END
