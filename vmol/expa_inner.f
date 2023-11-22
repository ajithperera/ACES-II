      SUBROUTINE EXPA_INNER(SUM,GP,TSF,POW,MXP2,NH4,I,I4,IR2M,IND,
     &                      NUCAB)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Arguments.
C-----------------------------------------------------------------------
      DOUBLE PRECISION SUM,GP,TSF,POW
      INTEGER MXP2,NH4,I,I4,IR2M,IND,NUCAB
C-----------------------------------------------------------------------
C     Local variables.
C-----------------------------------------------------------------------
      INTEGER IR2,I5,J
C-----------------------------------------------------------------------
      DIMENSION SUM(MXP2,3),GP(MXP2,NH4),TSF(10000),POW(3,NH4)
C
c      DO 860 IR2=IR2M,IND,2
c      I5 = I4+NUCAB*(IR2-1)
c      DO 850 J=1,NUCAB
c      SUM(J,I) = SUM(J,I) + GP(J,IR2)*TSF(J+I5)*POW(I,IND-IR2+1)
c  850 CONTINUE
c  860 CONTINUE
C
C     Above does not work at OSC with default optimization. With f90
C     compiler, default means scalar2 and vector2.
C
      DO 860 IR2=IR2M,IND,2
      I5 = I4+NUCAB*(IR2-1)
      CALL XDAB(SUM(1,I),POW(I,IND-IR2+1),GP(1,IR2),TSF(1+I5),NUCAB)
  860 CONTINUE
C
      RETURN
      END
