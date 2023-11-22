      SUBROUTINE GETSTAT(V,FIN1,FIN2,D_NEW_MAX,D_NEW_MIN,D_OLD_MAX,
     &                   D_OLD_MIN,IMAX,JMAX,IMIN,JMIN,IRREP,
     &                   IRREP_MIN,IRREP_MAX,ZQ)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER FIN1,FIN2
      DIMENSION V(FIN1,FIN2)
      DIMENSION ZQ(4)

      DO I=1,FIN1
      DO J=1,FIN2
         D_NEW_MAX = DABS(V(J,I))
         D_NEW_MIN = DABS(V(J,I))
         IF (D_NEW_MAX .GT. D_OLD_MAX) THEN
             IMAX = I
             JMAX = J
             D_OLD_MAX = D_NEW_MAX 
             IRREP_MAX = IRREP
         ENDIF
         IF (D_NEW_MIN .NE. 0.0D0) THEN
         IF (D_NEW_MIN .LT. D_OLD_MIN) THEN
             IMIN = I
             JMIN = J
             D_OLD_MIN = D_NEW_MIN
             IRREP_MIN = IRREP
         ENDIF
         ENDIF
         ZQ(3)= ZQ(3)+V(I,J)
         ZQ(4)= ZQ(4)+V(I,J)*V(J,I)
      ENDDO
      ENDDO
      
      RETURN
      END
