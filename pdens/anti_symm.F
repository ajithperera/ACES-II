      SUBROUTINE ANTI_SYM(WIN, NDIM)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C This is a simple routine to antisymmetrize a matrix.
C
      DIMENSION WIN(NDIM, NDIM)
     
      DO IDIM = 1, NDIM
         DO JDIM = 1, IDIM
    
            WIN(JDIM, IDIM) = -1.0D0*WIN(JDIM, IDIM)
       
         ENDDO
      ENDDO
C
      RETURN
      END
    
