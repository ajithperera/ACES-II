      SUBROUTINE GETAODENS(A, IUNIT, NAOBFNS)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NAOBFNS, NAOBFNS)
      COMMON /FLAGS/ IFLAGS(100)
C
      IF(IFLAGS(1).GE.10)THEN
         WRITE(6,1000)
 1000    FORMAT(' @PUTAODENS-I, Reading AO Density from AODENS. ')
      ENDIF
C
      NCOLS  =4 
      NTIMES = NAOBFNS/NCOLS
      NLEFT  = NAOBFNS - NTIMES*NCOLS
      IBEGIN = 0
      JBEGIN = 0
      JOFF   = 1
C
      IF (NTIMES .NE.0)THEN
         DO ITIMES=1,NTIMES
            DO I = 1, NAOBFNS
               READ(IUNIT,*) (A(IBEGIN+I, JBEGIN+J),
     &                        J = 1, NCOLS)
            ENDDO
            JBEGIN = JBEGIN + NCOLS
         ENDDO
      ENDIF
C
C Write the remainder.
C
      IF (NLEFT .NE.0) THEN
         DO I = 1, NAOBFNS
            READ(IUNIT,*) (A(IBEGIN+I, JBEGIN+J),
     &                     J = 1, NLEFT)
         ENDDO
      ENDIF
C
#ifdef _DEBUG_LVLM
      Write(6,*) "The AO density in GETAODENS"
      call output(A, 1, naobfns, 1, naobfns, naobfns, naobfns, 1)
#endif
C
      RETURN
      END
