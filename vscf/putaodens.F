      SUBROUTINE PUTAODENS(A, IUNIT, NAOBFNS)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NAOBFNS*NAOBFNS)
      COMMON /FLAGS/ IFLAGS(100)
C
      IF(IFLAGS(1).GE.10)THEN
         WRITE(6,1000)
 1000    FORMAT(' @PUTAODENS-I, Writing AO Density to AODENS. ')
      ENDIF
C
#ifdef _DEBUG_LVLM1
      Write(6,*) "The AO density in PUTAODENS"
      call output(a, 1, naobfns, 1, naobfns, naobfns, naobfns, 1)
#endif
      NCOLS  =4 
      NTIMES = NAOBFNS/NCOLS
      NLEFT  = NAOBFNS - NTIMES*NCOLS
      IOFF   = 0
C
      IF (NTIMES .NE.0)THEN
         DO ITIMES=1,NTIMES
            DO I = 1, NAOBFNS
               WRITE(IUNIT,1020) (A(IOFF + (J-1)*NAOBFNS + I),
     &                            J=1,NCOLS)
            ENDDO
            IOFF = IOFF + NAOBFNS*NCOLS
            WRITE(IUNIT,*) 
         ENDDO
 1020 FORMAT(4F20.10)
      ENDIF
C
C Write the remainder.
C
      IF (NLEFT .NE.0) THEN
         DO I = 1, NAOBFNS
            WRITE(IUNIT,1020) (A(IOFF + (J-1)*NAOBFNS + I),
     &                         J=1,NLEFT)
         ENDDO
      ENDIF
C
      RETURN
      END
