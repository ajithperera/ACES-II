      SUBROUTINE PSIZE(MEMMAX,NBLOCK,NSLICE,NUMP,NPASS,MXBUF,MXMAX,
     1                 LSQBUF,IPRINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ************************************************
C     **  COMPUTE BUFFER SIZE AND NUMBER OF PASSES  **
C     ************************************************
C
      DIMENSION IXBUF(10), IOOPS(10)
      COMMON /PTRFIL/ JOBIN, JOBOUT, LUMC, LUSCR, LUDA, LUPSO, LUPAO
      DATA DAWEIG /1.D0/, SQWEIG /1.D0/
      SQBUF = LSQBUF
      DO 10 I = 1,10
         IXBUF(I)= MIN(MXMAX,(MEMMAX*I - 2*NBLOCK)/NSLICE)
         XBUF  = IXBUF(I)
         DAOPS = 4*NUMP*DAWEIG
         DAOPS = DAOPS/XBUF
         SQOPS = 2*NUMP*I*SQWEIG
         SQOPS = SQOPS/SQBUF
         IOOPS(I) = DAOPS + SQOPS
10    CONTINUE
      IF (IPRINT .GT. 25) THEN
         WRITE(JOBOUT,'(A,/,A,/,(3I15))')
     1    '    I/O Requests for different numbers of passes ',
     2    '           Pass   Buffer Size        I/O Ops',
     3    (I, IXBUF(I), IOOPS(I), I = 1,10)
      ENDIF
      IMAX = 0
      IOMAX = 10**9
      DO 20 I = 1,10
         IF (IXBUF(I) .LT. IOMAX) THEN
            IMAX = I
            IOMAX = IOOPS(I)
         ENDIF
20    CONTINUE
      NPASS = IMAX
      MXBUF = IXBUF(I)
      RETURN
      END
