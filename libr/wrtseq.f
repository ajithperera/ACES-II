
C WRITES AN INTEGER VECTOR OF LENGTH IVEC TO SEQUENTIAL ACCESS
C UNIT IUNIT.

      SUBROUTINE WRTSEQ(IUNIT,IVEC,LENGTH)
      DIMENSION IVEC(LENGTH)
      WRITE(IUNIT) IVEC
      RETURN
      END