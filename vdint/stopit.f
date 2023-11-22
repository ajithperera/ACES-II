      SUBROUTINE STOPIT(SUB,PLACE,INT1,INT2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      CHARACTER*(*) SUB, PLACE
      WRITE (LUPRI,'(//3A)') ' Work space exceeded in subroutine ',
     *                         SUB,'.'
      IF (LEN(PLACE) .GT. 0) WRITE (LUPRI,'(/2A)') ' Location: ',PLACE
      WRITE (LUPRI,'(/A,I10)  ') ' Space required: ',MAX(INT1,INT2)
      WRITE (LUPRI,'( A,I10,/)') ' Space available:',MIN(INT1,INT2)
      STOP
      END
