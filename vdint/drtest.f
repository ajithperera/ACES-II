      SUBROUTINE DRTEST
C
C     Dummy routine to save memory. In order to use the
C     ".TEST  " input option under "*TWOEXP", you must UPDATE
C     this module with *DEFINE DRTEST.
C
C     15-Dec-1987 TUH+HJAaJ
C
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      ENTRY DTEST0
      ENTRY DTEST1
      ENTRY DTEST2
      WRITE (LUPRI,'(/A/A)')
     *   ' DRTEST has not been implemented in order to save memory.',
     *   ' Do UPDATE with *DEFINE DRTEST to activate test.'
      STOP 
      END
