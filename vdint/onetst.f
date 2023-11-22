      SUBROUTINE ONETST
C
C     TUH
C
C     Due to local memory hardware limits on CRAY 2 this subroutine
C     is not implemented. It would be possible to include this file
C     by compiling using the option -a static instead of -a stack.
C
C     tuh Jan. 28 1987
C     hjaaj 870521 - exclude ONETST unless requested in UPDATE
C                    in order to save memory.
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      ENTRY TSTINI
      ENTRY TSTLOP(IDTYP)
      ENTRY TSTAVE
      ENTRY TSTINT(KHKTA,KHKTB,IDENA,IDENB,ICENTA,ICENTB,ONECEN,
     *             NATOMC,LAST)
      ENTRY TSTNUM
      WRITE (LUPRI,'(/A/A)')
     *   ' ONETST has not been implemented in order to save memory.',
     *   ' Do UPDATE with *DEFINE ONETST to activate test.'
      STOP 
      END
