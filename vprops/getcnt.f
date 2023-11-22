      SUBROUTINE GETCNT(CENTER, VLIST, NOC, CENTR, COORD)
C
C     This routine does not compile under f90 on crunch. Apart from
C     that, it has some apparent problems. For example, in the caller
C     CENTER is the integer ILNMCM, but here it is treated as a double
C     precision variable. It is my understanding thta this routine is
C     only used in calculations which make use of VPINP (ie no supp-
C     orted ACES II calculations). Therefore, I am going to comment out
C     all of the original code and put in an error trap if this routine is
C     entered. At a later time, we may resurrect it.
C
CJDW 1/16/98.
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VLIST(NOC,10), CENTR(1), COORD(3)
C
      WRITE(6,*) '  @GETCNT-F, This routine is no longer supported. '
      WRITE(6,*) '             This calculation is aborting.        '
      WRITE(6,*) '             Please consult ACES II support.      '
      CALL ERREX
C
C Old code.
CC
CC....    SEE IF CENTER SPECIFIED BY COORD CAN BE MATCHED IN LIST
CC....    OF CENTERS
CC
C      DO 10 I = 1,NOC
C      DIST2 = 0.
C      DO 20 J = 1,3
C      DIST2 = DIST2 + (VLIST(I,J) - COORD(J))**2
C20    CONTINUE
C      IF (DIST2 .GT. 1.E-8) GOTO 10
C      CENTER = CENTR(I)
C      RETURN
C10    CONTINUE
C      DIST2 = 0.
C      DO 30 J = 1,3
C      DIST2 = DIST2 + COORD(J)**2
C30    CONTINUE
C      IF (DIST2 .LT. 1.E-8) CENTER = 6HORIGIN
      RETURN
      END
