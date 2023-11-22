      SUBROUTINE XFORMMOS(COLD,CNEW,W,ANGTYP,EVAL,NBASX,NBAS,IUHF)
      IMPLICIT NONE
      DOUBLE PRECISION COLD,CNEW,W,EVAL
      DOUBLE PRECISION DIFF,HALFRT2,HALF,HALFRT3,WSAVE
      DOUBLE PRECISION TWO,WTOL,THRE,C_TEMP
      INTEGER ANGTYP,NBASX,NBAS,IUHF,I,J,ISPIN,IMO,IBAS
      INTEGER IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      INTEGER IFLAGS
      LOGICAL PERMUT,NOMATCH,MATCHIJ
      DIMENSION COLD(NBASX,NBAS,IUHF+1),CNEW(NBASX),W(3,3),
     &          ANGTYP(NBASX),EVAL(NBAS*(IUHF+1))
      DIMENSION WSAVE(3,3)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/ IFLAGS(100)
      DATA TWO /2.0D+00/, THRE /3.0D+00/, WTOL /0.05D+00/
C
      IF(IFLAGS(1) .GE. 10)THEN
       write(6,*) ' @XFORMMOS-I, initial W '
       call output(w,1,3,1,3,3,3,1)
      ENDIF
C
c YAU : old
c     CALL ICOPY(9*IINTFP,W,1,WSAVE,1)
c YAU : new
      CALL DCOPY(9,W,1,WSAVE,1)
c YAU : end
      NOMATCH = .FALSE.
      DO 20 J=1,3
      DO 10 I=1,3
C
      MATCHIJ = .FALSE.
C
      DIFF = DABS(W(I,J))
      IF(DIFF.LT.WTOL)THEN
       W(I,J) = 0.0D+00
       MATCHIJ = .TRUE.
      ENDIF
C
      DIFF = DABS(W(I,J) + 1.0D+00)
      IF(DIFF.LT.WTOL)THEN
       W(I,J) = -1.0D+00
       MATCHIJ = .TRUE.
      ENDIF
C
      DIFF = DABS(W(I,J) - 1.0D+00)
      IF(DIFF.LT.WTOL)THEN
       W(I,J) =  1.0D+00
       MATCHIJ = .TRUE.
      ENDIF
C
      HALFRT2 = 0.5D+00 * DSQRT(2.0D+00)
C
      DIFF = DABS(W(I,J) + HALFRT2)
      IF(DIFF.LT.WTOL)THEN
       W(I,J) = -HALFRT2
       MATCHIJ = .TRUE.
      ENDIF
C
      DIFF = DABS(W(I,J) - HALFRT2)
      IF(DIFF.LT.WTOL)THEN
       W(I,J) =  HALFRT2
       MATCHIJ = .TRUE.
      ENDIF
C
      HALF    = 0.5D+00
C
      DIFF = DABS(W(I,J) + HALF)
      IF(DIFF.LT.WTOL)THEN
       W(I,J) =  -HALF
       MATCHIJ = .TRUE.
      ENDIF
C
      DIFF = DABS(W(I,J) - HALF)
      IF(DIFF.LT.WTOL)THEN
       W(I,J) =   HALF
       MATCHIJ = .TRUE.
      ENDIF
C
      HALFRT3 = 0.5D+00 * DSQRT(3.0D+00)
C
      DIFF = DABS(W(I,J) + HALFRT3)
      IF(DIFF.LT.WTOL)THEN
       W(I,J) =  -HALFRT3
       MATCHIJ = .TRUE.
      ENDIF
C
      DIFF = DABS(W(I,J) - HALFRT3)
      IF(DIFF.LT.WTOL)THEN
       W(I,J) =   HALFRT3
       MATCHIJ = .TRUE.
      ENDIF
C
C
C     If there is one mismatch, reset NOMATCH.
C
      IF(.NOT.MATCHIJ) NOMATCH = .TRUE.
C
   10 CONTINUE
   20 CONTINUE
C
C     If we didn't get a match, at least try to match W to a permutation
C     matrix.
C
      IF(NOMATCH)THEN
c YAU : old
c      CALL ICOPY(9*IINTFP,WSAVE,1,W,1)
c YAU : new
       CALL DCOPY(9,WSAVE,1,W,1)
c YAU : end
       CALL FINDPERM(W,PERMUT)
      ENDIF
C
C
      IF(IFLAGS(1) .GE. 10)THEN
       write(6,*) '  @XFORMMOS-I, modified W '
       call output(w,1,3,1,3,3,3,1)
      ENDIF
C
      DO 100 ISPIN=1,IUHF+1
      DO  90 IMO  =1,NBAS
      DO  80 IBAS =1,NBASX
C S
      IF(ANGTYP(IBAS) .EQ. 1) CNEW(IBAS) = COLD(IBAS,IMO,ISPIN)
C X
      IF(ANGTYP(IBAS) .EQ. 2)THEN
       CNEW(IBAS) = W(1,1) * COLD(IBAS  ,IMO,ISPIN) +
     &              W(2,1) * COLD(IBAS+1,IMO,ISPIN) +
     &              W(3,1) * COLD(IBAS+2,IMO,ISPIN)
      ENDIF
C Y
      IF(ANGTYP(IBAS) .EQ. 3)THEN
       CNEW(IBAS) = W(1,2) * COLD(IBAS-1,IMO,ISPIN) +
     &              W(2,2) * COLD(IBAS  ,IMO,ISPIN) +
     &              W(3,2) * COLD(IBAS+1,IMO,ISPIN)
      ENDIF
C Z
      IF(ANGTYP(IBAS) .EQ. 4)THEN
       CNEW(IBAS) = W(1,3) * COLD(IBAS-2,IMO,ISPIN) +
     &              W(2,3) * COLD(IBAS-1,IMO,ISPIN) +
     &              W(3,3) * COLD(IBAS  ,IMO,ISPIN)
      ENDIF
C XX
      IF(ANGTYP(IBAS) .EQ. 5)THEN
       CNEW(IBAS) =       W(1,1)*W(1,1) * COLD(IBAS  ,IMO,ISPIN) +
     &              TWO * W(1,1)*W(2,1) * COLD(IBAS+1,IMO,ISPIN) +
     &              TWO * W(1,1)*W(3,1) * COLD(IBAS+2,IMO,ISPIN) +
     &                    W(2,1)*W(2,1) * COLD(IBAS+3,IMO,ISPIN) +
     &              TWO * W(2,1)*W(3,1) * COLD(IBAS+4,IMO,ISPIN) +
     &                    W(3,1)*W(3,1) * COLD(IBAS+5,IMO,ISPIN)
      ENDIF
C XY
      IF(ANGTYP(IBAS) .EQ. 6)THEN
       CNEW(IBAS) =       W(1,1)*W(1,2)  * COLD(IBAS-1,IMO,ISPIN) +
     &   (W(1,1)*W(2,2) + W(2,1)*W(1,2)) * COLD(IBAS  ,IMO,ISPIN) +
     &   (W(1,1)*W(3,2) + W(3,1)*W(1,2)) * COLD(IBAS+1,IMO,ISPIN) +
     &                    W(2,1)*W(2,2)  * COLD(IBAS+2,IMO,ISPIN) +
     &   (W(2,1)*W(3,2) + W(3,1)*W(2,2)) * COLD(IBAS+3,IMO,ISPIN) +
     &                    W(3,1)*W(3,2)  * COLD(IBAS+4,IMO,ISPIN)
      ENDIF
C XZ
      IF(ANGTYP(IBAS) .EQ. 7)THEN
       CNEW(IBAS) =       W(1,1)*W(1,3)  * COLD(IBAS-2,IMO,ISPIN) +
     &   (W(1,1)*W(2,3) + W(2,1)*W(1,3)) * COLD(IBAS-1,IMO,ISPIN) +
     &   (W(1,1)*W(3,3) + W(3,1)*W(1,3)) * COLD(IBAS  ,IMO,ISPIN) +
     &                    W(2,1)*W(2,3)  * COLD(IBAS+1,IMO,ISPIN) +
     &   (W(2,1)*W(3,3) + W(3,1)*W(2,3)) * COLD(IBAS+2,IMO,ISPIN) +
     &                    W(3,1)*W(3,3)  * COLD(IBAS+3,IMO,ISPIN)
      ENDIF
C YY
      IF(ANGTYP(IBAS) .EQ. 8)THEN
       CNEW(IBAS) =       W(1,2)*W(1,2) * COLD(IBAS-3,IMO,ISPIN) +
     &              TWO * W(1,2)*W(2,2) * COLD(IBAS-2,IMO,ISPIN) +
     &              TWO * W(1,2)*W(3,2) * COLD(IBAS-1,IMO,ISPIN) +
     &                    W(2,2)*W(2,2) * COLD(IBAS  ,IMO,ISPIN) +
     &              TWO * W(2,2)*W(3,2) * COLD(IBAS+1,IMO,ISPIN) +
     &                    W(3,2)*W(3,2) * COLD(IBAS+2,IMO,ISPIN)
      ENDIF
C YZ
      IF(ANGTYP(IBAS) .EQ. 9)THEN
       CNEW(IBAS) =       W(1,2)*W(1,3)  * COLD(IBAS-4,IMO,ISPIN) +
     &   (W(1,2)*W(2,3) + W(2,2)*W(1,3)) * COLD(IBAS-3,IMO,ISPIN) +
     &   (W(1,2)*W(3,3) + W(3,2)*W(1,3)) * COLD(IBAS-2,IMO,ISPIN) +
     &                    W(2,2)*W(2,3)  * COLD(IBAS-1,IMO,ISPIN) +
     &   (W(2,2)*W(3,3) + W(3,2)*W(2,3)) * COLD(IBAS  ,IMO,ISPIN) +
     &                    W(3,2)*W(3,3)  * COLD(IBAS+1,IMO,ISPIN)
      ENDIF
C ZZ
      IF(ANGTYP(IBAS) .EQ.10)THEN
       CNEW(IBAS) =       W(1,3)*W(1,3) * COLD(IBAS-5,IMO,ISPIN) +
     &              TWO * W(1,3)*W(2,3) * COLD(IBAS-4,IMO,ISPIN) +
     &              TWO * W(1,3)*W(3,3) * COLD(IBAS-3,IMO,ISPIN) +
     &                    W(2,3)*W(2,3) * COLD(IBAS-2,IMO,ISPIN) +
     &              TWO * W(2,3)*W(3,3) * COLD(IBAS-1,IMO,ISPIN) +
     &                    W(3,3)*W(3,3) * COLD(IBAS  ,IMO,ISPIN)
      ENDIF
C XXX
      IF(ANGTYP(IBAS) .EQ. 11)THEN
       CNEW(IBAS) = W(1,1)*W(1,1)*W(1,1) * COLD(IBAS  ,IMO,ISPIN) +
     &       THRE * W(1,1)*W(1,1)*W(2,1) * COLD(IBAS+1,IMO,ISPIN) +
     &       THRE * W(1,1)*W(1,1)*W(3,1) * COLD(IBAS+2,IMO,ISPIN) +
     &       THRE * W(1,1)*W(2,1)*W(2,1) * COLD(IBAS+3,IMO,ISPIN) +
     &       THRE * W(1,1)*W(2,1)*W(3,1) * COLD(IBAS+4,IMO,ISPIN) +
     &       THRE * W(1,1)*W(3,1)*W(3,1) * COLD(IBAS+5,IMO,ISPIN) +
     &              W(2,1)*W(2,1)*W(2,1) * COLD(IBAS+6,IMO,ISPIN) +
     &       THRE * W(2,1)*W(2,1)*W(3,1) * COLD(IBAS+7,IMO,ISPIN) +
     &       THRE * W(2,1)*W(3,1)*W(3,1) * COLD(IBAS+8,IMO,ISPIN) +
     &              W(3,1)*W(3,1)*W(3,1) * COLD(IBAS+9,IMO,ISPIN)
      ENDIF
C XXY
      IF(ANGTYP(IBAS) .EQ. 12)THEN
       CNEW(IBAS) = W(1,1)*W(1,1)*W(1,2)  * COLD(IBAS-1,IMO,ISPIN) +
     &             (W(1,1)*W(1,1)*W(2,2)
     &        + TWO*W(1,1)*W(2,1)*W(1,2)) * COLD(IBAS,  IMO,ISPIN) +
     &             (W(1,1)*W(1,1)*W(3,2)
     &        + TWO*W(1,1)*W(3,1)*W(1,2)) * COLD(IBAS+1,IMO,ISPIN) +
     &             (W(1,1)*W(2,1)*W(2,2)
     &            + W(2,1)*W(1,1)*W(2,2)
     &            + W(2,1)*W(2,1)*W(1,2)) * COLD(IBAS+2,IMO,ISPIN) +
     &         (TWO*W(1,1)*W(2,1)*W(3,2)
     &        + TWO*W(1,1)*W(3,1)*W(2,2)
     &        + TWO*W(3,1)*W(2,1)*W(1,2)) * COLD(IBAS+3,IMO,ISPIN) +
     &         (TWO*W(1,1)*W(3,1)*W(3,2)
     &            + W(3,1)*W(3,1)*W(1,2)) * COLD(IBAS+4,IMO,ISPIN) +
     &              W(2,1)*W(2,1)*W(2,2)  * COLD(IBAS+5,IMO,ISPIN) +
     &             (W(2,1)*W(2,1)*W(3,2)
     &        + TWO*W(2,1)*W(3,1)*W(2,2)) * COLD(IBAS+6,IMO,ISPIN) +
     &         (TWO*W(2,1)*W(3,1)*W(3,2)
     &            + W(3,1)*W(3,1)*W(2,2)) * COLD(IBAS+7,IMO,ISPIN) +
     &              W(3,1)*W(3,1)*W(3,2)  * COLD(IBAS+8,IMO,ISPIN)
      ENDIF
C XXZ
      IF(ANGTYP(IBAS) .EQ. 13)THEN
       CNEW(IBAS) = W(1,1)*W(1,1)*W(1,3)  * COLD(IBAS-2,IMO,ISPIN) +
     &             (W(1,1)*W(1,1)*W(2,3)
     &        + TWO*W(1,1)*W(2,1)*W(1,3)) * COLD(IBAS-1,IMO,ISPIN) +
     &             (W(1,1)*W(1,1)*W(3,3)
     &        + TWO*W(1,1)*W(3,1)*W(1,3)) * COLD(IBAS,  IMO,ISPIN) +
     &         (TWO*W(1,1)*W(2,1)*W(2,3)
     &            + W(2,1)*W(2,1)*W(1,3)) * COLD(IBAS+1,IMO,ISPIN) +
     &         (TWO*W(1,1)*W(2,1)*W(3,3)
     &        + TWO*W(1,1)*W(3,1)*W(2,3)
     &        + TWO*W(3,1)*W(2,1)*W(1,3)) * COLD(IBAS+2,IMO,ISPIN) +
     &         (TWO*W(1,1)*W(3,1)*W(3,3)
     &            + W(3,1)*W(3,1)*W(1,3)) * COLD(IBAS+3,IMO,ISPIN) +
     &              W(2,1)*W(2,1)*W(2,3)  * COLD(IBAS+4,IMO,ISPIN) +
     &             (W(2,1)*W(2,1)*W(3,3)
     &        + TWO*W(3,1)*W(2,1)*W(2,3)) * COLD(IBAS+5,IMO,ISPIN) +
     &         (TWO*W(2,1)*W(3,1)*W(3,3)
     &            + W(3,1)*W(3,1)*W(2,3)) * COLD(IBAS+6,IMO,ISPIN) +
     &              W(3,1)*W(3,1)*W(3,3)  * COLD(IBAS+7,IMO,ISPIN)
      ENDIF
C XYY
      IF(ANGTYP(IBAS) .EQ. 14)THEN
       CNEW(IBAS) = W(1,1)*W(1,2)*W(1,2)  * COLD(IBAS-3,IMO,ISPIN) +
     &             (W(1,1)*W(1,2)*W(2,2)
     &        + TWO*W(2,1)*W(1,2)*W(1,2)) * COLD(IBAS-2,IMO,ISPIN) +
     &         (TWO*W(1,1)*W(1,2)*W(3,2)
     &            + W(3,1)*W(1,2)*W(1,2)) * COLD(IBAS-1,IMO,ISPIN) +
     &             (W(1,1)*W(2,2)*W(2,2)
     &        + TWO*W(2,1)*W(2,2)*W(1,2)) * COLD(IBAS,  IMO,ISPIN) +
     &         (TWO*W(1,1)*W(2,2)*W(3,2)
     &        + TWO*W(3,1)*W(2,2)*W(1,2)
     &        + TWO*W(2,1)*W(1,2)*W(3,2)) * COLD(IBAS+1,IMO,ISPIN) +
     &             (W(1,1)*W(3,2)*W(3,2)
     &        + TWO*W(3,1)*W(3,2)*W(1,2)) * COLD(IBAS+2,IMO,ISPIN) +
     &              W(2,1)*W(2,2)*W(2,2)  * COLD(IBAS+3,IMO,ISPIN) +
     &         (TWO*W(2,1)*W(2,2)*W(3,2)
     &            + W(3,1)*W(2,2)*W(1,2)) * COLD(IBAS+4,IMO,ISPIN) +
     &             (W(2,1)*W(3,2)*W(3,2)
     &        + TWO*W(3,1)*W(3,2)*W(2,2)) * COLD(IBAS+5,IMO,ISPIN) +
     &              W(3,1)*W(3,2)*W(3,2)  * COLD(IBAS+6,IMO,ISPIN)
      ENDIF
C XYZ
      IF(ANGTYP(IBAS) .EQ. 15)THEN
       CNEW(IBAS) = W(1,1)*W(1,2)*W(1,3)  * COLD(IBAS-4,IMO,ISPIN) +
     &             (W(1,1)*W(1,2)*W(2,3)
     &            + W(1,1)*W(2,2)*W(1,3)
     &            + W(2,1)*W(1,2)*W(1,3)) * COLD(IBAS-3,IMO,ISPIN) +
     &             (W(1,1)*W(1,2)*W(3,3)
     &            + W(1,1)*W(3,2)*W(1,3)
     &            + W(3,1)*W(1,2)*W(1,3)) * COLD(IBAS-2,IMO,ISPIN) +
     &             (W(1,1)*W(2,2)*W(2,3)
     &            + W(2,1)*W(2,2)*W(1,3)
     &            + W(2,1)*W(1,2)*W(2,3)) * COLD(IBAS-1,IMO,ISPIN) +
     &             (W(1,1)*W(2,2)*W(3,3)
     &            + W(1,1)*W(3,2)*W(2,3)
     &            + W(2,1)*W(3,2)*W(1,3)
     &            + W(2,1)*W(1,2)*W(3,3)
     &            + W(3,1)*W(1,2)*W(2,3)
     &            + W(3,1)*W(2,2)*W(1,3)) * COLD(IBAS,  IMO,ISPIN) +
     &             (W(1,1)*W(3,2)*W(3,3)
     &            + W(3,1)*W(3,2)*W(1,3)
     &            + W(3,1)*W(1,2)*W(3,3)) * COLD(IBAS+1,IMO,ISPIN) +
     &              W(2,1)*W(2,2)*W(2,3)  * COLD(IBAS+2,IMO,ISPIN) +
     &             (W(2,1)*W(2,2)*W(3,3)
     &            + W(2,1)*W(3,2)*W(2,3)
     &            + W(3,1)*W(2,2)*W(2,3)) * COLD(IBAS+3,IMO,ISPIN) +
     &             (W(2,1)*W(3,2)*W(3,3)
     &            + W(3,1)*W(3,2)*W(2,3)
     &            + W(3,1)*W(2,2)*W(3,3)) * COLD(IBAS+4,IMO,ISPIN) +
     &              W(3,1)*W(3,2)*W(3,3)  * COLD(IBAS+2,IMO,ISPIN)
      ENDIF
C XZZ
      IF(ANGTYP(IBAS) .EQ. 16)THEN
       CNEW(IBAS) = W(1,1)*W(1,3)*W(1,3)  * COLD(IBAS-5,IMO,ISPIN) +
     &         (TWO*W(1,1)*W(1,3)*W(2,3)
     &            + W(2,1)*W(1,3)*W(1,3)) * COLD(IBAS-4,IMO,ISPIN) +
     &         (TWO*W(1,1)*W(1,3)*W(3,3)
     &            + W(3,1)*W(1,3)*W(1,3)) * COLD(IBAS-3,IMO,ISPIN) +
     &             (W(1,1)*W(2,3)*W(2,3)
     &        + TWO*W(2,1)*W(2,3)*W(1,3)) * COLD(IBAS-2,IMO,ISPIN) +
     &         (TWO*W(1,1)*W(2,3)*W(3,3)
     &        + TWO*W(3,1)*W(2,3)*W(1,3)
     &        + TWO*W(2,1)*W(1,3)*W(3,3)) * COLD(IBAS-1,IMO,ISPIN) +
     &             (W(1,1)*W(3,3)*W(3,3)
     &        + TWO*W(3,1)*W(3,3)*W(1,3)) * COLD(IBAS,  IMO,ISPIN) +
     &              W(2,1)*W(2,3)*W(2,3)  * COLD(IBAS+1,IMO,ISPIN) +
     &         (TWO*W(2,1)*W(2,3)*W(3,3)
     &            + W(3,1)*W(2,3)*W(2,3)) * COLD(IBAS+2,IMO,ISPIN) +
     &             (W(2,1)*W(3,3)*W(3,3)
     &        + TWO*W(3,1)*W(3,3)*W(2,3)) * COLD(IBAS+3,IMO,ISPIN) +
     &              W(3,1)*W(3,3)*W(3,3)  * COLD(IBAS+4,IMO,ISPIN)
      ENDIF
C YYY
      IF(ANGTYP(IBAS) .EQ. 17)THEN
       CNEW(IBAS) = W(1,2)*W(1,2)*W(1,2) * COLD(IBAS-6,IMO,ISPIN) +
     &         THRE*W(1,2)*W(1,2)*W(2,2) * COLD(IBAS-5,IMO,ISPIN) +
     &         THRE*W(1,2)*W(1,2)*W(3,2) * COLD(IBAS-4,IMO,ISPIN) +
     &         THRE*W(1,2)*W(2,2)*W(2,2) * COLD(IBAS-3,IMO,ISPIN) +
     &         THRE*W(1,2)*W(2,2)*W(3,2) * COLD(IBAS-2,IMO,ISPIN) +
     &         THRE*W(1,2)*W(3,2)*W(3,2) * COLD(IBAS-1,IMO,ISPIN) +
     &              W(2,2)*W(2,2)*W(2,2) * COLD(IBAS,  IMO,ISPIN) +
     &         THRE*W(2,2)*W(2,2)*W(3,2) * COLD(IBAS+1,IMO,ISPIN) +
     &         THRE*W(2,2)*W(3,2)*W(3,2) * COLD(IBAS+2,IMO,ISPIN) +
     &              W(3,2)*W(3,2)*W(3,2) * COLD(IBAS+3,IMO,ISPIN)
      ENDIF
C YYZ
      IF(ANGTYP(IBAS) .EQ. 18)THEN
       CNEW(IBAS) = W(1,2)*W(1,2)*W(1,3)  * COLD(IBAS-7,IMO,ISPIN) +
     &             (W(1,2)*W(1,2)*W(2,3)
     &        + TWO*W(1,2)*W(2,2)*W(1,3)) * COLD(IBAS-6,IMO,ISPIN) +
     &             (W(1,2)*W(1,2)*W(3,3)
     &        + TWO*W(1,2)*W(3,2)*W(1,3)) * COLD(IBAS-5,IMO,ISPIN) +
     &         (TWO*W(1,2)*W(2,2)*W(2,3)
     &            + W(2,2)*W(2,2)*W(1,3)) * COLD(IBAS-4,IMO,ISPIN) +
     &         (TWO*W(1,2)*W(2,2)*W(3,3)
     &        + TWO*W(1,2)*W(3,2)*W(2,3)
     &        + TWO*W(3,2)*W(2,2)*W(1,3)) * COLD(IBAS-3,IMO,ISPIN) +
     &         (TWO*W(1,2)*W(3,2)*W(3,3)
     &            + W(3,2)*W(3,2)*W(1,3)) * COLD(IBAS-2,IMO,ISPIN) +
     &              W(2,2)*W(2,2)*W(2,3)  * COLD(IBAS-1,IMO,ISPIN) +
     &             (W(2,2)*W(2,2)*W(3,3)
     &        + TWO*W(2,2)*W(3,2)*W(2,3)) * COLD(IBAS,  IMO,ISPIN) +
     &         (TWO*W(2,2)*W(3,2)*W(3,3)
     &            + W(3,2)*W(3,2)*W(2,3)) * COLD(IBAS+1,IMO,ISPIN) +
     &              W(3,2)*W(3,2)*W(3,3)  * COLD(IBAS+2,IMO,ISPIN)
      ENDIF
C YZZ
      IF(ANGTYP(IBAS) .EQ. 19)THEN
       CNEW(IBAS) = W(1,2)*W(1,3)*W(1,3)  * COLD(IBAS-8,IMO,ISPIN) +
     &         (TWO*W(1,2)*W(1,3)*W(2,3)
     &            + W(2,2)*W(1,3)*W(1,3)) * COLD(IBAS-7,IMO,ISPIN) +
     &         (TWO*W(1,2)*W(1,3)*W(3,3)
     &            + W(3,2)*W(1,3)*W(1,3)) * COLD(IBAS-6,IMO,ISPIN) +
     &             (W(1,2)*W(2,3)*W(2,3)
     &        + TWO*W(2,2)*W(2,3)*W(1,3)) * COLD(IBAS-5,IMO,ISPIN) +
     &         (TWO*W(1,2)*W(2,3)*W(3,3)
     &        + TWO*W(3,2)*W(1,3)*W(2,3)
     &        + TWO*W(2,2)*W(1,3)*W(3,3)) * COLD(IBAS-4,IMO,ISPIN) +
     &             (W(1,2)*W(3,3)*W(3,3)
     &        + TWO*W(3,2)*W(3,3)*W(1,3)) * COLD(IBAS-3,IMO,ISPIN) +
     &              W(2,2)*W(2,3)*W(2,3)  * COLD(IBAS-2,IMO,ISPIN) +
     &             (W(3,2)*W(2,3)*W(2,3)
     &        + TWO*W(2,2)*W(2,3)*W(3,3)) * COLD(IBAS-1,IMO,ISPIN) +
     &             (W(2,2)*W(3,3)*W(3,3)
     &        + TWO*W(3,2)*W(3,3)*W(2,3)) * COLD(IBAS,  IMO,ISPIN) +
     &              W(3,2)*W(3,3)*W(3,3)  * COLD(IBAS+1,IMO,ISPIN)
      ENDIF
C ZZZ
      IF(ANGTYP(IBAS) .EQ. 20)THEN
       CNEW(IBAS) = W(1,3)*W(1,3)*W(1,3) * COLD(IBAS-9,IMO,ISPIN) +
     &         THRE*W(1,3)*W(1,3)*W(2,3) * COLD(IBAS-8,IMO,ISPIN) +
     &         THRE*W(1,3)*W(1,3)*W(3,3) * COLD(IBAS-7,IMO,ISPIN) +
     &         THRE*W(1,3)*W(2,3)*W(2,3) * COLD(IBAS-6,IMO,ISPIN) +
     &         THRE*W(1,3)*W(2,3)*W(3,3) * COLD(IBAS-5,IMO,ISPIN) +
     &         THRE*W(1,3)*W(3,3)*W(3,3) * COLD(IBAS-4,IMO,ISPIN) +
     &              W(2,3)*W(2,3)*W(2,3) * COLD(IBAS-3,IMO,ISPIN) +
     &         THRE*W(2,3)*W(2,3)*W(3,3) * COLD(IBAS-2,IMO,ISPIN) +
     &         THRE*W(2,3)*W(3,3)*W(3,3) * COLD(IBAS-1,IMO,ISPIN) +
     &              W(3,3)*W(3,3)*W(3,3) * COLD(IBAS,  IMO,ISPIN)
      ENDIF
C
      IF(ANGTYP(IBAS) .GE. 21 .AND. ANGTYP(IBAS) .LE. 35)THEN
       CALL XFORMMOS_G(COLD(1,IMO,ISPIN),CNEW,W,ANGTYP,IBAS)
      ENDIF
C
      IF(ANGTYP(IBAS) .GT. 35)THEN
       print *, '@XFORMMOS: AOBASMOS option is limited to g functions.'
       call errex
      ENDIF
C
   80 CONTINUE
c YAU : old
c     CALL ICOPY(NBASX*IINTFP,CNEW,1,COLD(1,IMO,ISPIN),1)
c YAU : new
      CALL DCOPY(NBASX,CNEW,1,COLD(1,IMO,ISPIN),1)
c YAU : end
   90 CONTINUE
  100 CONTINUE
C
C     Dump MOs in AO basis to JOBARC.
C
      IF(IUHF.EQ.0)THEN
       CALL PUTREC(20,'JOBARC','EVECAO_A',NBASX*NBAS*IINTFP,COLD(1,1,1))
      ELSE
       CALL PUTREC(20,'JOBARC','EVECAO_A',NBASX*NBAS*IINTFP,COLD(1,1,1))
       CALL PUTREC(20,'JOBARC','EVECAO_B',NBASX*NBAS*IINTFP,COLD(1,1,2))
      ENDIF
C
      RETURN
      END
