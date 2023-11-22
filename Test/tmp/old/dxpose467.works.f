
C TRANSPOSITION OF A RECTANGULAR MATRIX IN SITU.
C BY NORMAN BRENNER, MIT, 1/72. CF. ALG. 380, CACM, 5/70.
C TRANSPOSITION OF THE N1 BY N2 MATRIX A AMOUNTS TO
C REPLACING THE ELEMENT AT VECTOR POSITION I (0=ORIGIN)
C WITH THE ELEMENT AT POSITION N1*I (MOD N1*N2-1).
C EACH SUBCYCLE OF THIS PERMUTATION IS COMPLETED IN ORDER.
C MOVED IS A LOGICAL WORK ARRAY OF LENGTH NWORK.

      SUBROUTINE DXPOSE467_WORKS(A, N1, N2, N12, MOVED, NWORK)
      IMPLICIT NONE

      INTEGER N1, N2, N12, NWORK
      DOUBLE PRECISION A(N12)
C REALLY A(N1,N2) BUT N12 = N1*N2
      LOGICAL MOVED(NWORK)

      INTEGER IFACT(8), IPOWER(8), NEXP(8), IEXP(8)
      INTEGER I, N, M
      INTEGER I1, I2, I1MIN, I1MAX, I2MIN, I2MAX

      INTEGER IA1, IA2, MMIA1, MMIA2
      DOUBLE PRECISION ATEMP, BTEMP
      INTEGER NCOUNT, MMIST
      INTEGER IP, ITEST, ISTART, IDIV, ISOID
      INTEGER NPOWER

c ----------------------------------------------------------------------

      IF ((N1.LT.2).OR.(N2.LT.2)) RETURN

      N = N1
      M = N1*N2 - 1

C SQUARE MATRICES ARE DONE SEPARATELY FOR SPEED
      IF (N1.EQ.N2) THEN
         I1MIN = 2
         DO I1MAX = N, M, N
            I2 = I1MIN + N - 1
            DO I1 = I1MIN, I1MAX
               ATEMP = A(I1)
               A(I1) = A(I2)
               A(I2) = ATEMP
               I2 = I2 + N
            END DO
            I1MIN = I1MIN + N + 1
         END DO
         RETURN
      END IF

C MODULUS M IS FACTORED INTO PRIME POWERS. EIGHT FACTORS
C SUFFICE UP TO M = 2*3*5*7*11*13*17*19 = 9,767,520.
      CALL FACTOR(M,IFACT,IPOWER,NEXP,NPOWER)
      DO IP = 1, NPOWER
         IEXP(IP) = 0
      END DO

C GENERATE EVERY DIVISOR OF M LESS THAN M/2
      IDIV = 1
   50 CONTINUE
      IF (IDIV.GE.M/2) RETURN

C THE NUMBER OF ELEMENTS WHOSE INDEX IS DIVISIBLE BY IDIV
C AND BY NO OTHER DIVISOR OF M IS THE EULER TOTIENT
C FUNCTION, PHI(M/IDIV).
      NCOUNT = M/IDIV
      DO IP = 1, NPOWER
         IF (IEXP(IP).NE.NEXP(IP)) THEN
            NCOUNT = (NCOUNT/IFACT(IP))*(IFACT(IP)-1)
         END IF
      END DO

      DO I = 1, NWORK
         MOVED(I) = .FALSE.
      END DO

C THE STARTING POINT OF A SUBCYCLE IS DIVISIBLE ONLY BY IDIV
C AND MUST NOT APPEAR IN ANY OTHER SUBCYCLE.
      ISTART = IDIV

   80 CONTINUE
      MMIST = M - ISTART

      IF (ISTART.EQ.IDIV) GOTO 120
      IF (ISTART.GT.NWORK) GOTO 90
      IF (MOVED(ISTART)) GOTO 160

   90 CONTINUE
      ISOID = ISTART/IDIV
      DO IP = 1, NPOWER
         IF (IEXP(IP).NE.NEXP(IP)) THEN
            IF (MOD(ISOID,IFACT(IP)).EQ.0) GOTO 160
         END IF
      END DO

      IF (ISTART.LE.NWORK) GOTO 120
      ITEST = ISTART
      DO
         ITEST = MOD(N*ITEST,M)
         IF ((ITEST.LT.ISTART).OR.(ITEST.GT.MMIST)) GOTO 160
         IF ((ITEST.LE.ISTART).OR.(ITEST.GE.MMIST)) GOTO 110
      END DO
  110 CONTINUE
  120 CONTINUE
      ATEMP = A(ISTART+1)
      BTEMP = A(MMIST+1)
      IA1 = ISTART
  130 CONTINUE
      IA2 = MOD(N*IA1,M)
      MMIA1 = M - IA1
      MMIA2 = M - IA2
      IF (IA1.LE.NWORK)   MOVED(IA1)   = .TRUE.
      IF (MMIA1.LE.NWORK) MOVED(MMIA1) = .TRUE.
      NCOUNT = NCOUNT - 2

C MOVE TWO ELEMENTS, THE SECOND FROM THE NEGATIVE
C SUBCYCLE. CHECK FIRST FOR SUBCYCLE CLOSURE.
      IF (ISTART.EQ.IA2) THEN
         A(IA1+1)   = ATEMP
         A(MMIA1+1) = BTEMP
      ELSE
         IF (ISTART.EQ.MMIA2) THEN
            A(IA1+1)   = BTEMP
            A(MMIA1+1) = ATEMP
         ELSE
            A(IA1+1)   = A(IA2+1)
            A(MMIA1+1) = A(MMIA2+1)
            IA1 = IA2
            GOTO 130
         END IF
      END IF
  160 CONTINUE

      ISTART = ISTART + IDIV
      IF (NCOUNT.GT.0) GOTO 80
      DO IP = 1, NPOWER
         IF (IEXP(IP).NE.NEXP(IP)) THEN
            IEXP(IP) = IEXP(IP) + 1
            IDIV = IDIV*IFACT(IP)
            GOTO 50
         END IF
         IEXP(IP) = 0
         IDIV = IDIV/IPOWER(IP)
      END DO

      RETURN
      END

c ----------------------------------------------------------------------

C FACTOR N INTO ITS PRIME POWERS, NPOWER IN NUMBER.
C E.G., FOR N=1960=2**3 *5 *7**2, NPOWER=3, IFACT=3,5,7,
C IPOWER=8,5,49, AND NEXP=3,1,2.

      SUBROUTINE FACTOR(N,IFACT,IPOWER,NEXP,NPOWER)
      IMPLICIT NONE
      INTEGER N, IFACT(8), IPOWER(8), NEXP(8), NPOWER

      INTEGER IP, IFCUR, NPART, IDIV, IQUOT
      LOGICAL STILL_GOING

      IP = 0
      IFCUR = 0
      NPART = N
      IDIV = 2

      STILL_GOING = .TRUE.
      DO WHILE (STILL_GOING)
         STILL_GOING = .FALSE.
         IQUOT = NPART/IDIV
         IF (NPART.EQ.(IDIV*IQUOT)) THEN
            IF (IDIV.GT.IFCUR) THEN
               IP = IP + 1
               IFACT(IP) = IDIV
               IPOWER(IP) = IDIV
               IFCUR = IDIV
               NEXP(IP) = 1
            ELSE
               IPOWER(IP) = IDIV*IPOWER(IP)
               NEXP(IP) = NEXP(IP) + 1
            END IF
            NPART = IQUOT
            STILL_GOING = .TRUE.
         ELSE
            IF (IQUOT.GT.IDIV) THEN
               IF (IDIV.LE.2) THEN
                  IDIV = 3
               ELSE
                  IDIV = IDIV + 2
               END IF
               STILL_GOING = .TRUE.
            END IF
         END IF
      END DO

      IF (NPART.GT.1) THEN
         IF (NPART.GT.IFCUR) THEN
            IP = IP + 1
            IFACT(IP) = NPART
            IPOWER(IP) = NPART
            NEXP(IP) = 1
         ELSE
            IPOWER(IP) = NPART*IPOWER(IP)
            NEXP(IP) = NEXP(IP) + 1
         END IF
      END IF

      NPOWER = IP
      RETURN
      END
