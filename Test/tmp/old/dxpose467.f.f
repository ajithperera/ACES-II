
      SUBROUTINE DXPOSE467_F(A, N1, N2, N12, MOVED, NWORK)
      IMPLICIT NONE

      INTEGER N1, N2, N12, NWORK
      DOUBLE PRECISION A(N12)
      LOGICAL MOVED(NWORK)

      INTEGER IFACT(8), IPOWER(8), NEXP(8), IEXP(8)
      INTEGER I, N, M
      INTEGER I1, I2, I1MIN, I1MAX, I2MIN, I2MAX

      INTEGER IA1, IA2, MMIA1, MMIA2
      DOUBLE PRECISION ATEMP, BTEMP
      INTEGER NCOUNT, MMIST
      INTEGER IP, ITEST, ISTART, IDIV, ISOID
      INTEGER NPOWER
      LOGICAL REPEAT_MAIN_LOOP, REPEAT_INNER_LOOP, CYCLE_ELEMENTS
      LOGICAL EXIT_TMP_LOOP

c ---------------------------------------------

      IF ((N1.LT.2).OR.(N2.LT.2)) RETURN

      N = N1
      M = N1*N2 - 1

C   o SQUARE MATRICES
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

C   o MODULUS M FACTORED INTO PRIME POWERS
      CALL FACTOR(M,IFACT,IPOWER,NEXP,NPOWER)
      DO IP = 1, NPOWER
         IEXP(IP) = 0
      END DO

      IDIV = 1
      NCOUNT = M
      DO IP = 1, NPOWER
         IF (IEXP(IP).NE.NEXP(IP)) THEN
            NCOUNT = (NCOUNT/IFACT(IP))*(IFACT(IP)-1)
         END IF
      END DO
      DO I = 1, NWORK
         MOVED(I) = .FALSE.
      END DO
      ISTART = IDIV

      REPEAT_MAIN_LOOP = .TRUE.
      DO WHILE (REPEAT_MAIN_LOOP)
c ----------------------------------------------------------
      REPEAT_MAIN_LOOP = .FALSE.
      MMIST = M - ISTART

      IF (ISTART.NE.IDIV) THEN
         IF ((ISTART.GT.NWORK).OR.(.NOT.MOVED(ISTART))) THEN
            ISOID = ISTART/IDIV
            CYCLE_ELEMENTS = .TRUE.
            DO IP = 1, NPOWER
               IF (IEXP(IP).NE.NEXP(IP)) THEN
                  IF (MOD(ISOID,IFACT(IP)).EQ.0) CYCLE_ELEMENTS=.FALSE.
               END IF
            END DO
            IF ((ISTART.GT.NWORK).AND.CYCLE_ELEMENTS) THEN
               ITEST = ISTART
               EXIT_TMP_LOOP = .FALSE.
               DO WHILE (.NOT.EXIT_TMP_LOOP)
                  ITEST = MOD(N*ITEST,M)
                  IF ((ITEST.LT.ISTART).OR.(ITEST.GT.MMIST)) THEN
                     CYCLE_ELEMENTS = .FALSE.
                     EXIT_TMP_LOOP = .TRUE.
                  ELSE
                     IF ((ITEST.LE.ISTART).OR.(ITEST.GE.MMIST)) THEN
                        CYCLE_ELEMENTS = .TRUE.
                        EXIT_TMP_LOOP = .TRUE.
                     END IF
                  END IF
               END DO
            END IF
         ELSE
            CYCLE_ELEMENTS = .FALSE.
         END IF
      ELSE
         CYCLE_ELEMENTS = .TRUE.
      END IF

      IF (CYCLE_ELEMENTS) THEN
      ATEMP = A(ISTART+1)
      BTEMP = A(MMIST+1)
      IA1 = ISTART
      REPEAT_INNER_LOOP = .TRUE.
      DO WHILE (REPEAT_INNER_LOOP)
         REPEAT_INNER_LOOP = .FALSE.
         IA2 = MOD(N*IA1,M)
         MMIA1 = M - IA1
         MMIA2 = M - IA2
         IF (IA1.LE.NWORK)   MOVED(IA1)   = .TRUE.
         IF (MMIA1.LE.NWORK) MOVED(MMIA1) = .TRUE.
         NCOUNT = NCOUNT - 2
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
               REPEAT_INNER_LOOP = .TRUE.
            END IF
         END IF
c   o END DO WHILE (REPEAT_INNER_LOOP)
      END DO
c   o END IF (CYCLE_ELEMENTS)
      END IF

      ISTART = ISTART + IDIV
      IF (NCOUNT.LE.0) THEN
         DO IP = 1, NPOWER
            IF (.NOT.REPEAT_MAIN_LOOP) THEN
            IF (IEXP(IP).EQ.NEXP(IP)) THEN
               IEXP(IP) = 0
               IDIV = IDIV/IPOWER(IP)
            ELSE
               IEXP(IP) = IEXP(IP) + 1
               IDIV = IDIV*IFACT(IP)
               IF (IDIV.GE.M/2) RETURN
               NCOUNT = M/IDIV
               DO I = 1, NPOWER
                  IF (IEXP(I).NE.NEXP(I)) THEN
                     NCOUNT = (NCOUNT/IFACT(I))*(IFACT(I)-1)
                  END IF
               END DO
               DO I = 1, NWORK
                  MOVED(I) = .FALSE.
               END DO
               ISTART = IDIV
               REPEAT_MAIN_LOOP = .TRUE.
            END IF
            END IF
         END DO
      ELSE
         REPEAT_MAIN_LOOP = .TRUE.
      END IF

c ----------------------------------------------------------
c   o END DO WHILE (REPEAT_MAIN_LOOP)
      END DO

      RETURN
      END

c ---------------------------------------------------------

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

