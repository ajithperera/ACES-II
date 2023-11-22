
C This routine factors N into its prime powers, lnumber in number.
C E.G., for N=1960=(2**3)*(5**1)*(7**2), lnumber=3, lfactor=2,5,7,
C lpower=8,5,49, and lexponent=3,1,2.

      SUBROUTINE prime_factors(N,lfactor,lpower,lexponent,lnumber)
      IMPLICIT NONE
      INTEGER N, lfactor(8), lpower(8), lexponent(8), lnumber

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
               lfactor(IP) = IDIV
               lpower(IP) = IDIV
               IFCUR = IDIV
               lexponent(IP) = 1
            ELSE
               lpower(IP) = IDIV*lpower(IP)
               lexponent(IP) = lexponent(IP) + 1
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
            lfactor(IP) = NPART
            lpower(IP) = NPART
            lexponent(IP) = 1
         ELSE
            lpower(IP) = NPART*lpower(IP)
            lexponent(IP) = lexponent(IP) + 1
         END IF
      END IF

      lnumber = IP
      RETURN
      END
