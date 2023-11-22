
      SUBROUTINE XFORMMOS_G(COLD,CNEW,W,ANGTYP,IBAS)
      IMPLICIT NONE
      INTEGER ANGTYP(*),IBAS
      DOUBLE PRECISION COLD(*),CNEW(*),W(3,3)

c SHELL  = s p d  f  g  h  i  j   k   l   m   n   o
c ANGMOM = 0 1 2  3  4  5  6  7   8   9  10  11  12
c NCOMB  = 1 3 6 10 15 21 28 36  45  55  66  78  91
c RUNOFF = 0 1 4 10 20 35 56 84 120 165 220 286 364
      INTEGER    ANGMOM,   NCOMB,    RUNOFF
      PARAMETER (ANGMOM=4, NCOMB=15, RUNOFF=20)

      INTEGER IX, PERM(ANGMOM,NCOMB), I, J, K, L, IOFF
      DOUBLE PRECISION UNIQ(NCOMB), D_TEMP(3,3,3,3)
      DOUBLE PRECISION C_TEMP, W34, W234

      DATA ((PERM(I,J),I=1,ANGMOM),J=1,NCOMB) /
     & 1,1,1,1, ! XXXX
     & 1,1,1,2, ! XXXY
     & 1,1,1,3, ! XXXZ
     & 1,1,2,2, ! XXYY
     & 1,1,2,3, ! XXYZ
     & 1,1,3,3, ! XXZZ
     & 1,2,2,2, ! XYYY
     & 1,2,2,3, ! XYYZ
     & 1,2,3,3, ! XYZZ
     & 1,3,3,3, ! XZZZ
     & 2,2,2,2, ! YYYY
     & 2,2,2,3, ! YYYZ
     & 2,2,3,3, ! YYZZ
     & 2,3,3,3, ! YZZZ
     & 3,3,3,3  ! ZZZZ
     & /

      IX = ANGTYP(IBAS)-RUNOFF
      IOFF = IBAS-IX
      print *, 'function ',ibas,' has angmom index ',angtyp(ibas)
     &       , ' (',IX,') with offset ',IOFF

c   o unique values
      DO I = 1, NCOMB
         UNIQ(I) = COLD(IOFF+I)
      END DO

c   o permutations (let the compiler order these for efficiency...)
      D_TEMP(1,1,1,1) = UNIQ(1)

      D_TEMP(1,1,1,2) = UNIQ(2)
      D_TEMP(1,1,2,1) = UNIQ(2)
      D_TEMP(1,2,1,1) = UNIQ(2)
      D_TEMP(2,1,1,1) = UNIQ(2)

      D_TEMP(1,1,1,3) = UNIQ(3)
      D_TEMP(1,1,3,1) = UNIQ(3)
      D_TEMP(1,3,1,1) = UNIQ(3)
      D_TEMP(3,1,1,1) = UNIQ(3)

      D_TEMP(1,1,2,2) = UNIQ(4)
      D_TEMP(1,2,1,2) = UNIQ(4)
      D_TEMP(1,2,2,1) = UNIQ(4)
      D_TEMP(2,1,1,2) = UNIQ(4)
      D_TEMP(2,1,2,1) = UNIQ(4)
      D_TEMP(2,2,1,1) = UNIQ(4)

      D_TEMP(1,1,2,3) = UNIQ(5)
      D_TEMP(1,2,1,3) = UNIQ(5)
      D_TEMP(1,2,3,1) = UNIQ(5)
      D_TEMP(2,1,1,3) = UNIQ(5)
      D_TEMP(2,1,3,1) = UNIQ(5)
      D_TEMP(2,3,1,1) = UNIQ(5)
      D_TEMP(1,1,3,2) = UNIQ(5)
      D_TEMP(1,3,1,2) = UNIQ(5)
      D_TEMP(1,3,2,1) = UNIQ(5)
      D_TEMP(3,1,1,2) = UNIQ(5)
      D_TEMP(3,1,2,1) = UNIQ(5)
      D_TEMP(3,2,1,1) = UNIQ(5)

      D_TEMP(1,1,3,3) = UNIQ(6)
      D_TEMP(1,3,1,3) = UNIQ(6)
      D_TEMP(1,3,3,1) = UNIQ(6)
      D_TEMP(3,1,1,3) = UNIQ(6)
      D_TEMP(3,1,3,1) = UNIQ(6)
      D_TEMP(3,3,1,1) = UNIQ(6)

      D_TEMP(1,2,2,2) = UNIQ(7)
      D_TEMP(2,1,2,2) = UNIQ(7)
      D_TEMP(2,2,1,2) = UNIQ(7)
      D_TEMP(2,2,2,1) = UNIQ(7)

      D_TEMP(2,2,1,3) = UNIQ(8)
      D_TEMP(2,1,2,3) = UNIQ(8)
      D_TEMP(2,1,3,2) = UNIQ(8)
      D_TEMP(1,2,2,3) = UNIQ(8) ! unique D index
      D_TEMP(1,2,3,2) = UNIQ(8)
      D_TEMP(1,3,2,2) = UNIQ(8)
      D_TEMP(2,2,3,1) = UNIQ(8)
      D_TEMP(2,3,2,1) = UNIQ(8)
      D_TEMP(2,3,1,2) = UNIQ(8)
      D_TEMP(3,2,2,1) = UNIQ(8)
      D_TEMP(3,2,1,2) = UNIQ(8)
      D_TEMP(3,1,2,2) = UNIQ(8)

      D_TEMP(3,3,1,2) = UNIQ(9)
      D_TEMP(3,1,3,2) = UNIQ(9)
      D_TEMP(3,1,2,3) = UNIQ(9)
      D_TEMP(1,3,3,2) = UNIQ(9)
      D_TEMP(1,3,2,3) = UNIQ(9)
      D_TEMP(1,2,3,3) = UNIQ(9) ! unique D index
      D_TEMP(3,3,2,1) = UNIQ(9)
      D_TEMP(3,2,3,1) = UNIQ(9)
      D_TEMP(3,2,1,3) = UNIQ(9)
      D_TEMP(2,3,3,1) = UNIQ(9)
      D_TEMP(2,3,1,3) = UNIQ(9)
      D_TEMP(2,1,3,3) = UNIQ(9)

      D_TEMP(1,3,3,3) = UNIQ(10)
      D_TEMP(3,1,3,3) = UNIQ(10)
      D_TEMP(3,3,1,3) = UNIQ(10)
      D_TEMP(3,3,3,1) = UNIQ(10)

      D_TEMP(2,2,2,2) = UNIQ(11)

      D_TEMP(2,2,2,3) = UNIQ(12)
      D_TEMP(2,2,3,2) = UNIQ(12)
      D_TEMP(2,3,2,2) = UNIQ(12)
      D_TEMP(3,2,2,2) = UNIQ(12)

      D_TEMP(2,2,3,3) = UNIQ(13)
      D_TEMP(2,3,2,3) = UNIQ(13)
      D_TEMP(2,3,3,2) = UNIQ(13)
      D_TEMP(3,2,2,3) = UNIQ(13)
      D_TEMP(3,2,3,2) = UNIQ(13)
      D_TEMP(3,3,2,2) = UNIQ(13)

      D_TEMP(2,3,3,3) = UNIQ(14)
      D_TEMP(3,2,3,3) = UNIQ(14)
      D_TEMP(3,3,2,3) = UNIQ(14)
      D_TEMP(3,3,3,2) = UNIQ(14)

      D_TEMP(3,3,3,3) = UNIQ(15)

      C_TEMP = 0.d0
      do l = 1, 3
      do k = 1, 3
         W34  = W(k,PERM(3,IX)) * W(l,PERM(4,IX))
      do j = 1, 3
         W234 = W(j,PERM(2,IX)) * W34
         C_TEMP = C_TEMP + W234 * (   D_TEMP(1,j,k,l) * W(1,PERM(1,IX))
     &                              + D_TEMP(2,j,k,l) * W(2,PERM(1,IX))
     &                              + D_TEMP(3,j,k,l) * W(3,PERM(1,IX))
     &                            )
      end do
      end do
      end do
      CNEW(IBAS) = C_TEMP

      return
      end

