C  Copyright (c) 1997-1999, 2003 Massachusetts Institute of Technology
C 
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
C  USA
         SUBROUTINE  OED__PVP_1D_AB_INTEGRALS
     +
     +                    ( SHELLA,SHELLB,
     +                      NGEXCEN,
     +                      EXP2A,EXP2B,
     +                      INT1D,
     +
     +                               OUT1D )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__OVL_1D_DERV_INTEGRALS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation performs a single derivation step on
C                the input 1D overlap integrals:
C
C            I'(n,a,b) = delta (DERA,1) *
C                              [-a*I(n,a-1,b) + 2*expa(n)*I(n,a+1,b)]
C                      + delta (DERB,1) *
C                              [-b*I(n,a,b-1) + 2*expb(n)*I(n,a,b+1)]
C
C                and returns the result in a separate array.
C
C                The derivatives of the 1D integrals are calculated for
C                all the present set of exponent pairs. The values of
C                DERA and DERB can be only 1 or 0 and only one of
C                them can be equal to 1. If both of them were equal
C                to 1 this would indicate an atomic case and the above
C                derivative expression would be exactly equal to zero.
C
C
C                  Input:
C
C                    SHELLx      =  maximum shell type for centers
C                                   x = A and B after differentiation
C                    DERx        =  indicator, if differentiation is
C                                   to be performed with respect to
C                                   centers x = A and B. Two possible
C                                   values: 0 = no differentiation,
C                                   1 = differentiate, i.e. d/dx
C                    NGEXCEN,    =  current # of exponent pairs
C                                   corresponding to the contracted
C                                   shell pairs A,B
C                    EXP2x       =  all NEXP exponents x 2 for both
C                                   centers x = A and B in the
C                                   appropriate order
C                    INT1D       =  all input 1D overlap integrals
C                                   before differentiation.
C
C
C                  Output:
C
C                    OUT1D       =  all differentiated 1D overlap
C                                   integrals
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT  NONE

         INTEGER   CASE1D
         INTEGER   A,B,N,F
         INTEGER   AM,AP,BM,BP
         INTEGER   NGEXCEN
         INTEGER   SHELLA,SHELLB,SHELLP

         DOUBLE PRECISION  FOURAB
         DOUBLE PRECISION  ZERO,ONE,TWO

         DOUBLE PRECISION  EXP2A  (1:NGEXCEN)
         DOUBLE PRECISION  EXP2B  (1:NGEXCEN)

         DOUBLE PRECISION  INT1D (1:NGEXCEN,0:SHELLA+1,0:SHELLB+1)
         DOUBLE PRECISION  OUT1D (1:NGEXCEN,0:SHELLA  ,0:SHELLB  )

         PARAMETER  (ZERO   =  0.D0)
         PARAMETER  (ONE    =  1.D0)
         PARAMETER  (TWO    =  2.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...the case A = s-shell and B = s-shell.
C
C
         WRITE (*,*) "SHELLA,SHELLB: ",SHELLA,SHELLB
         WRITE (*,*) "NGEXCEN:       ",NGEXCEN

         SHELLP = SHELLA + SHELLB

         DO N = 1,NGEXCEN
             OUT1D (N,0,0) = EXP2A (N) * EXP2B (N) * INT1D (N,1,1)
         END DO
C
C
C             ...jump according to the 4 different cases that can arise:
C
C                  A-shell = s- or higher angular momentum
C                  B-shell = s- or higher angular momentum
C
C                each leading to specific simplifications.
C
C
         CASE1D = 2 * MIN (1,SHELLA) + MIN (1,SHELLB) + 1

         GOTO (1,2,3,4) CASE1D
C
C
C             ...the case A = s-shell and B = s-shell (done already,
C                immediate return).
C
C
    1    RETURN
C
C
C             ...the cases A = s-shell and B >= p-shell.
C
C
    2    CONTINUE

         WRITE (*,*) "CASE 2 in OED__PVP_1D_INTEGRALS"

         DO B = 1,SHELLB
            BP = B + 1
            BM = B - 1
            DO N = 1,NGEXCEN
               OUT1D (N,0,B) = EXP2A(N) * EXP2B(N) * INT1D(N,1,BP)
     +                              - B * EXP2A(N) * INT1D(N,1,BM)
            END DO
         END DO

         RETURN
C
C
C             ...the cases A >= p-shell and B = s-shell.
C
C
    3    CONTINUE

         WRITE (*,*) "CASE 3 in OED__PVP_1D_INTEGRALS"

         DO A = 1,SHELLA
            AP = A + 1
            AM = A - 1
            DO N = 1,NGEXCEN
               OUT1D (N,A,0) = EXP2A(N) * EXP2B(N) * INT1D(N,AP,1)
     +                              - A * EXP2B(N) * INT1D(N,AM,1)
            END DO
         END DO

         RETURN
C
C
C             ...the cases A >= p-shell and B >= p-shell.
C
C
    4    CONTINUE

         WRITE (*,*) "CASE 4 in OED__PVP_1D_INTEGRALS"

         DO B = 1,SHELLB
            BM = B - 1
            BP = B + 1
            DO A = 1,SHELLA
               AM = A - 1
               AP = A + 1
               DO N = 1,NGEXCEN
                 OUT1D (N,A,B) = EXP2A (N) * EXP2B (N) * INT1D (N,AP,BP)
     +                                - B  * EXP2A (N) * INT1D (N,AP,BM)
     +                                - A  * EXP2B (N) * INT1D (N,AM,BP)
     +                                        + A * B  * INT1D (N,AM,BM)
               END DO
            END DO
         END DO

         RETURN
C
C
C             ...ready!
C
C
         RETURN
         END
