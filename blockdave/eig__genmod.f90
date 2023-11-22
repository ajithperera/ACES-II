        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 25 21:54:18 2021
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EIG__genmod
          INTERFACE 
            SUBROUTINE EIG(A,B,JUNK,N,SORT)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(N,N)
              REAL(KIND=8) :: B(N,N)
              INTEGER(KIND=4) :: JUNK
              INTEGER(KIND=4) :: SORT
            END SUBROUTINE EIG
          END INTERFACE 
        END MODULE EIG__genmod
