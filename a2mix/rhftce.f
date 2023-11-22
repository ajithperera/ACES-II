      SUBROUTINE RHFTCE(C, A, E, FAC, ITYP, ITM, LMN)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C$$$      SAVE IP

      DIMENSION LMN(27)
      DIMENSION FAC(9,9),XH(9),YH(9),ZH(9)
      COMMON /HIGHL/ LMNVAL(3,84), ANORM(84)
C
      DIMENSION C(27)
      DIMENSION A(3),E(3)
C
      IBTAND(I,J) = IAND(I,J)
      IBTOR(I,J)  = IOR(I,J)
      IBTXOR(I,J) = IEOR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTNOT(I)   = NOT(I)
C
      AEX = E(1)-A(1)
      AEY = E(2)-A(2)
      AEZ = E(3)-A(3)
C
      Q = 1.00D0
      X = 1.00D0
C
      L = LMNVAL(1,ITYP)+1
      M = LMNVAL(2,ITYP)+1
      N = LMNVAL(3,ITYP)+1
C
      DO L1=1,L
         XH(L1)=X
         X=X*AEX
      ENDDO
C      
      Y=1.0D0
      DO  M1=1,M
         YH(M1)=Y
         Y=Y*AEY
      ENDDO
C
      Z=1.0D0
      DO N1=1,N
         ZH(N1)=Z
         Z=Z*AEZ
      ENDDO
C
      ITM=0
C
      DO L1=1,L
C
         X=FAC(L,L1)*XH(L1)*Q
C
         DO  M1=1,M
C
           Y=FAC(M,M1)*YH(M1)*X
C
            DO N1=1,N
C
               Z=FAC(N,N1)*ZH(N1)*Y
               ITM=ITM+1
C
               LMN(ITM)=IBTOR(IBTSHL((L-L1),20),
     &                  IBTOR(IBTSHL((M-M1),10),(N-N1)))
                 C(ITM)=Z
C
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END



