      SUBROUTINE SETRHF(FAC, IOFFST, IP)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DIMENSION IP(20), FAC(9, 9), IOFFST(15)
C
      IBTAND(I,J) = IAND(I,J)
      IBTOR(I,J)  = IOR(I,J)
      IBTXOR(I,J) = IEOR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTNOT(I)   = NOT(I)
c
      IOFFST(1) = 0
      IOFFST(3) = 1
      IOFFST(6) = 4 
      IOFFST(10) = 10
      IOFFST(15) = 20
C
      CALL SETLMN
C
      FAC(1,1)=1.0
      FAC(2,1)=1.0
      FAC(2,2)=1.0
C
      DO  I=3,9
         FAC(I,1)=1.0
         FAC(I,I)=1.0
         JE=I-1
         DO J=2,JE
            FAC(I,J)=FAC(I-1,J-1)+FAC(I-1,J)
         ENDDO
      ENDDO
C
      IP(1)=0
      IP(2)=IBTSHL(1,20)
      IP(3)=IBTSHL(1,10)
      IP(4)=1
      IP(5)=IBTSHL(2,20)
      IP(6)=IBTSHL(2,10)
      IP(7)=2
      IP(8)=IBTOR(IBTSHL(1,20),IBTSHL(1,10))
      IP(9)=IBTOR(IBTSHL(1,20),1)
      IP(10)=IBTOR(IBTSHL(1,10),1)
      IP(11)=IBTSHL(3,20)
      IP(12)=IBTSHL(3,10)
      IP(13)=3
      IP(14)=IBTOR(IBTSHL(2,20),IBTSHL(1,10))
      IP(15)=IBTOR(IBTSHL(2,20),1)
      IP(16)=IBTOR(IBTSHL(1,20),IBTSHL(2,10))
      IP(17)=IBTOR(IBTSHL(2,10),1)
      IP(18)=IBTOR(IBTSHL(1,20),2)
      IP(19)=IBTOR(IBTSHL(1,10),2)
      IP(20)=IBTOR(IBTSHL(1,20),IBTOR(IBTSHL(1,10),1))
C
      RETURN
      END
