      SUBROUTINE RDFPH(NO,NU,FHP,FPH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION FHP(NO,NU),FPH(NU,NO)
      NOU=NO*NU
      IASV=5*NO+2*NU+4
      READ(NALL4,REC=IASV)FHP
      DO 1 I=1,NO
      DO 1 A=1,NU
      FPH(A,I)=FHP(I,A)
 1    CONTINUE
      RETURN
      END
      SUBROUTINE RDGEN(NTAP,IREC,N,A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N)
      READ(NTAP,REC=IREC)A
      RETURN
      END
      SUBROUTINE RDIJKL(K,L,M,N,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      call ienter(22)
      NO2=NO*NO
      NO3=NO2*NO
      NO4=NO3*NO
      NU2=NU*NU
      NU3=NU2*NU
      NU4=NU3*NU
      CALL IPERM24(K,L,M,N,IPRM)
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     *23,24)IPRM
 1    CONTINUE
      IAS=IT4(K,L,M,N)
      CALL RDT4(IAS,NU,YT)
      GO TO 100
 2    CONTINUE
      IAS=IT4(K,L,N,M)
      CALL RDT4(IAS,NU,Y)
      DO 51 A=1,NU
      DO 51 B=1,NU
      DO 51 C=1,NU
      DO 51 D=1,NU
      YT(A,B,D,C)=Y(A,B,C,D)
 51   CONTINUE
      GO TO 100
 3    CONTINUE
      IAS=IT4(K,M,L,N)
      CALL RDT4(IAS,NU,Y)
      DO 52 A=1,NU
      DO 52 B=1,NU
      DO 52 C=1,NU
      DO 52 D=1,NU
      YT(A,C,B,D)=Y(A,B,C,D)
 52   CONTINUE
      GO TO 100
 4    CONTINUE
      IAS=IT4(K,N,L,M)
      CALL RDT4(IAS,NU,Y)
      DO 53 A=1,NU
      DO 53 B=1,NU
      DO 53 C=1,NU
      DO 53 D=1,NU
      YT(A,C,D,B)=Y(A,B,C,D)
 53   CONTINUE
      GO TO 100
 5    CONTINUE
      IAS=IT4(K,M,N,L)
      CALL RDT4(IAS,NU,Y)
      DO 54 A=1,NU
      DO 54 B=1,NU
      DO 54 C=1,NU
      DO 54 D=1,NU
      YT(A,D,B,C)=Y(A,B,C,D)
 54   CONTINUE
      GO TO 100
 6    CONTINUE
      IAS=IT4(K,N,M,L)
      CALL RDT4(IAS,NU,Y)
      DO 55 A=1,NU
      DO 55 B=1,NU
      DO 55 C=1,NU
      DO 55 D=1,NU
      YT(A,D,C,B)=Y(A,B,C,D)
 55   CONTINUE
      GO TO 100
 7    CONTINUE
      IAS=IT4(M,K,L,N)
      CALL RDT4(IAS,NU,Y)
      DO 56 A=1,NU
      DO 56 B=1,NU
      DO 56 C=1,NU
      DO 56 D=1,NU
      YT(B,C,A,D)=Y(A,B,C,D)
 56   CONTINUE
      GO TO 100
 8    CONTINUE
      IAS=IT4(N,K,L,M)
      CALL RDT4(IAS,NU,Y)
      DO 57 A=1,NU
      DO 57 B=1,NU
      DO 57 C=1,NU
      DO 57 D=1,NU
      YT(B,C,D,A)=Y(A,B,C,D)
 57   CONTINUE
      GO TO 100
 9    CONTINUE
      IAS=IT4(M,K,N,L)
      CALL RDT4(IAS,NU,Y)
      DO 58 A=1,NU
      DO 58 B=1,NU
      DO 58 C=1,NU
      DO 58 D=1,NU
      YT(B,D,A,C)=Y(A,B,C,D)
 58   CONTINUE
      GO TO 100
 10   CONTINUE
      IAS=IT4(N,K,M,L)
      CALL RDT4(IAS,NU,Y)
      DO 59 A=1,NU
      DO 59 B=1,NU
      DO 59 C=1,NU
      DO 59 D=1,NU
      YT(B,D,C,A)=Y(A,B,C,D)
 59   CONTINUE
      GO TO 100
 11   CONTINUE
      IAS=IT4(M,N,K,L)
      CALL RDT4(IAS,NU,Y)
      DO 60 A=1,NU
      DO 60 B=1,NU
      DO 60 C=1,NU
      DO 60 D=1,NU
      YT(C,D,A,B)=Y(A,B,C,D)
 60   CONTINUE
      GO TO 100
 12   CONTINUE
      IAS=IT4(N,M,K,L)
      CALL RDT4(IAS,NU,Y)
      DO 61 A=1,NU
      DO 61 B=1,NU
      DO 61 C=1,NU
      DO 61 D=1,NU
      YT(C,D,B,A)=Y(A,B,C,D)
 61   CONTINUE
      GO TO 100
 13   CONTINUE
      IAS=IT4(L,K,M,N)
      CALL RDT4(IAS,NU,Y)
      DO 150 A=1,NU
      DO 150 B=1,NU
      DO 150 C=1,NU
      DO 150 D=1,NU
      YT(B,A,C,D)=Y(A,B,C,D)
 150   CONTINUE
      GO TO 100
 14   CONTINUE
      IAS=IT4(L,K,N,M)
      CALL RDT4(IAS,NU,Y)
      DO 151 A=1,NU
      DO 151 B=1,NU
      DO 151 C=1,NU
      DO 151 D=1,NU
      YT(B,A,D,C)=Y(A,B,C,D)
 151   CONTINUE
      GO TO 100
 15   CONTINUE
      IAS=IT4(L,M,K,N)
      CALL RDT4(IAS,NU,Y)
      DO 152 A=1,NU
      DO 152 B=1,NU
      DO 152 C=1,NU
      DO 152 D=1,NU
      YT(C,A,B,D)=Y(A,B,C,D)
 152   CONTINUE
      GO TO 100
 16   CONTINUE
      IAS=IT4(L,N,K,M)
      CALL RDT4(IAS,NU,Y)
      DO 153 A=1,NU
      DO 153 B=1,NU
      DO 153 C=1,NU
      DO 153 D=1,NU
      YT(C,A,D,B)=Y(A,B,C,D)
 153   CONTINUE
      GO TO 100
 17   CONTINUE
      IAS=IT4(L,M,N,K)
      CALL RDT4(IAS,NU,Y)
      DO 154 A=1,NU
      DO 154 B=1,NU
      DO 154 C=1,NU
      DO 154 D=1,NU
      YT(D,A,B,C)=Y(A,B,C,D)
 154   CONTINUE
      GO TO 100
 18   CONTINUE
      IAS=IT4(L,N,M,K)
      CALL RDT4(IAS,NU,Y)
      DO 155 A=1,NU
      DO 155 B=1,NU
      DO 155 C=1,NU
      DO 155 D=1,NU
      YT(D,A,C,B)=Y(A,B,C,D)
 155   CONTINUE
      GO TO 100
 19   CONTINUE
      IAS=IT4(M,L,K,N)
      CALL RDT4(IAS,NU,Y)
      DO 156 A=1,NU
      DO 156 B=1,NU
      DO 156 C=1,NU
      DO 156 D=1,NU
      YT(C,B,A,D)=Y(A,B,C,D)
 156   CONTINUE
      GO TO 100
 20   CONTINUE
      IAS=IT4(N,L,K,M)
      CALL RDT4(IAS,NU,Y)
      DO 157 A=1,NU
      DO 157 B=1,NU
      DO 157 C=1,NU
      DO 157 D=1,NU
      YT(C,B,D,A)=Y(A,B,C,D)
 157   CONTINUE
      GO TO 100
 21   CONTINUE
      IAS=IT4(M,L,N,K)
      CALL RDT4(IAS,NU,Y)
      DO 158 A=1,NU
      DO 158 B=1,NU
      DO 158 C=1,NU
      DO 158 D=1,NU
      YT(D,B,A,C)=Y(A,B,C,D)
 158   CONTINUE
      GO TO 100
 22   CONTINUE
      IAS=IT4(N,L,M,K)
      CALL RDT4(IAS,NU,Y)
      DO 159 A=1,NU
      DO 159 B=1,NU
      DO 159 C=1,NU
      DO 159 D=1,NU
      YT(D,B,C,A)=Y(A,B,C,D)
 159   CONTINUE
      GO TO 100
 23   CONTINUE
      IAS=IT4(M,N,L,K)
      CALL RDT4(IAS,NU,Y)
      DO 160 A=1,NU
      DO 160 B=1,NU
      DO 160 C=1,NU
      DO 160 D=1,NU
      YT(D,C,A,B)=Y(A,B,C,D)
 160   CONTINUE
      GO TO 100
 24   CONTINUE
      IAS=IT4(N,M,L,K)
      CALL RDT4(IAS,NU,Y)
      DO 161 A=1,NU
      DO 161 B=1,NU
      DO 161 C=1,NU
      DO 161 D=1,NU
      YT(D,C,B,A)=Y(A,B,C,D)
 161   CONTINUE
      GO TO 100
 100  CONTINUE
      call iexit(22)
      RETURN
      END
      SUBROUTINE RDIJKLT3(K,L,M,N,NO,NU,YT,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      DATA TWO/2.0D+0/
      NO2=NO*NO
      NO3=NO2*NO
      NO4=NO3*NO
      NU2=NU*NU
      NU3=NU2*NU
      NU4=NU3*NU
      CALL IPERM24(K,L,M,N,IPRM)
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     *23,24)IPRM
 1    CONTINUE
      IAS=IT4(K,L,M,N)
      CALL RDT4(IAS,NU,Y)
      DO 50 A=1,NU
      DO 50 B=1,NU
      DO 50 C=1,NU
      DO 50 D=1,NU
      YT(A,B,C,D)=Y(A,B,C,D)
 50   CONTINUE
      GO TO 100
 2    CONTINUE
      IAS=IT4(K,L,N,M)
      CALL RDT4(IAS,NU,Y)
      DO 51 A=1,NU
      DO 51 B=1,NU
      DO 51 C=1,NU
      DO 51 D=1,NU
      YT(A,B,D,C)=Y(A,B,C,D)
 51   CONTINUE
      GO TO 100
 3    CONTINUE
      IAS=IT4(K,M,L,N)
      CALL RDT4(IAS,NU,Y)
      DO 52 A=1,NU
      DO 52 B=1,NU
      DO 52 C=1,NU
      DO 52 D=1,NU
      YT(A,C,B,D)=Y(A,B,C,D)
 52   CONTINUE
      GO TO 100
 4    CONTINUE
      IAS=IT4(K,N,L,M)
      CALL RDT4(IAS,NU,Y)
      DO 53 A=1,NU
      DO 53 B=1,NU
      DO 53 C=1,NU
      DO 53 D=1,NU
      YT(A,C,D,B)=Y(A,B,C,D)
 53   CONTINUE
      GO TO 100
 5    CONTINUE
      IAS=IT4(K,M,N,L)
      CALL RDT4(IAS,NU,Y)
      DO 54 A=1,NU
      DO 54 B=1,NU
      DO 54 C=1,NU
      DO 54 D=1,NU
      YT(A,D,B,C)=Y(A,B,C,D)
 54   CONTINUE
      GO TO 100
 6    CONTINUE
      IAS=IT4(K,N,M,L)
      CALL RDT4(IAS,NU,Y)
      DO 55 A=1,NU
      DO 55 B=1,NU
      DO 55 C=1,NU
      DO 55 D=1,NU
      YT(A,D,C,B)=Y(A,B,C,D)
 55   CONTINUE
      GO TO 100
 7    CONTINUE
      IAS=IT4(M,K,L,N)
      CALL RDT4(IAS,NU,Y)
      DO 56 A=1,NU
      DO 56 B=1,NU
      DO 56 C=1,NU
      DO 56 D=1,NU
      YT(B,C,A,D)=Y(A,B,C,D)
 56   CONTINUE
      GO TO 100
 8    CONTINUE
      IAS=IT4(N,K,L,M)
      CALL RDT4(IAS,NU,Y)
      DO 57 A=1,NU
      DO 57 B=1,NU
      DO 57 C=1,NU
      DO 57 D=1,NU
      YT(B,C,D,A)=Y(A,B,C,D)
 57   CONTINUE
      GO TO 100
 9    CONTINUE
      IAS=IT4(M,K,N,L)
      CALL RDT4(IAS,NU,Y)
      DO 58 A=1,NU
      DO 58 B=1,NU
      DO 58 C=1,NU
      DO 58 D=1,NU
      YT(B,D,A,C)=Y(A,B,C,D)
 58   CONTINUE
      GO TO 100
 10   CONTINUE
      IAS=IT4(N,K,M,L)
      CALL RDT4(IAS,NU,Y)
      DO 59 A=1,NU
      DO 59 B=1,NU
      DO 59 C=1,NU
      DO 59 D=1,NU
      YT(B,D,C,A)=Y(A,B,C,D)
 59   CONTINUE
      GO TO 100
 11   CONTINUE
      IAS=IT4(M,N,K,L)
      CALL RDT4(IAS,NU,Y)
      DO 60 A=1,NU
      DO 60 B=1,NU
      DO 60 C=1,NU
      DO 60 D=1,NU
      YT(C,D,A,B)=Y(A,B,C,D)
 60   CONTINUE
      GO TO 100
 12   CONTINUE
      IAS=IT4(N,M,K,L)
      CALL RDT4(IAS,NU,Y)
      DO 61 A=1,NU
      DO 61 B=1,NU
      DO 61 C=1,NU
      DO 61 D=1,NU
      YT(C,D,B,A)=Y(A,B,C,D)
 61   CONTINUE
      GO TO 100
 13   CONTINUE
      IAS=IT4(L,K,M,N)
      CALL RDT4(IAS,NU,Y)
      DO 150 A=1,NU
      DO 150 B=1,NU
      DO 150 C=1,NU
      DO 150 D=1,NU
      YT(B,A,C,D)=Y(A,B,C,D)
 150   CONTINUE
      GO TO 100
 14   CONTINUE
      IAS=IT4(L,K,N,M)
      CALL RDT4(IAS,NU,Y)
      DO 151 A=1,NU
      DO 151 B=1,NU
      DO 151 C=1,NU
      DO 151 D=1,NU
      YT(B,A,D,C)=Y(A,B,C,D)
 151   CONTINUE
      GO TO 100
 15   CONTINUE
      IAS=IT4(L,M,K,N)
      CALL RDT4(IAS,NU,Y)
      DO 152 A=1,NU
      DO 152 B=1,NU
      DO 152 C=1,NU
      DO 152 D=1,NU
      YT(C,A,B,D)=Y(A,B,C,D)
 152   CONTINUE
      GO TO 100
 16   CONTINUE
      IAS=IT4(L,N,K,M)
      CALL RDT4(IAS,NU,Y)
      DO 153 A=1,NU
      DO 153 B=1,NU
      DO 153 C=1,NU
      DO 153 D=1,NU
      YT(C,A,D,B)=Y(A,B,C,D)
 153   CONTINUE
      GO TO 100
 17   CONTINUE
      IAS=IT4(L,M,N,K)
      CALL RDT4(IAS,NU,Y)
      DO 154 A=1,NU
      DO 154 B=1,NU
      DO 154 C=1,NU
      DO 154 D=1,NU
      YT(D,A,B,C)=Y(A,B,C,D)
 154   CONTINUE
      GO TO 100
 18   CONTINUE
      IAS=IT4(L,N,M,K)
      CALL RDT4(IAS,NU,Y)
      DO 155 A=1,NU
      DO 155 B=1,NU
      DO 155 C=1,NU
      DO 155 D=1,NU
      YT(D,A,C,B)=Y(A,B,C,D)
 155   CONTINUE
      GO TO 100
 19   CONTINUE
      IAS=IT4(M,L,K,N)
      CALL RDT4(IAS,NU,Y)
      DO 156 A=1,NU
      DO 156 B=1,NU
      DO 156 C=1,NU
      DO 156 D=1,NU
      YT(C,B,A,D)=Y(A,B,C,D)
 156   CONTINUE
      GO TO 100
 20   CONTINUE
      IAS=IT4(N,L,K,M)
      CALL RDT4(IAS,NU,Y)
      DO 157 A=1,NU
      DO 157 B=1,NU
      DO 157 C=1,NU
      DO 157 D=1,NU
      YT(C,B,D,A)=Y(A,B,C,D)
 157   CONTINUE
      GO TO 100
 21   CONTINUE
      IAS=IT4(M,L,N,K)
      CALL RDT4(IAS,NU,Y)
      DO 158 A=1,NU
      DO 158 B=1,NU
      DO 158 C=1,NU
      DO 158 D=1,NU
      YT(D,B,A,C)=Y(A,B,C,D)
 158   CONTINUE
      GO TO 100
 22   CONTINUE
      IAS=IT4(N,L,M,K)
      CALL RDT4(IAS,NU,Y)
      DO 159 A=1,NU
      DO 159 B=1,NU
      DO 159 C=1,NU
      DO 159 D=1,NU
      YT(D,B,C,A)=Y(A,B,C,D)
 159   CONTINUE
      GO TO 100
 23   CONTINUE
      IAS=IT4(M,N,L,K)
      CALL RDT4(IAS,NU,Y)
      DO 160 A=1,NU
      DO 160 B=1,NU
      DO 160 C=1,NU
      DO 160 D=1,NU
      YT(D,C,A,B)=Y(A,B,C,D)
 160   CONTINUE
      GO TO 100
 24   CONTINUE
      IAS=IT4(N,M,L,K)
      CALL RDT4(IAS,NU,Y)
      DO 161 A=1,NU
      DO 161 B=1,NU
      DO 161 C=1,NU
      DO 161 D=1,NU
      YT(D,C,B,A)=Y(A,B,C,D)
 161   CONTINUE
      GO TO 100
 100  CONTINUE
      do 300 a=1,nu
      do 300 b=1,nu
      do 300 c=1,nu
      do 300 d=1,nu
      Y(A,B,C,D)=TWO*YT(A,B,C,D)-YT(A,B,D,C)-YT(A,D,C,B)-YT(D,B,C,A)
  300 CONTINUE
      RETURN
      END
      SUBROUTINE RDIJKM(I,J,K,M,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      common/nn/no2,no3,no4,nu2,nu3,nu4,nou,no2u,no3u,nou2,nou3,no2u2
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      call ienter(23)
      CALL IPERM4(I,J,K,M,IPRM)
      GO TO (1,2,3,4)IPRM
 1    CONTINUE
      IAS=IT4(I,J,K,M)
      CALL RDT4(IAS,NU,YT)
      GO TO 100
 2    CONTINUE
      IAS=IT4(I,J,M,K)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,1)
      GO TO 100
 3    CONTINUE
      IAS=IT4(I,M,J,K)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,3)
      GO TO 100
 4    CONTINUE
      IAS=IT4(M,I,J,K)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,21)
      GO TO 100
 100  CONTINUE
      call iexit(23)
      RETURN
      END
      SUBROUTINE RDIJKMsm(I,J,K,M,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      common/nn/no2,no3,no4,nu2,nu3,nu4,nou,no2u,no3u,nou2,nou3,no2u2
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      data two/2.0d0/
      call ienter(23)
      CALL IPERM4(I,J,K,M,IPRM)
      GO TO (1,2,3,4)IPRM
 1    CONTINUE
      IAS=IT4(I,J,K,M)
      CALL RDT4(IAS,NU,Y)
      DO 50 A=1,NU
      DO 50 B=1,NU
      DO 50 C=1,NU
      DO 50 D=1,NU
      YT(A,B,C,D)=Y(A,B,C,D)*two-Y(A,B,D,C)-Y(A,D,C,B)-Y(D,B,C,A)
 50   CONTINUE
      GO TO 100
 2    CONTINUE
      IAS=IT4(I,J,M,K)
      CALL RDT4(IAS,NU,Y)
      DO 51 A=1,NU
      DO 51 B=1,NU
      DO 51 C=1,NU
      DO 51 D=1,NU
      YT(A,B,D,C)=Y(A,B,C,D)*two-Y(A,B,D,C)-Y(A,C,B,D)-Y(C,B,A,D)
 51   CONTINUE
      GO TO 100
 3    CONTINUE
      IAS=IT4(I,M,J,K)
      CALL RDT4(IAS,NU,Y)
      DO 52 A=1,NU
      DO 52 B=1,NU
      DO 52 C=1,NU
      DO 52 D=1,NU
      YT(A,C,D,B)=Y(A,B,C,D)*TWO-Y(A,D,C,B)-Y(A,C,B,D)-Y(B,A,C,D)
 52   CONTINUE
      GO TO 100
 4    CONTINUE
      IAS=IT4(M,I,J,K)
      CALL RDT4(IAS,NU,Y)
      DO 53 A=1,NU
      DO 53 B=1,NU
      DO 53 C=1,NU
      DO 53 D=1,NU
      YT(B,C,D,A)=Y(A,B,C,D)*TWO-Y(D,B,C,A)-Y(C,B,A,D)-Y(B,A,C,D)
 53   CONTINUE
      GO TO 100
 100  CONTINUE
      call iexit(23)
      RETURN
      END
      SUBROUTINE RDIJMK(I,J,K,M,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      common/nn/no2,no3,no4,nu2,nu3,nu4,nou,no2u,no3u,nou2,nou3,no2u2
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      call ienter(24)
      CALL IPERM4(I,J,K,M,IPRM)
      GO TO (1,2,3,4)IPRM
 1    CONTINUE
      IAS=IT4(I,J,K,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,1)
      GO TO 100
 2    CONTINUE
      IAS=IT4(I,J,M,K)
      CALL RDT4(IAS,NU,YT)
      GO TO 100
 3    CONTINUE
      IAS=IT4(I,M,J,K)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,7)
      GO TO 100
 4    CONTINUE
      IAS=IT4(M,I,J,K)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,13)
 100  CONTINUE
      call iexit(24)
      RETURN
      END
      SUBROUTINE RDIJMKsm(I,J,K,M,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      common/nn/no2,no3,no4,nu2,nu3,nu4,nou,no2u,no3u,nou2,nou3,no2u2
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      data two/2.0d+0/
      call ienter(24)
      nu4=nu*nu*nu*nu
      CALL IPERM4(I,J,K,M,IPRM)
      GO TO (1,2,3,4)IPRM
 1    CONTINUE
      IAS=IT4(I,J,K,M)
      CALL RDT4(IAS,NU,Y)
      DO 50 A=1,NU
      DO 50 B=1,NU
      DO 50 C=1,NU
      DO 50 D=1,NU
      YT(A,B,C,D)=Y(A,B,C,D)*two-y(a,b,D,C)-Y(D,B,C,A)-Y(A,D,C,B)
 50   CONTINUE
      GO TO 100
 2    CONTINUE
      IAS=IT4(I,J,M,K)
      CALL RDT4(IAS,NU,Y)
      DO 51 A=1,NU
      DO 51 B=1,NU
      DO 51 C=1,NU
      DO 51 D=1,NU
      YT(A,B,C,D)=Y(A,B,D,C)*two-y(a,b,C,d)-Y(D,B,A,C)-Y(A,D,B,C)
 51   CONTINUE
      GO TO 100
 3    CONTINUE
      IAS=IT4(I,M,J,K)
      CALL RDT4(IAS,NU,Y)
      DO 52 A=1,NU
      DO 52 B=1,NU
      DO 52 C=1,NU
      DO 52 D=1,NU
      YT(A,B,C,D)=Y(A,D,B,C)*TWO-Y(A,C,B,D)-Y(D,A,B,C)-Y(A,B,D,C)
 52   CONTINUE
      GO TO 100
 4    CONTINUE
      IAS=IT4(M,I,J,K)
      CALL RDT4(IAS,NU,Y)
      DO 53 A=1,NU
      DO 53 B=1,NU
      DO 53 C=1,NU
      DO 53 D=1,NU
      YT(A,B,C,D)=Y(D,A,B,C)*TWO-Y(C,A,B,D)-Y(A,D,B,C)-Y(B,A,D,C)
 53   CONTINUE
      GO TO 100
 100  CONTINUE
      call iexit(24)
      RETURN
      END
      SUBROUTINE RDIJMN(K,L,M,N,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      common/nn/no2,no3,no4,nu2,nu3,nu4,nou,no2u,no3u,nou2,nou3,no2u2
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      call ienter(25)
      CALL IPERM(K,L,M,N,IPRM)
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12)IPRM
 1    CONTINUE
      IAS=IT4(K,L,M,N)
      CALL RDT4(IAS,NU,YT)
      GO TO 100
 2    CONTINUE
      IAS=IT4(K,L,N,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,1)
      GO TO 100
 3    CONTINUE
      IAS=IT4(K,M,L,N)
      CALL RDT4(IAS,NU,Yt)
      call mtrans(yt,nu,7)
      GO TO 100
 4    CONTINUE
      IAS=IT4(K,N,L,M)
      CALL RDT4(IAS,NU,Yt)
      call mtrans(yt,nu,3)
      GO TO 100
 5    CONTINUE
      IAS=IT4(K,M,N,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,2)
      GO TO 100
 6    CONTINUE
      IAS=IT4(K,N,M,L)
      CALL RDT4(IAS,NU,Yt)
      call mtrans(yt,nu,10)
      GO TO 100
 7    CONTINUE
      IAS=IT4(M,K,L,N)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,13)
      GO TO 100
 8    CONTINUE
      IAS=IT4(N,K,L,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,21)
      GO TO 100
 9    CONTINUE
      IAS=IT4(M,K,N,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,23)
      GO TO 100
 10   CONTINUE
      IAS=IT4(N,K,M,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,16)
      GO TO 100
 11   CONTINUE
      IAS=IT4(M,N,K,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,9)
      GO TO 100
 12   CONTINUE
      IAS=IT4(N,M,K,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,17)
      GO TO 100
 100  CONTINUE
      call iexit(25)
      RETURN
      END
      SUBROUTINE RDIJMNT3(K,L,M,N,NO,NU,YT,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      common/nn/no2,no3,no4,nu2,nu3,nu4,nou,no2u,no3u,nou2,nou3,no2u2
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      DATA TWO/2.0D+0/
      CALL IPERM(K,L,M,N,IPRM)
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12)IPRM
 1    CONTINUE
      IAS=IT4(K,L,M,N)
      CALL RDT4(IAS,NU,YT)
      GO TO 100
 2    CONTINUE
      IAS=IT4(K,L,N,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,1)
      GO TO 100
 3    CONTINUE
      IAS=IT4(K,M,L,N)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,7)
      GO TO 100
 4    CONTINUE
      IAS=IT4(K,N,L,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,3)
      GO TO 100
 5    CONTINUE
      IAS=IT4(K,M,N,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,2)
      GO TO 100
 6    CONTINUE
      IAS=IT4(K,N,M,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,10)
      GO TO 100
 7    CONTINUE
      IAS=IT4(M,K,L,N)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,13)
      GO TO 100
 8    CONTINUE
      IAS=IT4(N,K,L,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,21)
      GO TO 100
 9    CONTINUE
      IAS=IT4(M,K,N,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,23)
      GO TO 100
 10   CONTINUE
      IAS=IT4(N,K,M,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,16)
      GO TO 100
 11   CONTINUE
      IAS=IT4(M,N,K,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,9)
      GO TO 100
 12   CONTINUE
      IAS=IT4(N,M,K,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,17)
      GO TO 100
 100  CONTINUE
      do 80 d=1,nu
      do 80 c=1,nu
      do 80 b=1,nu
      do 80 a=1,nu
      Y(A,B,C,D)=TWO*YT(A,B,C,D)-YT(A,B,D,C)-YT(A,D,C,B)-YT(D,B,C,A)
 80   CONTINUE
      RETURN
      END
      SUBROUTINE RDILKM(I,J,K,M,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      common/nn/no2,no3,no4,nu2,nu3,nu4,nou,no2u,no3u,nou2,nou3,no2u2
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      call ienter(26)
      CALL IPERM4(I,J,K,M,IPRM)
      GO TO (1,2,3,4)IPRM
 1    CONTINUE
      IAS=IT4(I,J,K,M)
      CALL RDT4(IAS,NU,Y)
      DO 50 A=1,NU
      DO 50 B=1,NU
      DO 50 C=1,NU
      DO 50 D=1,NU
      YT(A,C,B,D)=Y(A,B,C,D)
 50   CONTINUE
      GO TO 100
 2    CONTINUE
      IAS=IT4(I,J,M,K)
      CALL RDT4(IAS,NU,Y)
      DO 51 A=1,NU
      DO 51 B=1,NU
      DO 51 C=1,NU
      DO 51 D=1,NU
      YT(A,D,B,C)=Y(A,B,C,D)
 51   CONTINUE
      GO TO 100
 3    CONTINUE
      IAS=IT4(I,M,J,K)
      CALL RDT4(IAS,NU,Y)
      DO 52 A=1,NU
      DO 52 B=1,NU
      DO 52 C=1,NU
      DO 52 D=1,NU
      YT(A,D,C,B)=Y(A,B,C,D)
 52   CONTINUE
      GO TO 100
 4    CONTINUE
      IAS=IT4(M,I,J,K)
      CALL RDT4(IAS,NU,Y)
      DO 53 A=1,NU
      DO 53 B=1,NU
      DO 53 C=1,NU
      DO 53 D=1,NU
      YT(B,D,C,A)=Y(A,B,C,D)
 53   CONTINUE
      GO TO 100
 100  CONTINUE
      call iexit(26)
      RETURN
      END
      SUBROUTINE RDIMJK(I,J,K,M,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      common/nn/no2,no3,no4,nu2,nu3,nu4,nou,no2u,no3u,nou2,nou3,no2u2
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      call ienter(27)
      CALL IPERM4(I,J,K,M,IPRM)
      GO TO (1,2,3,4)IPRM
 1    CONTINUE
      IAS=IT4(I,J,K,M)
      CALL RDT4(IAS,NU,Y)
      DO 50 A=1,NU
      DO 50 B=1,NU
      DO 50 C=1,NU
      DO 50 D=1,NU
      YT(A,B,C,D)=Y(A,C,D,B)
 50   CONTINUE
      GO TO 100
 2    CONTINUE
      IAS=IT4(I,J,M,K)
      CALL RDT4(IAS,NU,Y)
      DO 51 A=1,NU
      DO 51 B=1,NU
      DO 51 C=1,NU
      DO 51 D=1,NU
      YT(A,B,C,D)=Y(A,C,B,D)
 51   CONTINUE
      GO TO 100
 3    CONTINUE
      IAS=IT4(I,M,J,K)
      CALL RDT4(IAS,NU,Y)
      DO 52 A=1,NU
      DO 52 B=1,NU
      DO 52 C=1,NU
      DO 52 D=1,NU
      YT(A,B,C,D)=Y(A,B,C,D)
 52   CONTINUE
      GO TO 100
 4    CONTINUE
      IAS=IT4(M,I,J,K)
      CALL RDT4(IAS,NU,Y)
      DO 53 A=1,NU
      DO 53 B=1,NU
      DO 53 C=1,NU
      DO 53 D=1,NU
      YT(A,B,C,D)=Y(B,A,C,D)
 53   CONTINUE
      GO TO 100
 100  CONTINUE
      call iexit(27)
      RETURN
      END
      SUBROUTINE RDIMJKN(I,J,K,M,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      common/nn/no2,no3,no4,nu2,nu3,nu4,nou,no2u,no3u,nou2,nou3,no2u2
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      call ienter(27)
      CALL IPERM4(I,J,K,M,IPRM)
      GO TO (1,2,3,4)IPRM
 1    CONTINUE
      IAS=IT4(I,J,K,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,2)
      GO TO 100
 2    CONTINUE
      IAS=IT4(I,J,M,K)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,7)
      GO TO 100
 3    CONTINUE
      IAS=IT4(I,M,J,K)
      CALL RDT4(IAS,NU,YT)
      GO TO 100
 4    CONTINUE
      IAS=IT4(M,I,J,K)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,6)
      GO TO 100
 100  CONTINUE
      call iexit(27)
      RETURN
      END
      SUBROUTINE RDIMJKsm(I,J,K,M,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      DATA TWO/2.0D+0/
      call ienter(27)
      NO2=NO*NO
      NO3=NO2*NO
      NO4=NO3*NO
      NU2=NU*NU
      NU3=NU2*NU
      NU4=NU3*NU
      CALL IPERM4(I,J,K,M,IPRM)
      GO TO (1,2,3,4)IPRM
 1    CONTINUE
      IAS=IT4(I,J,K,M)
      CALL RDT4(IAS,NU,Y)
      DO 50 A=1,NU
      DO 50 B=1,NU
      DO 50 C=1,NU
      DO 50 D=1,NU
      YT(A,B,C,D)=Y(A,B,C,D)*TWO-Y(A,D,C,B)-Y(A,B,D,C)-Y(D,B,C,A)
 50   CONTINUE
      GO TO 100
 2    CONTINUE
      IAS=IT4(I,J,M,K)
      CALL RDT4(IAS,NU,Y)
      DO 51 A=1,NU
      DO 51 B=1,NU
      DO 51 C=1,NU
      DO 51 D=1,NU
      YT(A,B,C,D)=Y(A,b,d,c)*TWO-Y(A,D,B,C)-Y(A,B,C,D)-Y(D,B,A,C)
 51   CONTINUE
      GO TO 100
 3    CONTINUE
      IAS=IT4(I,M,J,K)
      CALL RDT4(IAS,NU,Y)
      DO 52 A=1,NU
      DO 52 B=1,NU
      DO 52 C=1,NU
      DO 52 D=1,NU
      YT(A,B,C,D)=Y(a,d,b,c)*TWO-Y(A,B,D,C)-Y(A,C,B,D)-Y(D,A,B,C)
 52   CONTINUE
      GO TO 100
 4    CONTINUE
      IAS=IT4(M,I,J,K)
      CALL RDT4(IAS,NU,Y)
      DO 53 A=1,NU
      DO 53 B=1,NU
      DO 53 C=1,NU
      DO 53 D=1,NU
      YT(A,B,C,D)=Y(d,a,b,c)*TWO-Y(B,A,D,C)-Y(C,A,B,D)-Y(A,D,B,C)
 53   CONTINUE
      GO TO 100
 100  CONTINUE
      call iexit(27)
      RETURN
      END
      SUBROUTINE RDIMKN(K,L,M,N,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      common/nn/no2,no3,no4,nu2,nu3,nu4,nou,no2u,no3u,nou2,nou3,no2u2
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      call ienter(28)
      CALL IPERM(K,L,M,N,IPRM)
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12)IPRM
 1    CONTINUE
      IAS=IT4(K,L,M,N)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,7)
      GO TO 100
 2    CONTINUE
      IAS=IT4(K,L,N,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,2)
      GO TO 100
 3    CONTINUE
      IAS=IT4(K,M,L,N)
      CALL RDT4(IAS,NU,YT)
      GO TO 100
 4    CONTINUE
      IAS=IT4(K,N,L,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,10)
      GO TO 100
 5    CONTINUE
      IAS=IT4(K,M,N,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,1)
      GO TO 100
 6    CONTINUE
      IAS=IT4(K,N,M,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,3)
      GO TO 100
 7    CONTINUE
      IAS=IT4(M,K,L,N)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,6)
      GO TO 100
 8    CONTINUE
      IAS=IT4(N,K,L,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,16)
      GO TO 100
 9    CONTINUE
      IAS=IT4(M,K,N,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,19)
      GO TO 100
 10   CONTINUE
      IAS=IT4(N,K,M,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,21)
      GO TO 100
 11   CONTINUE
      IAS=IT4(M,N,K,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,22)
      GO TO 100
 12   CONTINUE
      IAS=IT4(N,M,K,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,5)
      GO TO 100
 100  CONTINUE
      call iexit(28)
      RETURN
      END
      SUBROUTINE RDIMNL(K,L,M,N,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      common/nn/no2,no3,no4,nu2,nu3,nu4,nou,no2u,no3u,nou2,nou3,no2u2
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      call ienter(29)
      CALL IPERM(K,L,M,N,IPRM)
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12)IPRM
 1    CONTINUE
      IAS=IT4(K,L,M,N)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,3)
      GO TO 100
 2    CONTINUE
      IAS=IT4(K,L,N,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,10)
      GO TO 100
 3    CONTINUE
      IAS=IT4(K,M,L,N)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,1)
      GO TO 100
 4    CONTINUE
      IAS=IT4(K,N,L,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,2)
      GO TO 100
 5    CONTINUE
      IAS=IT4(K,M,N,L)
      CALL RDT4(IAS,NU,YT)
      GO TO 100
 6    CONTINUE
      IAS=IT4(K,N,M,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,7)
      GO TO 100
 7    CONTINUE
      IAS=IT4(M,K,L,N)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,19)
      GO TO 100
 8    CONTINUE
      IAS=IT4(N,K,L,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,23)
      GO TO 100
 9    CONTINUE
      IAS=IT4(M,K,N,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,6)
      GO TO 100
 10   CONTINUE
      IAS=IT4(N,K,M,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,13)
      GO TO 100
 11   CONTINUE
      IAS=IT4(M,N,K,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,12)
      GO TO 100
 12   CONTINUE
      IAS=IT4(N,M,K,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,8)
      GO TO 100
 100  CONTINUE
      call iexit(29)
      RETURN
      END
      SUBROUTINE RDKJML(I,J,K,M,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      call ienter(30)
      NO2=NO*NO
      NO3=NO2*NO
      NO4=NO3*NO
      NU2=NU*NU
      NU3=NU2*NU
      NU4=NU3*NU
      CALL IPERM4(I,J,K,M,IPRM)
      GO TO (1,2,3,4)IPRM
 1    CONTINUE
      IAS=IT4(I,J,K,M)
      CALL RDT4(IAS,NU,Y)
      DO 50 A=1,NU
      DO 50 B=1,NU
      DO 50 C=1,NU
      DO 50 D=1,NU
      YT(B,A,D,C)=Y(A,B,C,D)
 50   CONTINUE
      GO TO 100
 2    CONTINUE
      IAS=IT4(I,J,M,K)
      CALL RDT4(IAS,NU,Y)
      DO 51 A=1,NU
      DO 51 B=1,NU
      DO 51 C=1,NU
      DO 51 D=1,NU
      YT(B,A,C,D)=Y(A,B,C,D)
 51   CONTINUE
      GO TO 100
 3    CONTINUE
      IAS=IT4(I,M,J,K)
      CALL RDT4(IAS,NU,Y)
      DO 52 A=1,NU
      DO 52 B=1,NU
      DO 52 C=1,NU
      DO 52 D=1,NU
      YT(C,A,B,D)=Y(A,B,C,D)
 52   CONTINUE
      GO TO 100
 4    CONTINUE
      IAS=IT4(M,I,J,K)
      CALL RDT4(IAS,NU,Y)
      DO 53 A=1,NU
      DO 53 B=1,NU
      DO 53 C=1,NU
      DO 53 D=1,NU
      YT(C,B,A,D)=Y(A,B,C,D)
 53   CONTINUE
      GO TO 100
 100  CONTINUE
      call iexit(30)
      RETURN
      END
      SUBROUTINE RDLJKM(I,J,K,M,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      call ienter(31)
      NO2=NO*NO
      NO3=NO2*NO
      NO4=NO3*NO
      NU2=NU*NU
      NU3=NU2*NU
      NU4=NU3*NU
      CALL IPERM4(I,J,K,M,IPRM)
      GO TO (1,2,3,4)IPRM
 1    CONTINUE
      IAS=IT4(I,J,K,M)
      CALL RDT4(IAS,NU,Y)
      DO 50 A=1,NU
      DO 50 B=1,NU
      DO 50 C=1,NU
      DO 50 D=1,NU
      YT(C,A,B,D)=Y(A,B,C,D)
 50   CONTINUE
      GO TO 100
 2    CONTINUE
      IAS=IT4(I,J,M,K)
      CALL RDT4(IAS,NU,Y)
      DO 51 A=1,NU
      DO 51 B=1,NU
      DO 51 C=1,NU
      DO 51 D=1,NU
      YT(D,A,B,C)=Y(A,B,C,D)
 51   CONTINUE
      GO TO 100
 3    CONTINUE
      IAS=IT4(I,M,J,K)
      CALL RDT4(IAS,NU,Y)
      DO 52 A=1,NU
      DO 52 B=1,NU
      DO 52 C=1,NU
      DO 52 D=1,NU
      YT(D,A,C,B)=Y(A,B,C,D)
 52   CONTINUE
      GO TO 100
 4    CONTINUE
      IAS=IT4(M,I,J,K)
      CALL RDT4(IAS,NU,Y)
      DO 53 A=1,NU
      DO 53 B=1,NU
      DO 53 C=1,NU
      DO 53 D=1,NU
      YT(D,B,C,A)=Y(A,B,C,D)
 53   CONTINUE
      GO TO 100
 100  CONTINUE
      call iexit(31)
      RETURN
      END
      SUBROUTINE RDMIJK(I,J,K,M,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      call ienter(32)
      NO2=NO*NO
      NO3=NO2*NO
      NO4=NO3*NO
      NU2=NU*NU
      NU3=NU2*NU
      NU4=NU3*NU
      CALL IPERM4(I,J,K,M,IPRM)
      GO TO (1,2,3,4)IPRM
 1    CONTINUE
      IAS=IT4(I,J,K,M)
      CALL RDT4(IAS,NU,Y)
      DO 50 A=1,NU
      DO 50 B=1,NU
      DO 50 C=1,NU
      DO 50 D=1,NU
      YT(A,B,C,D)=Y(B,C,D,A)
 50   CONTINUE
      GO TO 100
 2    CONTINUE
      IAS=IT4(I,J,M,K)
      CALL RDT4(IAS,NU,Y)
      DO 51 A=1,NU
      DO 51 B=1,NU
      DO 51 C=1,NU
      DO 51 D=1,NU
      YT(A,B,C,D)=Y(b,c,a,d)
 51   CONTINUE
      GO TO 100
 3    CONTINUE
      IAS=IT4(I,M,J,K)
      CALL RDT4(IAS,NU,Y)
      DO 52 A=1,NU
      DO 52 B=1,NU
      DO 52 C=1,NU
      DO 52 D=1,NU
      YT(A,B,C,D)=Y(b,a,c,D)
 52   CONTINUE
      GO TO 100
 4    CONTINUE
      IAS=IT4(M,I,J,K)
      CALL RDT4(IAS,NU,Y)
      DO 53 A=1,NU
      DO 53 B=1,NU
      DO 53 C=1,NU
      DO 53 D=1,NU
      YT(A,B,C,D)=Y(A,B,C,D)
 53   CONTINUE
      GO TO 100
 100  CONTINUE
      call iexit(32)
      RETURN
      END
      SUBROUTINE RDMIJKS(I,J,K,M,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      common/nn/no2,no3,no4,nu2,nu3,nu4,nou,no2u,no3u,nou2,nou3,no2u2
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      call ienter(32)
      CALL IPERM4(I,J,K,M,IPRM)
      GO TO (1,2,3,4)IPRM
 1    CONTINUE
      IAS=IT4(I,J,K,M)
      CALL RDT4(IAS,NU,YT)
 50   CONTINUE
      GO TO 100
 2    CONTINUE
      IAS=IT4(I,J,M,K)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,1)
 51   CONTINUE
      GO TO 100
 3    CONTINUE
      IAS=IT4(I,M,J,K)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,3)
 52   CONTINUE
      GO TO 100
 4    CONTINUE
      IAS=IT4(M,I,J,K)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,21)
      GO TO 100
 100  CONTINUE
      call iexit(32)
      RETURN
      END
      SUBROUTINE RDMIJKSM(I,J,K,M,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      DATA TWO/2.0D+0/
      call ienter(32)
      NO2=NO*NO
      NO3=NO2*NO
      NO4=NO3*NO
      NU2=NU*NU
      NU3=NU2*NU
      NU4=NU3*NU
      CALL IPERM4(I,J,K,M,IPRM)
      GO TO (1,2,3,4)IPRM
 1    CONTINUE
      IAS=IT4(I,J,K,M)
      CALL RDT4(IAS,NU,Y)
      DO 50 A=1,NU
      DO 50 B=1,NU
      DO 50 C=1,NU
      DO 50 D=1,NU
      YT(A,B,C,D)=Y(A,B,C,D)*TWO-Y(D,B,C,A)-Y(A,D,C,B)-Y(A,B,D,C)
 50   CONTINUE
      GO TO 100
 2    CONTINUE
      IAS=IT4(I,J,M,K)
      CALL RDT4(IAS,NU,Y)
      DO 51 A=1,NU
      DO 51 B=1,NU
      DO 51 C=1,NU
      DO 51 D=1,NU
      YT(A,B,C,D)=Y(A,B,D,C)*TWO-Y(D,B,A,C)-Y(A,D,B,C)-Y(A,B,C,D)
 51   CONTINUE
      GO TO 100
 3    CONTINUE
      IAS=IT4(I,M,J,K)
      CALL RDT4(IAS,NU,Y)
      DO 52 A=1,NU
      DO 52 B=1,NU
      DO 52 C=1,NU
      DO 52 D=1,NU
      YT(A,B,C,D)=Y(A,D,B,C)*TWO-Y(D,A,B,C)-Y(A,B,D,C)-Y(A,C,B,D)
 52   CONTINUE
      GO TO 100
 4    CONTINUE
      IAS=IT4(M,I,J,K)
      CALL RDT4(IAS,NU,Y)
      DO 53 A=1,NU
      DO 53 B=1,NU
      DO 53 C=1,NU
      DO 53 D=1,NU
      YT(A,B,C,D)=Y(D,A,B,C)*TWO-Y(A,D,B,C)-Y(B,A,D,C)-Y(C,A,B,D)
 53   CONTINUE
      GO TO 100
 100  CONTINUE
      call iexit(32)
      RETURN
      END
      SUBROUTINE RDMJKN(K,L,M,N,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      common/nn/no2,no3,no4,nu2,nu3,nu4,nou,no2u,no3u,nou2,nou3,no2u2
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      call ienter(33)
      CALL IPERM(K,L,M,N,IPRM)
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12)IPRM
 1    CONTINUE
      IAS=IT4(K,L,M,N)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,12)
      GO TO 100
 2    CONTINUE
      IAS=IT4(K,L,N,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,20)
      GO TO 100
 3    CONTINUE
      IAS=IT4(K,M,L,N)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,6)
      GO TO 100
 4    CONTINUE
      IAS=IT4(K,N,L,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,15)
      GO TO 100
 5    CONTINUE
      IAS=IT4(K,M,N,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,19)
      GO TO 100
 6    CONTINUE
      IAS=IT4(K,N,M,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,22)
      GO TO 100
 7    CONTINUE
      IAS=IT4(M,K,L,N)
      CALL RDT4(IAS,NU,YT)
      GO TO 100
 8    CONTINUE
      IAS=IT4(N,K,L,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,11)
      GO TO 100
 9    CONTINUE
      IAS=IT4(M,K,N,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,1)
      GO TO 100
 10   CONTINUE
      IAS=IT4(N,K,M,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,5)
      GO TO 100
 11   CONTINUE
      IAS=IT4(M,N,K,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,3)
      GO TO 100
 12   CONTINUE
      IAS=IT4(N,M,K,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,21)
      GO TO 100
 100  CONTINUE
      call iexit(33)
      RETURN
      END
      SUBROUTINE RDMJNL(K,L,M,N,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      common/nn/no2,no3,no4,nu2,nu3,nu4,nou,no2u,no3u,nou2,nou3,no2u2
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      call ienter(34)
      CALL IPERM(K,L,M,N,IPRM)
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12)IPRM
 1    CONTINUE
      IAS=IT4(K,L,M,N)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,22)
      GO TO 100
 2    CONTINUE
      IAS=IT4(K,L,N,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,15)
      GO TO 100
 3    CONTINUE
      IAS=IT4(K,M,L,N)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,19)
      GO TO 100
 4    CONTINUE
      IAS=IT4(K,N,L,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,20)
      GO TO 100
 5    CONTINUE
      IAS=IT4(K,M,N,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,6)
      GO TO 100
 6    CONTINUE
      IAS=IT4(K,N,M,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,12)
      GO TO 100
 7    CONTINUE
      IAS=IT4(M,K,L,N)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,1)
      GO TO 100
 8    CONTINUE
      IAS=IT4(N,K,L,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,14)
      GO TO 100
 9    CONTINUE
      IAS=IT4(M,K,N,L)
      CALL RDT4(IAS,NU,YT)
      GO TO 100
 10   CONTINUE
      IAS=IT4(N,K,M,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,8)
      GO TO 100
 11   CONTINUE
      IAS=IT4(M,N,K,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,7)
      GO TO 100
 12   CONTINUE
      IAS=IT4(N,M,K,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,13)
      GO TO 100
 100  CONTINUE
      call iexit(34)
      RETURN
      END
      SUBROUTINE RDMNKL(K,L,M,N,NO,NU,Y,YT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      common/nn/no2,no3,no4,nu2,nu3,nu4,nou,no2u,no3u,nou2,nou3,no2u2
      DIMENSION Y(NU,NU,NU,NU),YT(NU,NU,NU,NU)
      call ienter(35)
      CALL IPERM(K,L,M,N,IPRM)
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12)IPRM
 1    CONTINUE
      IAS=IT4(K,L,M,N)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,9)
      GO TO 100
 2    CONTINUE
      IAS=IT4(K,L,N,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,4)
      GO TO 100
 3    CONTINUE
      IAS=IT4(K,M,L,N)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,23)
      GO TO 100
 4    CONTINUE
      IAS=IT4(K,N,L,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,14)
      GO TO 100
 5    CONTINUE
      IAS=IT4(K,M,N,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,13)
      GO TO 100
 6    CONTINUE
      IAS=IT4(K,N,M,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,8)
      GO TO 100
 7    CONTINUE
      IAS=IT4(M,K,L,N)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,2)
      GO TO 100
 8    CONTINUE
      IAS=IT4(N,K,L,M)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,20)
      GO TO 100
 9    CONTINUE
      IAS=IT4(M,K,N,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,7)
      GO TO 100
 10   CONTINUE
      IAS=IT4(N,K,M,L)
      CALL RDT4(IAS,NU,Yt)
      call mtrans(yt,nu,12)
      GO TO 100
 11   CONTINUE
      IAS=IT4(M,N,K,L)
      CALL RDT4(IAS,NU,YT)
      GO TO 100
 12   CONTINUE
      IAS=IT4(N,M,K,L)
      CALL RDT4(IAS,NU,YT)
      call mtrans(yt,nu,6)
 100  CONTINUE
      call iexit(35)
      RETURN
      END
      SUBROUTINE RDOV4(INO,NO,NU,TI,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION TI(NU,NU,NU),V(NU,NU,NU,NU)
      call ienter(112)
      NNO=NO
      NNU=NU
      INN=INO
      IF(INO.EQ.1.OR.INO.EQ.11)THEN
      NNO=NU
      NNU=NO
      ENDIF
      NLAST=4*NNO+NNU
      IF (INO.GE.10) THEN
      INN=INO-10
      NLAST=6*NNO+3*NNU+7
      ENDIF
      INO1=1-INN
      DO 1 A=1,NU
      IASV=NLAST+INO1*NNO+A
      READ(NALL4,REC=IASV)TI
      DO 1 B=1,NU
      DO 1 C=1,NU
      DO 1 D=1,NU
      V(A,B,C,D)=TI(B,C,D)
 1    CONTINUE
      NU4=NU*NU*NU*NU
      if(ino.eq.0)then
      endif
      call iexit(112)
      RETURN
      END
      SUBROUTINE RDOV4A(INO,NO,NU,TI,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION TI(NU,NU,NU),V(NU,NU,NU,NU)
      call ienter(112)
      NNO=NO
      NNU=NU
      INN=INO
      NLAST=0
      NALL4=45
      IF(INO.EQ.1)THEN
      NNO=NU
      NNU=NO
      ENDIF
      INO1=1-INN
      DO 1 A=1,NU
      IASV=NLAST+INO1*NNO+A
      READ(NALL4,REC=IASV)TI
      DO 1 B=1,NU
      DO 1 C=1,NU
      DO 1 D=1,NU
      V(A,B,C,D)=TI(B,C,D)
 1    CONTINUE
      call iexit(112)
      RETURN
      END
      SUBROUTINE RDT32(K,L,NO,NU,TI,T3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      DIMENSION TI(NU,NU,NU),T3(NO,NU,NU,NU)
      DATA TWO/2.0D+0/
      NU3=NU*NU*NU
      DO 100 I=1,NO
      IF (I.EQ.K.AND.I.EQ.L)THEN
      CALL ZEROMA(TI,1,NU3)
      GO TO 90
      ENDIF
      CALL RDVT3ONW(I,K,L,NU,TI)
 90   CONTINUE
      DO 110 A=1,NU
      DO 110 B=1,NU
      DO 110 C=1,NU
      T3(I,A,B,C)=TWO*TI(A,B,C)-TI(B,A,C)-TI(C,B,A)
 110  CONTINUE
 100  CONTINUE
      RETURN
      END
      SUBROUTINE RDT3A(NO,NU,TI,T3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      DIMENSION TI(NU,NU,NU),T3(NO,NU,NU,NU,NO,NO)
      NO2=NO*NO
      DO 1 I=1,NO
      DO 1 J=1,NO
      DO 1 K=1,NO
      IF (I.EQ.J.AND.J.EQ.K) THEN
      CALL ZEROMA(TI,1,NU3)
      GOTO 9
      ENDIF
      IAS=NO2*(I-1)+NO*(J-1)+K
      CALL RDVT3(IAS,NU,TI)
 9    CONTINUE
      DO 11 A=1,NU
      DO 11 B=1,NU
      DO 11 C=1,NU
      T3(I,A,B,C,J,K)=TI(A,B,C)
 11   CONTINUE
 1    CONTINUE
      RETURN
      END

      SUBROUTINE RDT4(IASV,NU,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NEWT4/NT4,NO4,LT4,NALL4,LL4
      DIMENSION V(NU,NU,NU,NU)
      call ienter(131)
      READ(NO4,REC=IASV)V
      call iexit(131)
      RETURN
      END
      SUBROUTINE RDT4N(IASV,NU,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NEWT4/NT4,NO4,LT4,NALL4,LL4
      DIMENSION V(NU,NU,NU,NU)
      READ(NT4,REC=IASV)V
      RETURN
      END
      SUBROUTINE RDV43(A,INO,NO,NU,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION V(NU,NU,NU)
      INN=INO
      NLAST=4*NO+NU
      IF (INO.EQ.10) THEN
      INN=INO-10
      NLAST=6*NO+3*NU+7
      ENDIF
      IASV=NLAST+NO+A
      READ(NALL4,REC=IASV)V
      CALL TRANT3(V,NU,5)
      RETURN
      END
      SUBROUTINE RDVA(IAS,NU,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A
      DIMENSION V(NU,NU,NU,NU)
      call ienter(36)
      DO 10 A=1,NU
      IASV=IAS+A
      CALL RDVT3(IASV,NU,V(1,1,1,A))
   10 CONTINUE
      call iexit(36)
      RETURN
      END



      SUBROUTINE RDVEM1(IN,NO,NU,TI,VEM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION TI(NU,NU,NU),VEM(NO,NU,NU,NU)
      DATA ZERO/0.0D+0/
      call ienter(110)
      INO=IABS(IN)
      NNO=NO
      NNU=NU
      IF (INO.EQ.1.OR.INO.EQ.11.OR.INO.EQ.21)THEN
      NNO=NU
      NNU=NO
      ENDIF
      IF(INO.GE.20) THEN 
      INO=INO-20
      NLAST=9*NNO+4*NNU+7
      ELSE 
      IF (INO.GE.10) THEN
      INO=INO-10
      NLAST=5*NNO+2*NNU+7
      ELSE
      NLAST=3*NNO
      ENDIF
      ENDIF
      DO 1 I=1,NO
      IASV=NLAST+NU*INO+I
      READ(NALL4,REC=IASV)TI
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 C=1,NU
      X=ZERO
      IF (IN.LT.0)X=VEM(I,C,A,B)
      VEM(I,C,A,B)=X+TI(A,B,C)
 1    CONTINUE
      NOU3=NO*NU*NU*NU
      call iexit(110)
      RETURN
      END
      SUBROUTINE RDVEM1A(IN,NO,NU,TI,VEM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      DIMENSION TI(NU,NU,NU),VEM(NO,NU,NU,NU)
      DATA ZERO/0.0D+0/
      call ienter(110)
      NALL4=45
      nu3=nu*nu*Nu
      nou3=no*nu3
      NNO=NO
      NNU=NU
      IF (IN.EQ.1.OR.IN.EQ.11)THEN
      NNO=NU
      NNU=NO
      ENDIF
      I11=0
      IF(IN.GE.10)I11=NO+NU
      INO=1
      IF(IN.EQ.1.OR.IN.EQ.11)INO=0
      NLAST=3*NNO+NNU+I11
      DO 1 I=1,NO
      IASV=NLAST+NU*INO+I
      READ(NALL4,REC=IASV)TI
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 C=1,NU
      VEM(I,C,A,B)=TI(A,B,C)
 1    CONTINUE
      call iexit(110)
      return
      END
      SUBROUTINE RDVEM2(IN,NO,NU,TI,VEM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION TI(NU,NU,NU),VEM(NU,NO,NU,NU)
      DATA ZERO/0.0D+0/
      call ienter(110)
      INO=IABS(IN)
      NNO=NO
      NNU=NU
      IF (INO.EQ.1.OR.INO.EQ.11.OR.INO.EQ.21)THEN
      NNO=NU
      NNU=NO
      ENDIF
      IF(INO.GE.20) THEN 
      INO=INO-20
      NLAST=9*NNO+4*NNU+7
      ELSE 
      IF (INO.GE.10) THEN
      INO=INO-10
      NLAST=5*NNO+2*NNU+7
      ELSE
      NLAST=3*NNO
      ENDIF
      ENDIF
      DO 1 I=1,NO
      IASV=NLAST+NU*INO+I
      READ(NALL4,REC=IASV)TI
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 C=1,NU
      X=ZERO
      IF(IN.LT.0)X=VEM(C,I,A,B)
      VEM(C,I,A,B)=X+TI(A,B,C)
 1    CONTINUE
      NOU3=NO*NU*NU*NU
      call iexit(110)
      RETURN
      END
      SUBROUTINE RDVEM2A(IN,NO,NU,TI,VEM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      DIMENSION TI(NU,NU,NU),VEM(NU,NO,NU,NU)
      DATA ZERO/0.0D+0/
      call ienter(110)
      NALL4=45
      NNO=NO
      NNU=NU
      IF (IN.EQ.1.OR.IN.EQ.11)THEN
      NNO=NU
      NNU=NO
      ENDIF
      I11=0
      IF(IN.GE.10)I11=NO+NU
      INO=1
      IF(IN.EQ.1.OR.IN.EQ.11)INO=0
      NLAST=3*NNO+NNU+I11
      DO 1 I=1,NO
      IASV=NLAST+NU*INO+I
      READ(NALL4,REC=IASV)TI
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 C=1,NU
      VEM(C,I,A,B)=TI(A,B,C)
 1    CONTINUE
      call iexit(110)
      return
      END
      SUBROUTINE RDVEM3(IN,NO,NU,TI,VEM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION TI(NU,NU,NU),VEM(NU,NU,NO,NU)
      DATA ZERO/0.0D+0/
      call ienter(110)
      INO=IABS(IN)
      NNO=NO
      NNU=NU
      IF (INO.EQ.1.OR.INO.EQ.11.OR.INO.EQ.21)THEN
      NNO=NU
      NNU=NO
      ENDIF
      IF(INO.GE.20) THEN 
      INO=INO-20
      NLAST=9*NNO+4*NNU+7
      ELSE 
      IF (INO.GE.10) THEN
      INO=INO-10
      NLAST=5*NNO+2*NNU+7
      ELSE
      NLAST=3*NNO
      ENDIF
      ENDIF
      DO 1 I=1,NO
      IASV=NLAST+NU*INO+I
      READ(NALL4,REC=IASV)TI
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 C=1,NU
      X=ZERO
      IF (IN.LT.0)X=VEM(A,B,I,C)
      VEM(A,B,I,C)=X+TI(A,B,C)
 1    CONTINUE
      NOU3=NO*NU*NU*NU
      call iexit(110)
      RETURN
      END
      SUBROUTINE RDVEM4(IN,NO,NU,TI,VEM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION TI(NU,NU,NU),VEM(NU,NU,NU,NO)
      DATA ZERO/0.0D+0/
      call ienter(110)
      INO=IABS(IN)
      NNO=NO
      NNU=NU
      IF (INO.EQ.1.OR.INO.EQ.11.OR.INO.EQ.21)THEN
      NNO=NU
      NNU=NO
      ENDIF
      IF(INO.GE.20) THEN 
      INO=INO-20
      NLAST=9*NNO+4*NNU+7
      ELSE 
      IF (INO.GE.10) THEN
      INO=INO-10
      NLAST=5*NNO+2*NNU+7
      ELSE
      NLAST=3*NNO
      ENDIF
      ENDIF
      DO 1 I=1,NO
      IASV=NLAST+NU*INO+I
      READ(NALL4,REC=IASV)TI
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 C=1,NU
      X=ZERO
      IF(IN.LT.0)X=VEM(A,B,C,I)
      VEM(A,B,C,I)=X+TI(A,B,C)
 1    CONTINUE
      NOU3=NO*NU*NU*NU
      call iexit(110)
      RETURN
      END
      SUBROUTINE RDVEM4A(IN,NO,NU,TI,VEM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      DIMENSION TI(NU,NU,NU),VEM(NU,NU,NU,NO)
      DATA ZERO/0.0D+0/
      call ienter(110)
      NALL4=45
      nu3=nu*nu*Nu
      nou3=no*nu3
      NNO=NO
      NNU=NU
      IF (IN.EQ.1.OR.IN.EQ.11)THEN
      NNO=NU
      NNU=NO
      ENDIF
      I11=0
      IF(IN.GE.10)I11=NO+NU
      INO=1
      IF(IN.EQ.1.OR.IN.EQ.11)INO=0
      NLAST=3*NNO+NNU+I11
      DO 1 I=1,NO
      IASV=NLAST+NU*INO+I
      READ(NALL4,REC=IASV)TI
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 C=1,NU
      VEM(A,B,C,I)=TI(A,B,C)
 1    CONTINUE
      call iexit(110)
      return
      END

      SUBROUTINE RDVT3(IASV,NU,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NEWTAPE/NT2T4,NT4INT
      DIMENSION V(NU,NU,NU)
      call ienter(37)
      READ(NT4INT,REC=IASV)V
      call iexit(37)
      RETURN
      END
      SUBROUTINE RDVT3I(IASV,NU,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/NEWT3/NT3INT,LT3INT
      DIMENSION V(NU,NU,NU)
      call ienter(38)
      READ(NT3INT,REC=IASV)V
      call iexit(38)
      RETURN
      END
      SUBROUTINE RDVT3IMN(NO,NU,TI,VT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A
      COMMON/NEWT3/NT3INT,LT3INT
      DIMENSION TI(NO,NO,NO),VT(NU,NO,NO,NO,NO,NO)
      NO2=NO*NO
      NO3=NO2*NO
      DO 10 A=1,NU
      DO 10 M=1,NO
      DO 10 N=1,NO
      IASV=2*NO3+NO2*(A-1)+NO*(M-1)+N
      READ(NT3INT,REC=IASV)TI
      DO 15 I=1,NO
      DO 15 J=1,NO
      DO 15 K=1,NO
      VT(A,M,N,I,J,K)=TI(I,J,K)
 15   CONTINUE
 10   CONTINUE
      RETURN
      END
      SUBROUTINE RDVT3N(IASV,NU,V)
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      COMMON/NT3T3/NO3,NT3,LT3
      DIMENSION V(NU,NU,NU)
      call  ienter(39)
      READ(NT3,REC=IASV)V
      call iexit(39)
      RETURN
      END
      SUBROUTINE RDVT3O(IASV,NU,V)
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      COMMON/NT3T3/NO3,NT3,LT3
      DIMENSION V(NU,NU,NU)
      call ienter(40)
      READ(NO3,REC=IASV)V
      call iexit(40)
      RETURN
      END
      SUBROUTINE RDVT3ONW(I,J,K,NU,V)
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      common/nn/no2,nno3,no4,nu2,nu3,nu4,nou,no2u,no3u,nou2,nou3,no2u2
      COMMON/NT3T3/NO3,NT3,LT3
      DIMENSION V(NU,NU,NU)
      call ienter(41)
 489  format('i,j,k:',i3,5x,3i6,3x,3i3)
      if(i.eq.j.and.j.eq.k)then
         call zeroma(v,1,nu3)
         return
      endif
      IF (I.GT.J)THEN
      IF(J.GT.K) THEN
         IAS=IT3(I,J,K)
      READ(NO3,REC=IAS)V
      ELSE
      IF(I.GT.K)THEN
      IAS=IT3(I,K,J)
      READ(NO3,REC=IAS)V
      CALL TRANt3(V,NU,1)
      ELSE
      IAS=IT3(K,I,J)
      READ(NO3,REC=IAS)V
      CALL TRANt3(V,NU,5)
      ENDIF
      ENDIF
      ELSE
      IF (I.GT.K)THEN
      IAS=IT3(J,I,K)
      READ(NO3,REC=IAS)V
      CALL TRANt3(V,NU,2)
      ELSE
      IF(J.GT.K) THEN
      IAS=IT3(J,K,I)
      READ(NO3,REC=IAS)V
      CALL TRANt3(V,NU,4)
      ELSE
      IAS=IT3(K,J,I)
      READ(NO3,REC=IAS)V
      CALL TRANt3(V,NU,3)
      ENDIF
      ENDIF
      ENDIF
      call iexit(41)
      RETURN
      END

      SUBROUTINE RDVT3NNW(I,J,K,NU,V)
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      COMMON/NT3T3/NT3,NO3,LT3
      DIMENSION V(NU,NU,NU)
      call ienter(41)
 489  format('i,j,k:',i3,5x,3i6,3x,3i3)
      IF (I.GT.J)THEN
      IF(J.GT.K) THEN
         IAS=IT3(I,J,K)
      READ(NO3,REC=IAS)V
      ELSE
      IF(I.GT.K)THEN
      IAS=IT3(I,K,J)
      READ(NO3,REC=IAS)V
      CALL TRANt3(V,NU,1)
      ELSE
      IAS=IT3(K,I,J)
      READ(NO3,REC=IAS)V
      CALL TRANt3(V,NU,5)
      ENDIF
      ENDIF
      ELSE
      IF (I.GT.K)THEN
      IAS=IT3(J,I,K)
      READ(NO3,REC=IAS)V
      CALL TRANt3(V,NU,2)
      ELSE
      IF(J.GT.K) THEN
      IAS=IT3(J,K,I)
      READ(NO3,REC=IAS)V
      CALL TRANt3(V,NU,4)
      ELSE
      IAS=IT3(K,J,I)
      READ(NO3,REC=IAS)V
      CALL TRANt3(V,NU,3)
      ENDIF
      ENDIF
      ENDIF
      call iexit(41)
      RETURN
      END
      SUBROUTINE WRVT3NNW(I,J,K,NU,V)
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      COMMON/NT3T3/NT3,NO3,LT3
      DIMENSION V(NU,NU,NU)
      call ienter(41)
      nu3=nu*nu*nu
 489  format('i,j,k:',i3,5x,3i6,3x,3i3)
      IF (I.GT.J)THEN
         IF(J.GT.K) THEN
            IAS=IT3(I,J,K)
            write(NO3,REC=IAS)V
         ELSE
            IF(I.GT.K)THEN
               IAS=IT3(I,K,J)
               CALL TRANt3(V,NU,1)
               write(NO3,REC=IAS)V
            ELSE
               IAS=IT3(K,I,J)
               CALL TRANt3(V,NU,5)
               write(NO3,REC=IAS)V
            ENDIF
         ENDIF
      ELSE
         IF (I.GT.K)THEN
            IAS=IT3(J,I,K)
            CALL TRANt3(V,NU,2)
            write(NO3,REC=IAS)V
         ELSE
            IF(J.GT.K) THEN
               IAS=IT3(J,K,I)
               CALL TRANt3(V,NU,4)
               write(NO3,REC=IAS)V
            ELSE
               IAS=IT3(K,J,I)
               CALL TRANt3(V,NU,3)
               write(NO3,REC=IAS)V
            ENDIF
         ENDIF
      ENDIF
      call iexit(41)
      RETURN
      END




      SUBROUTINE RO2HHP(INO,NO,NU,TI,T2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON /NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION TI(NU,NU,NO),T2(NO,NO,NU,NU)
      call ienter(111)
      NLAST=0
      INN=INO
      NOU2=NO*NU*NU
      IF(INO.GE.10) THEN
      INN=INO-10
      NLAST=6*NO+4*NU+7
      ENDIF
      DO 1 I=1,NO
      IASV=NLAST+I+INN*NO
      READ(NALL4,REC=IASV)TI
      DO 1 J=1,NO
      DO 1 A=1,NU
      DO 1 B=1,NU
      T2(I,J,A,B)=TI(A,B,J)
 1    CONTINUE
      NO2U2=NO*NO*NU*NU
      call iexit(111)
      RETURN
      END
      SUBROUTINE RO2HPH(INO,NO,NU,TI,T2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON /NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION TI(NU,NU,NO),T2(NO,NU,NO,NU)
      call ienter(111)
      NLAST=0
      NOU2=NO*NU*NU
      INN=INO
      IF(INO.GE.10) THEN
      INN=INO-10
      NLAST=6*NO+4*NU+7
      ENDIF
      DO 1 I=1,NO
      IASV=NLAST+I+INN*NO
      READ(NALL4,REC=IASV)TI
      DO 1 J=1,NO
      DO 1 A=1,NU
      DO 1 B=1,NU
      T2(I,A,J,B)=TI(A,B,J)
 1    CONTINUE
      NO2U2=NO*NO*NU*NU
      call iexit(111)
      RETURN
      END
      SUBROUTINE RO2HPP(INO,NO,NU,TI,T2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON /NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION TI(NU,NU,NO),T2(NO,NU,NU,NO)
      call ienter(111)
      NLAST=0
      INN=INO
      IF(INO.GE.10) THEN
      INN=INO-10
      NLAST=6*NO+4*NU+7
      ENDIF
      DO 1 I=1,NO
      IASV=NLAST+I+INN*NO
      READ(NALL4,REC=IASV)TI
      DO 1 J=1,NO
      DO 1 A=1,NU
      DO 1 B=1,NU
      T2(I,A,B,J)=TI(A,B,J)
 1    CONTINUE
      NO2U2=NO*NO*NU*NU
      RETURN
      END
      SUBROUTINE RDINTVO(IP,NO,NU,TI,T2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      DIMENSION TI(NU,NU,NO),T2(NO,NU,NU,NO)
      NTVOE=60
      NLAST=NO*IP
      DO 1 I=1,NO
      IASV=NLAST+I
      READ(NTVOE,REC=IASV)TI
      DO 1 J=1,NO
      DO 1 A=1,NU
      DO 1 B=1,NU
      T2(I,A,B,J)=TI(A,B,J)
 1    CONTINUE
      RETURN
      END
      SUBROUTINE RDINTPP(IP,N,TI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      DIMENSION TI(N,N)
      NTVOE=59
      IASV=IP+1
      READ(NTVOE,REC=IASV)TI
      RETURN
      END
      SUBROUTINE WRINTPP(IP,N,TI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      DIMENSION TI(N,N)
      NTVOE=59
      IASV=IP+1
      WRITE(NTVOE,REC=IASV)TI
      RETURN
      END
      SUBROUTINE RINTHH(IP,N,TI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      DIMENSION TI(N,N)
      NHH=58
      IF(IP.EQ.1)THEN
         n2=n*n
         IR=1
         call rdd(nhh,ir,n2,ti)
         READ(NHH,REC=IR)TI
      ELSE
         DO 10 I=1,N
            IR=I+1
         call rdd(nhh,ir,n,ti(1,i))
 10      CONTINUE
      ENDIF
      RETURN
      END
      SUBROUTINE VRINTHH(IP,N,TI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      DIMENSION TI(N,N)
      NPP=58
      IF(IP.EQ.1)THEN
         n2=n*n
         IR=1
            call wrt(npp,ir,n2,ti)
      ELSE
         DO 10 I=1,N
            IR=I+1
            call wrt(npp,ir,n,ti(1,i))
 10      CONTINUE
      ENDIF
      RETURN
      END
      SUBROUTINE RINTPP(IP,N,TI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      DIMENSION TI(N,N)
      NPP=59
      n2=n*n
      IF(IP.EQ.1)THEN
         IR=1
         call rdd(npp,ir,n2,ti)
      ELSE
         DO 10 I=1,N
            IR=I+1
         call rdd(npp,ir,n,ti(1,i))
 10      CONTINUE
      ENDIF
      RETURN
      END
      SUBROUTINE VRINTPP(IP,N,TI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      DIMENSION TI(N,N)
      NPP=59
      n2=n*n
      IF(IP.EQ.1)THEN
         IR=1
            call wrt(npp,ir,n2,ti)
         WRITE(NPP,REC=IR)TI
      ELSE
         DO 10 I=1,N
            IR=I+1
            call wrt(npp,ir,n,ti(1,i))
 10      CONTINUE
      ENDIF
      RETURN
      END
      SUBROUTINE VRINTPPtmp(IP,N,TI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      DIMENSION TI(N,N)
      NPP=57
      n2=n*n
      IF(IP.EQ.1)THEN
         IR=1
            call wrt(npp,ir,n2,ti)
         WRITE(NPP,REC=IR)TI
      ELSE
         DO 10 I=1,N
            IR=I+1
            call wrt(npp,ir,n,ti(1,i))
 10      CONTINUE
      ENDIF
      RETURN
      END
      subroutine wrt(npp,ir,n,t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension t(n)
      write(npp,rec=ir)t
      return
      end
      subroutine rdd(npp,ir,n,t)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension t(n)
      read(npp,rec=ir)t
      return
      end
      SUBROUTINE WRINTVO(IP,NO,NU,TI,T2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      DIMENSION TI(NU,NU,NO),T2(NO,NU,NU,NO)
      NTVOE=60
      NLAST=NO*IP
      DO 1 I=1,NO
      IASV=NLAST+I
      DO 2 J=1,NO
      DO 2 A=1,NU
      DO 2 B=1,NU
      TI(A,B,J)=T2(I,A,B,J)
 2    CONTINUE
      WRITE(NTVOE,REC=IASV)TI
 1    CONTINUE
      RETURN
      END
      SUBROUTINE VRINTVO(IP,NO,NU,TI,T2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      DIMENSION TI(NU,NU,NO),T2(NO,NU,NU,NO)
      NTVOE=60
      NLAST=NO*(IP-1)
      DO 1 I=1,NO
      IASV=NLAST+I
      DO 2 J=1,NO
      DO 2 A=1,NU
      DO 2 B=1,NU
      TI(A,B,J)=T2(I,A,B,J)
 2    CONTINUE
      WRITE(NTVOE,REC=IASV)TI
 1    CONTINUE
      RETURN
      END
      SUBROUTINE RINTVO(IP,NO,NU,TI,T2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      DIMENSION TI(NU,NU,NO),T2(NO,NU,NU,NO)
      data zero/0.0d+0/
      NTVOE=60
      IVO=IABS(IP)
      NLAST=NO*(IVO-1)
      DO 1 I=1,NO
      IASV=NLAST+I
      READ(NTVOE,REC=IASV)TI
      DO 1 J=1,NO
      DO 1 A=1,NU
      DO 1 B=1,NU
         X=ZERO
         IF (IP.LT.0)X=T2(I,A,B,J)
      T2(I,A,B,J)=X+TI(A,B,J)
 1    CONTINUE
      call iexit(111)
      RETURN
      END
      SUBROUTINE RO2PHP(INO,NO,NU,TI,T2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON /NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION TI(NU,NU,NO),T2(NU,NO,NU,NO)
      call ienter(111)
      NLAST=0
      INN=INO
      IF(INO.GE.10) THEN
      INN=INO-10
      NLAST=6*NO+4*NU+7
      ENDIF
      DO 1 I=1,NO
      IASV=NLAST+I+INN*NO
      READ(NALL4,REC=IASV)TI
      DO 1 J=1,NO
      DO 1 A=1,NU
      DO 1 B=1,NU
      T2(A,I,B,J)=TI(A,B,J)
 1    CONTINUE
      NO2U2=NO*NO*NU*NU
      call iexit(111)
      RETURN
      END
      SUBROUTINE RO2PPH(ION,NO,NU,TI,T2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON /NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION TI(NU,NU,NO),T2(NU,NU,NO,NO)
      call ienter(111)
      INO=IABS(ION)
      NTP=NALL4
      NLAST=0
      INN=INO
      IF(INO.GE.10) THEN
      INN=INO-10
      NLAST=6*NO+4*NU+7
      ENDIF
      IF (ION.LT.0) THEN
      NTP=45
      NLAST=NO+NU
      INN=INO-1
      ENDIF
      DO 1 I=1,NO
      IASV=NLAST+I+INN*NO
      READ(NTP,REC=IASV)TI
      DO 1 J=1,NO
      DO 1 A=1,NU
      DO 1 B=1,NU
      T2(A,B,I,J)=TI(A,B,J)
 1    CONTINUE
      NO2U2=NO*NO*NU*NU
      call iexit(111)
      RETURN
      END
      SUBROUTINE ROVOEP(IS,NO,NU,TI,VOE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      COMMON/NEWT4/NT4,NOT4,LT4,NALL4,LL4
      DIMENSION TI(NU,NU,NO),VOE(NO,NU,NU,NO)
      NLAST=10*NO+5*NU+7
      NLAST=NLAST+(IS-1)*NO
      nou2=no*nu*nu
      no2u2=no*nou2
      DO 2 I=1,NO
      IASV=NLAST+I
      READ (NALL4,REC=IASV)TI
      DO 1 A=1,NU
      DO 1 B=1,NU
      DO 1 J=1,NO
      VOE(I,A,B,J)=TI(A,B,J)
 1    CONTINUE
 2    CONTINUE
      RETURN
      END
      SUBROUTINE RDT32S(K,L,NO,NU,TI,T3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      DIMENSION TI(NU,NU,NU),T3(NO,NU,NU,NU)
      DATA TWO/2.0D+0/
      NU3=NU*NU*NU
      DO 100 I=1,NO
      IF (I.EQ.K.AND.I.EQ.L)THEN
      CALL ZEROMA(TI,1,NU3)
      GO TO 90
      ENDIF
      CALL RDVT3ONW(I,K,L,NU,TI)
 90   CONTINUE
      DO 110 A=1,NU
      DO 110 B=1,NU
      DO 110 C=1,NU
      T3(I,A,B,C)=TWO*TI(A,B,C)-Ti(c,b,a)-ti(b,a,c)
 110  CONTINUE
 100  CONTINUE
      RETURN
      END
