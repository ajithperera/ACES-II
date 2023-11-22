
c The Broyden-Fletcher-Goldfarb-Shanno update.

      SUBROUTINE PCCD_BFGS(V,H,SCRATCH,STEP,TBT,N)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION V(N),H(N,N),SCRATCH(7*N*N),TBT(3*N*N),STEP(N)
      DATA DONE,DNULL,EPS /1.0D0,0.0D0,3.0D-08/
C
C N      : Order of the Hessian matrix (and other array)
C V      : The current gradient 
C SCRATCH: Length of 7*N (the first N elements are the previous gradient)
C TBT    : Scratch array length 3*N
C STEP   : Is defined as x(k) = H^-1(k)x(K)

      N2 =1+N
      N3 =1+2*N
      N4 =1+3*N
      NX2=1+N*N
C
C Scratch(N2) = Delta(G);  (G(k-1)-G(k)), computing the gradient difference. 
C
      CALL VADD(SCRATCH(N2),V,SCRATCH(1),N,-1.D0)

      Write(6,"(a)") "-----------printing begins @-pccd-bfgs-----------"
      write(6,"(a,5(F15.10))") "G(k)       :",(v(i),i=1,n)
      write(6,"(a,5(F15.10))") "G(k-1)     :",(scratch(i),i=1,n)
      write(6,"(a,5(F15.10))") "x(k-1)-x(k):",(step(i),i=1,n)
      write(6,"(a,5(F15.10))") "G(k)-G(k-1):",(scratch(n2-1+i),i=1,n)
      write(6,*) "The current Hessian"
      call output(h,1,n,1,n,n,n,1)

C Compute H(k)*Delta(G)

      Call  DGEMM("N","N",N,N,N,DONE,H,N,SCRATCH(N2),N,DNULL,
     +            SCRATCH(N3),N)

      write(6,"(a)") " The H*DG"
      call output(scratch(n3),1,n,1,1,n,1,1)

C Calculate the dot products for the denominators 

      FAC    = DDOT(N,STEP,1,SCRATCH(N2),1)      
      FAE    = DDOT(N,SCRATCH(N3),1,SCRATCH(N2),1)
      SUMDG  = DDOT(N,SCRATCH(N2),1,SCRATCH(N2),1)
      SUMDG  = SUMDG**2
      SUMXI  = DDOT(N,STEP,1,STEP,1)
      SUMXI  = SUMXI**2
      CUTOFF = DSQRT(EPS*SUMDG*SUMXI)
      Write(6,"(a)") " The denominators"
      Write(6,"(a,F15.10)") "FAC    :",FAC
      Write(6,"(a,F15.10)") "FAE    :",FAE
      Write(6,"(a,F15.10)") "SUMDGF :",SUMDG
      Write(6,"(a,F15.10)") "SUMXIF :",SUMXI
      Write(6,"(a,F15.10)") "CUTOFF :",CUTOF

C Skip update if FAC is not sufficiently positive, otherwise implement the
C BFGS update formula 

      IF (FAC .GT. CUTOFF) THEN
          FAC = DONE/FAC
          FAD = DONE/FAE
          DO I = 1, N
             SCRATCH(N2+I-1) = FAC*STEP(I)-FAD*SCRATCH(N3+I-1)
          ENDDO
          DO I = 1, N
             DO J = 1, N
                H(I,J) = H(I,J) + FAC*STEP(I)*STEP(J) -
     +                   FAD*SCRATCH(N3+I-1)*SCRATCH(N3+J-1) +
     +                   FAE*SCRATCH(N2+I-1)*SCRATCH(N2+J-1)
                H(J,I) = H(I,J)
             ENDDO
          ENDDO
      ENDIF 

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "-----The BFGS updated Hessian-----"
      Call output(H,1,N,1,N,N,N,1)
#endif
      Write(6,"(a)") "-----------printing ends @-pccd-bfgs-----------"

      RETURN
      END

