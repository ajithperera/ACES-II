      FUNCTION SS(NNA,LLA,MM,NNB,LLB,ALPHA,BETA,D,Y,Z)
      IMPLICIT NONE
C-----------------------------------------------------------------------
      DOUBLE PRECISION SS,ALPHA,BETA,D,Y,Z
      INTEGER NNA,LLA,MM,NNB,LLB
C-----------------------------------------------------------------------
      DOUBLE PRECISION A,B,P,PT,X
      INTEGER NA,LA,M,NB,LB,K,I,J,LAMBDA
C-----------------------------------------------------------------------
      DIMENSION D(0:3,0:3,0:3),Y(0:8,0:8,203),Z(0:16,0:8,0:8)
      DIMENSION A(0:16),B(0:16)
C-----------------------------------------------------------------------
C
      NA=NNA
      LA=LLA
      M=MM
      NB=NNB
      LB=LLB
      P =(ALPHA + BETA)/2.0D+00
      PT=(ALPHA - BETA)/2.0D+00
      X =0.0D+00
      M =IABS(M)
C
      IF( LB.LT.LA .OR. (LB.EQ.LA .AND. NB.LT.NA))THEN
       K =NA
       NA=NB
       NB=K
       K =LA
       LA=LB
       LB=K
       PT=-PT
       write(6,*) ' @SS-I, switch was done ! '
      ENDIF
C
      CALL AINTGS(A,P ,NA+NB)
      CALL BINTGS(B,PT,NA+NB)
      write(6,*) ' a values ',a
      write(6,*) ' b values ',b
C
      LAMBDA=(5-M)*(24-10*M+M**2)*(83-30*M+3*M**2)    /120 +
     &       (30-9*LA+LA**2-2*NA)*(28-9*LA+LA**2-2*NA)/  8 +
     &       (30-9*LB+LB**2-2*NB)                     /  2
C
      DO 20 J=0,NA+NB
      DO 10 I=0,NA+NB
C
      X=X+Y(I,J,LAMBDA)*A(I)*B(J)
      if(m.eq.1) write(6,*) i,j,a(i),b(j),y(i,j,lambda)
C
   10 CONTINUE
   20 CONTINUE
C
      SS = D(M,LLA,LLB) * X
C
C     Sign change in certain cases involving pz.
C
C     pz on both.
c      IF(LLA.EQ.LLB .AND. LLA.EQ.1 .AND. MM.EQ.0) SS = -SS
C     s  on A, pz on B.
c      IF(LLA.EQ.0   .AND. LLB.EQ.1 .AND. MM.EQ.0) SS = -SS
C
      IF(LLB .EQ. 1 .AND. MM .EQ. 0) SS = -SS
C
      RETURN
      END
