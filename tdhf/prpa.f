      SUBROUTINE PRPA(S,T,ZVEC,YVEC,REDVEC,WVAL,APVEC,AMVEC,APVEC1
     X ,AMVEC1,OMEGA,NVOSYM,NREDUC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION APVEC(NVOSYM),AMVEC(NVOSYM),APVEC1(NVOSYM),
     X AMVEC1(NVOSYM),S(NVOSYM,NREDUC),T(NVOSYM,NREDUC),
     X ZVEC(NVOSYM,NREDUC),YVEC(NVOSYM,NREDUC),REDVEC(NREDUC,2)
     X ,WVAL(NVOSYM)
      DATA ZERO,ONE,TWO,FOUR/0.D0,1.D0,2.D0,4.D0/ 
      WRITE(6,*) ' Now we are in PRPA ', NREDUC 
      SQRT2= SQRT(TWO)
C .... U(+)+U(-) = -S*(omega+w)**(-1)*(Z**T*h(+)+Y**T*h(-))
C .....            +S*(omega-w)**(-1)*(Y**T*h(+)+Z**T*f(-)  .........
C .... U(+)-U(-) = -T*(omega+w)**(-1)*(Z**T*h(+)+Y**T*h(-))
C .....            -T*(omega-w)**(-1)*(Y**T*h(+)+Z**T*f(-)  .........
C     write(6,*) ' S '
C     CALL OUTMXD(S,NVOSYM,NVOSYM,NVOSYM)    
C     write(6,*) ' T '
C     CALL OUTMXD(T,NVOSYM,NVOSYM,NVOSYM)    
C     write(6,*) ' Z '
C     CALL OUTMXD(ZVEC,NVOSYM,NVOSYM,NVOSYM)    
C     write(6,*) ' Y '
C     CALL OUTMXD(YVEC,NVOSYM,NVOSYM,NVOSYM)    
      DO 40 J=1,NREDUC
      REDV1=ZERO
      REDV2=ZERO
      DO 41 I=1,NVOSYM
      REDV1= REDV1 + ZVEC(I,J)*APVEC1(I) +YVEC(I,J)*AMVEC1(I)
      REDV2= REDV2 + ZVEC(I,J)*AMVEC1(I) +YVEC(I,J)*APVEC1(I)
   41 CONTINUE
      REDVEC(J,1) =  -REDV1/(OMEGA+WVAL(J)) + REDV2/(OMEGA-WVAL(J))
      REDVEC(J,2) =  -REDV1/(OMEGA+WVAL(J)) - REDV2/(OMEGA-WVAL(J))
   40 CONTINUE
      DO 42 I=1,NVOSYM
      APVEC(I)=ZERO
      AMVEC(I)=ZERO
      DO 43 J=1,NREDUC
      APVEC(I) = APVEC(I) + S(I,J)*REDVEC(J,1)
      AMVEC(I) = AMVEC(I) + T(I,J)*REDVEC(J,2)
   43 CONTINUE
      APVEC(I)=APVEC(I)/SQRT2
      AMVEC(I)=AMVEC(I)/SQRT2
   42 CONTINUE
      RETURN
      END
