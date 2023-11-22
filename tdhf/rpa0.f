      SUBROUTINE RPA0(S,WVAL,REDV,APVEC,APVEC1,NVOSYM,NREDUC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION APVEC(NVOSYM),APVEC1(NVOSYM),WVAL(NVOSYM),
     X S(NVOSYM,NREDUC),REDV(NREDUC)
      DATA ZERO,ONE,TWO,FOUR/0.D0,1.D0,2.D0,4.D0/ 
      WRITE(6,*) ' Now we are in RPA0'
      SQRT2= SQRT(TWO)
C     write(6,*) ' S '
C     CALL OUTMXD(S,NVOSYM,NVOSYM,NVOSYM)    
C .... U+U =  S*(+w)**(-1)*(Z**T*h+Y**T*h)-S*(-w)**(-1)*(Y**T*h+Z**T*h)
C ....    =  2*S*w**(-1)*(Z+Y)*h=  2*S*w**(-1)*S**T*h
C ...   U =  S*w**(-1)*S**T*h
      DO 40 J=1,NREDUC
      REDV(J)=ZERO
      DO 41 I=1,NVOSYM
      REDV(J)= REDV(J) + S(I,J)*APVEC1(I) 
   41 CONTINUE
   40 CONTINUE
      DO 42 I=1,NVOSYM
      APVEC(I)=ZERO
      DO 43 J=1,NREDUC
      APVEC(I) = APVEC(I) + S(I,J)*REDV(J)/WVAL(J)
   43 CONTINUE
   42 CONTINUE
      RETURN
      END
