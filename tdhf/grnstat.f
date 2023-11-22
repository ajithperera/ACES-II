      SUBROUTINE GRNSTAT(ABP,ABM,ABPM,UVAL,UVEC,WVAL,WVEC,PORT,
     X APVEC,AMVEC,APVEC1,NVOSYM,NREDUC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ABPM(NVOSYM,NVOSYM),PORT(NVOSYM,NVOSYM),WVAL(NVOSYM)
     X ,WVEC(NVOSYM,NVOSYM),ABP(NVOSYM,NVOSYM)
     X ,ABM(NVOSYM,NVOSYM),UVAL(NVOSYM),UVEC(NVOSYM,NVOSYM)
     X ,APVEC(NVOSYM),AMVEC(NVOSYM),APVEC1(NVOSYM)
C  ..... SCR(NVO+1,5)
C     DIMENSION TEST1(35,35),TEST2(35,35)
      COMMON/THRES/ TOLPER,DEGEN,EPSI
      COMMON/IPRNT/ IOPEV,IOPDA,IOPU,IOPFE,IOPPR,IWRPA
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      DATA ZERO,ONE,TWO,FOUR/0.D0,1.D0,2.D0,4.D0/ 
      WRITE(6,*) ' We are in GRNSTAT '
      SQRT2= SQRT(TWO)
      write(6,*) ' NVOSYM,NREDUC = ',NVOSYM,NREDUC
C     write(6,*) ' ABP in GRNSTAT'
C     CALL OUTMXD(ABP,NVOSYM,NVOSYM,NVOSYM)    
C     WRITE(6,*) ' (A+B)*(A-B) matrix in GRNSTAT'
C     CALL OUTMXD(ABPM,NVOSYM,NVOSYM,NVOSYM)    
C     WRITE(6,*) ' UVEC in GRNSTAT'
C     CALL OUTMXD(UVEC,NVOSYM,NVOSYM,NVOSYM)    
C     WRITE(6,*) ' WVEC in GRNSTAT'
C     CALL OUTMXD(WVEC,NVOSYM,NVOSYM,NVOSYM)    
C ........................................................................
      IF(IWRPA.NE.0) THEN
      DO 1 I=1,NREDUC
      IF(WVAL(I).LT.ZERO) THEN
      WRITE(6,*) ' instable at ',I,'th vector'
      ELSE
      EN = SQRT(WVAL(I))
      WRITE(6,*) I,'th excitation in RPA = ',EN,WVAL(I)
      END IF
    1 CONTINUE
      END IF
      DO 17 I=1,NVOSYM
      DO 17 J=1,NREDUC
      ABM(I,J)=WVEC(I,J)/WVAL(J)
   17 CONTINUE
C ..... PORT = WVEC*(WVAL)**(-1)*(WVEC)t  ......
      CALL XGEMM('N','T',NVOSYM,NVOSYM,NREDUC,ONE,ABM,NVOSYM,
     X WVEC,NVOSYM,ZERO,PORT,NVOSYM)
C     WRITE(6,*) ' inversed matrix in GRNSTAT '
C     CALL OUTMXD(PORT,NVOSYM,NVOSYM,NVOSYM)
C   ...... inverse owtta .....
C   .....  which is a inversed matrix: PORT
C   ..... Calculate perturbed vectors using input vectors 
C   .....  the Porlarization Green's Function
C   .....  which is a inversed matrix: PORT
C   ..... AMVEC =  (A+B)*(A-B)*APVEC1, APVEC1= - h: constant rhs  .....
      CALL XGEMM('N','N',NVOSYM,1,NVOSYM,ONE,ABPM,NVOSYM,
     X APVEC1,NVOSYM,ZERO,AMVEC,NVOSYM)
C ....  APVEC = (A+B)**(-1/2)*AMVEC
C             = -(A+B)**(1/2)*(A-B)*h.....
      CALL XGEMM('T','N',NVOSYM,1,NVOSYM,ONE,UVEC,NVOSYM,
     X AMVEC,NVOSYM,ZERO,APVEC,NVOSYM)
C ..... AMVEC = PORT*APVEC = (A+B)**(1/2)*U  .............
      CALL MXM(PORT,NVOSYM,APVEC,NVOSYM,AMVEC,1) 
C ...... APVEC = (A+B)**(-1/2)*AMVEC = U  .......................
      CALL XGEMM('N','N',NVOSYM,1,NVOSYM,ONE,UVEC,NVOSYM,
     X AMVEC,NVOSYM,ZERO,APVEC,NVOSYM)
      RETURN
      END
