      SUBROUTINE FINVMX(ABP,ABM,ABPM,UVAL,UVEC,WVAL,WVEC,PORT,SCR,
     X APVEC,AMVEC,APVEC1,AMVEC1,OMEGA,NFIRST,NVOSYM,NREDUC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ABPM(NVOSYM,NVOSYM),PORT(NVOSYM,NVOSYM),WVAL(NVOSYM)
     X ,WVEC(NVOSYM,NVOSYM),SCR(1),ABP(NVOSYM,NVOSYM)
     X ,ABM(NVOSYM,NVOSYM),UVAL(NVOSYM),UVEC(NVOSYM,NVOSYM)
     X ,APVEC(NVOSYM),AMVEC(NVOSYM),APVEC1(NVOSYM),AMVEC1(NVOSYM)
C  ..... SCR(NVO+1,5)
C     DIMENSION TEST1(35,35),TEST2(35,35)
      COMMON/THRES/ TOLPER,DEGEN,EPSI
      COMMON/IPRNT/ IOPEV,IOPDA,IOPU,IOPFE,IOPPR,IWRPA
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      DATA ZERO,ONE,TWO,FOUR/0.D0,1.D0,2.D0,4.D0/ 
      WRITE(6,*) ' We are in FINVMX with NFIRST and ',NFIRST
      SQRT2= SQRT(TWO)
      write(6,*) ' NVOSYM,NREDUC = ',NVOSYM,NREDUC
C ........................................................................
      IF(NFIRST.EQ.0) THEN
      CALL MXM(ABP,NVOSYM,ABM,NVOSYM,ABPM,NVOSYM) 
C     write(6,*) ' ABP in FINVMX'
C     CALL OUTMXD(ABP,NVOSYM,NVOSYM,NVOSYM)    
C     write(6,*) ' ABM '
C     CALL OUTMXD(ABM,NVOSYM,NVOSYM,NVOSYM)    
C     WRITE(6,*) ' (A+B)*(A-B) matrix '
C     CALL OUTMXD(ABPM,NVOSYM,NVOSYM,NVOSYM)    
c YAU : old
c     CALL ICOPY(NVOSYM*NVOSYM*IINTFP,ABP,1,PORT,1)
c YAU : new
      CALL DCOPY(NVOSYM*NVOSYM,ABP,1,PORT,1)
c YAU : end
C     write(6,*) ' copied ? , NVOSYM,IINTFP',NVOSYM,IINTFP
C     CALL OUTMXD(PORT,NVOSYM,NVOSYM,NVOSYM)    
C     CALL RS(NVOSYM,NVOSYM,PORT,UVAL,1,UVEC,SCR,SCR(NVOSYM+1),IERR)
      CALL MHOUSE(NVOSYM,NVOSYM,NVOSYM+1,PORT,UVEC,SCR)
C     CALL MHOUSE(NVOSYM,NVOSYM,PORT,UVEC)
C     WRITE(6,*) ' IERR in (A+B) = ',IERR
      DO 14 I=1,NVOSYM
   14 UVAL(I)=PORT(I,I)
C  ...... Ordering of UVAL and UVEC ...............
      DO 30 I=1,NVOSYM
      DO 30 K=I,NVOSYM
      IF(UVAL(K).LT.UVAL(I)) THEN
      Y = UVAL(I)
      UVAL(I) = UVAL(K)
      UVAL(K) = Y
      DO 31 J=1,NVOSYM
      Y = UVEC(J,I)
      UVEC(J,I) = UVEC(J,K)
   31 UVEC(J,K) = Y
      END IF
   30 CONTINUE
C  ................................................
      DO 15 I=1,NVOSYM
      DO 15 J=1,NVOSYM
   15 UVEC(I,J)= SQRT(UVAL(J))*UVEC(I,J)
C ....... U*X**(1/2) ........
      CALL XGEMM('T','N',NVOSYM,NVOSYM,NVOSYM,ONE,UVEC,NVOSYM,ABM
     X ,NVOSYM,ZERO,PORT,NVOSYM)
      CALL XGEMM('N','N',NVOSYM,NVOSYM,NVOSYM,ONE,PORT,NVOSYM,UVEC
     X ,NVOSYM,ZERO,ABM,NVOSYM)
C     ABM: W**2 = (U*X**(1/2))t*(A-B)*(U*X**(1/2)) 
C     WRITE(6,*) ' W**2 = (U*X**(1/2))t*(A-B)*(U*X**(1/2)) matrix '
C     CALL OUTMXD(ABM,NVOSYM,NVOSYM,NVOSYM)    
      DO 16 I=1,NVOSYM
      DO 16 J=1,NVOSYM
   16 UVEC(I,J)= UVEC(I,J)/UVAL(J)
C  ......... U*X**(-1/2) ........
C ..............................................
C     DO 18 I=1,NVOSYM
C     DO 18 J=1,NVOSYM
C     TEST1(I,J)=ABM(I,J)
C     IF(I.EQ.J) TEST1(I,J) = TEST1(I,J) -OMEGA*OMEGA
C  18  CONTINUE
C     CALL RS(NSIZVO,NVOSYM,ABM,WVAL,1,WVEC,SCR,SCR(NVOSYM+1),IERR)
      CALL MHOUSE(NVOSYM,NVOSYM,NVOSYM+1,ABM,WVEC,SCR)
C     CALL MHOUSE(NVOSYM,NVOSYM,ABM,WVEC)
      DO 13 I=1,NVOSYM
   13 WVAL(I)=ABM(I,I)
C  ...... Ordering of WVAL and WVEC ...............
      DO 20 I=1,NVOSYM
      DO 20 K=I,NVOSYM
      IF(WVAL(K).LT.WVAL(I)) THEN
      Y = WVAL(I)
      WVAL(I) = WVAL(K)
      WVAL(K) = Y
      DO 21 J=1,NVOSYM
      Y = WVEC(J,I)
      WVEC(J,I) = WVEC(J,K)
   21 WVEC(J,K) = Y
      END IF
   20 CONTINUE
C  ................................................
C     WRITE(6,*) ' IERR in W**2 = ',IERR
      DO 1 I=1,NVOSYM
      IF(WVAL(I).LT.ZERO) THEN
      WRITE(6,*) ' instable at ',I,'th vector'
      ELSE
      WVAL(I) = SQRT(WVAL(I))
      END IF
    1 CONTINUE
      IF(IWRPA.NE.0) THEN
      DO 2 I=1,NREDUC
      WEV=27.2113957D0*WVAL(I)
      WRITE(6,*) I,'th excitation in RPA = ',WVAL(I),' a.u. ',
     X WEV,' eV'
    2 CONTINUE
      END IF
C  ..... end of IWRPA.NE.0 ...
      END IF
C   .. end of NFIRST=0 .............................................
C     WRITE(6,*) ' UVEC in FINVMX'
C     CALL OUTMXD(UVEC,NVOSYM,NVOSYM,NVOSYM)    
C     WRITE(6,*) ' WVEC in FINVMX'
C     CALL OUTMXD(WVEC,NVOSYM,NVOSYM,NVOSYM)    
C     IF(OMEGA.EQ.0.D0) THEN
C     OMEGA2 = 0.D0
C     ELSE
      OMEGA2=OMEGA*OMEGA
C     END IF
      DO 17 I=1,NVOSYM
      DO 17 J=1,NREDUC
      WVAL2=WVAL(J)*WVAL(J)
      ABM(I,J)=WVEC(I,J)/(WVAL2-OMEGA2)
   17 CONTINUE
C ..... PORT = WVEC*(WVAL-OMEGA**2)**(-1)*(WVEC)t  ......
      CALL XGEMM('N','T',NVOSYM,NVOSYM,NREDUC,ONE,ABM,NVOSYM,
     X WVEC,NVOSYM,ZERO,PORT,NVOSYM)
C     WRITE(6,*) ' inversed matrix in FINVMX '
C     CALL OUTMXD(PORT,NVOSYM,NVOSYM,NVOSYM)
C   ...... inverse owtta .....
C   .....  which is a inversed matrix: PORT
C   ..... Calculate perturbed vectors using input vectors 
C   .....  APVEC1 and AMVEC1 together with 
C   .....  the Porlarization Green's Function
C   .....  which is a inversed matrix: PORT
C   .....  AMVEC = (A+B)*(h(+)-h(-))*omega ............
      CALL XGEMM('N','N',NVOSYM,1,NVOSYM,OMEGA,ABP,NVOSYM,
     X AMVEC1,NVOSYM,ZERO,AMVEC,NVOSYM)
C   ..... AMVEC = (A+B)*(h(+)-h(-))*omega - (A+B)*(A-B)*(h(+)+h(-)) .....
      CALL XGEMM('N','N',NVOSYM,1,NVOSYM,-ONE,ABPM,NVOSYM,
     X APVEC1,NVOSYM,ONE,AMVEC,NVOSYM)
C ....  APVEC = (A+B)**(-1/2)*AMVEC
C   = (A+B)**(-1/2)*(h(+)-h(-))*omega - (A+B)**(-1/2)*(A-B)*(h(+)+h(-)) .....
      CALL XGEMM('T','N',NVOSYM,1,NVOSYM,ONE,UVEC,NVOSYM,
     X AMVEC,NVOSYM,ZERO,APVEC,NVOSYM)
C ..... AMVEC = PORT*APVEC = (A+B)**(1/2)*(U(+)+U(-)) ..................
      CALL MXM(PORT,NVOSYM,APVEC,NVOSYM,AMVEC,1) 
C ...... APVEC = (A+B)**(-1/2)*AMVEC = U(+)+U(-) .......................
      CALL XGEMM('N','N',NVOSYM,1,NVOSYM,ONE,UVEC,NVOSYM,
     X AMVEC,NVOSYM,ZERO,APVEC,NVOSYM)
C  ...... AMVEC = (A+B)*(U(+)+U(-)) .....
      CALL MXM(ABP,NVOSYM,APVEC,NVOSYM,AMVEC,1) 
C   .........................................
      DO 49 I=1,NVOSYM
      IF(ABS(OMEGA).GT.0.00000001) THEN
C  ..... AMVEC = -((A+B)*(U(+)+U(-))+(h(+)+h(-)))/(sqrt(2)*omega)
C              = (U(+)-U(-))/sqrt(2) 
      AMVEC(I)=-(AMVEC(I)+APVEC1(I))/OMEGA/SQRT2
      ELSE
      AMVEC(I)=ZERO
      END IF
C ..... APVEC = (U(+)+U(-))/sqrt(2)
      APVEC(I)=APVEC(I)/SQRT2
   49 CONTINUE
C...  APVEC:(U(+)+U(-))/sqrt(2) AMVEC:((U(+)-U(-))/sqrt(2) 
C     WRITE(6,*) ' EIG will be called '
C     CALL EIG(A,FVEC,NSIZVO,NSIZVO,0)
C .................................
C     IF(NFIRST.EQ.0) THEN
C  ....... test for inversed matrix .....
C     WRITE(6,*) ' test for inversion '
C     CALL MXM(TEST1,NVO,PORT,NVO,TEST2,NVO)
C     CALL OUTMXD(TEST2,NVO,NVO,NVO)
C     END IF
      RETURN
      END
