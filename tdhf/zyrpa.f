      SUBROUTINE ZYRPA(ABP,ABM,UVAL,UVEC,WVAL,WVEC,PORT,SCR,
     X S,T,ZVEC,YVEC,NVOSYM,NREDUC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PORT(NVOSYM,NVOSYM),WVAL(NVOSYM)
     X ,WVEC(NVOSYM,NVOSYM),SCR(1),ABP(NVOSYM,NVOSYM)
     X ,ABM(NVOSYM,NVOSYM),UVAL(NVOSYM),UVEC(NVOSYM,NVOSYM)
     X  ,S(NVOSYM,NREDUC),T(NVOSYM,NREDUC),
     X ZVEC(NVOSYM,NREDUC),YVEC(NVOSYM,NREDUC)
      COMMON/THRES/ TOLPER,DEGEN,EPSI
      COMMON/IPRNT/ IOPEV,IOPDA,IOPU,IOPFE,IOPPR,IWRPA
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      DATA ZERO,ONE,TWO,FOUR/0.D0,1.D0,2.D0,4.D0/ 
CSSS      WRITE(6,*) ' We are in ZYRPA '
CSSS      write(6,*) ' NVOSYM,NREDUC = ',NVOSYM,NREDUC
C ........................................................................
CSSS      write(6,*) "A + B in TDSOLVE" 
CSSS      CALL output(ABP, 1, NVOSYM, 1, NVOSYM, NVOSYM, NVOSYM, 1)
CSSS      write(6,*) "A - B in TDSOLVE" 
CSSS      CALL output(ABM, 1, NVOSYM, 1, NVOSYM, NVOSYM, NVOSYM, 1)

C     CALL OUTMXD(ABP,NVOSYM,NVOSYM,NVOSYM)    
C     write(6,*) ' ABM '
C     CALL OUTMXD(ABM,NVOSYM,NVOSYM,NVOSYM)    
      CALL ICOPY(NVOSYM*NVOSYM*IINTFP,ABP,1,PORT,1)
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
C ....... Now UVEC = U*X**(1/2) ........
CSSS      Write(6,*) "The sqrt of A+B"
CSSS      CALL output(UVEC, 1, NVOSYM, 1, NVOSYM, NVOSYM, NVOSYM, 1)
      
      CALL XGEMM('T','N',NVOSYM,NVOSYM,NVOSYM,ONE,UVEC,NVOSYM,ABM
     X ,NVOSYM,ZERO,PORT,NVOSYM)
      CALL XGEMM('N','N',NVOSYM,NVOSYM,NVOSYM,ONE,PORT,NVOSYM,UVEC
     X ,NVOSYM,ZERO,ABM,NVOSYM)
C     ABM: W**2 = (U*X**(1/2))t*(A-B)*(U*X**(1/2)) 
C     WRITE(6,*) ' W**2 = (U*X**(1/2))t*(A-B)*(U*X**(1/2)) matrix '
C     CALL OUTMXD(ABM,NVOSYM,NVOSYM,NVOSYM)    
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
      IF(IWRPA.GT.0) THEN
      DO 2 I=1,NREDUC
      WEV= 27.2113957D0*WVAL(I)
      WCMI=WEV*8065.54093D0
      WCM=1.D0/WCMI
      WRITE(6,*) I,'th excitation in RPA = ',WVAL(I),' a.u. ',
     X WEV,' eV ',WCMI,' /cm ',WCM,' cm'
    2 CONTINUE
      END IF
C  ..... end of IWRPA.NE.0 ...
CSSS      Write(6,*) "The eigenvectors"
CSSS      CALL output(WVEC, 1, NVOSYM, 1, NVOSYM, NVOSYM, NVOSYM, 1)
      DO 17 J=1,NREDUC
      WSQRT= SQRT(WVAL(J))
      DO 17 I=1,NVOSYM
      ZVEC(I,J)=WVEC(I,J)*WSQRT
      YVEC(I,J)=WVEC(I,J)/WSQRT
   17 CONTINUE
CSSS      Write(6,*) "The ZVEC"
CSSS      CALL output(ZVEC, 1, NVOSYM, 1, NVOSYM, NVOSYM, NVOSYM, 1)
CSSS      Write(6,*) "The YVEC"
CSSS      CALL output(YVEC, 1, NVOSYM, 1, NVOSYM, NVOSYM, NVOSYM, 1)

C ....    T  = U*X**(1/2)*V*w**(-1/2) .........
      CALL XGEMM('N','N',NVOSYM,NREDUC,NVOSYM,ONE,UVEC,NVOSYM,YVEC
     X ,NVOSYM,ZERO,T,NVOSYM)
      DO 16 I=1,NVOSYM
      DO 16 J=1,NVOSYM
   16 UVEC(I,J)= UVEC(I,J)/UVAL(J)
C  ......... Now UVEC =  U*X**(-1/2) ........
C     write(6,*) ' U*X**(-1/2) '
C     CALL OUTMXD(UVEC,NVOSYM,NVOSYM,NVOSYM)    
C ....    S  = U*X**(-1/2)*V*w**(1/2) .........
      CALL XGEMM('N','N',NVOSYM,NREDUC,NVOSYM,ONE,UVEC,NVOSYM,ZVEC
     X ,NVOSYM,ZERO,S,NVOSYM)
CSSS      Write(6,*) "The TVEC"
CSSS      CALL output(T, 1, NVOSYM, 1, NVOSYM, NVOSYM, NVOSYM, 1)
CSSS      Write(6,*) "The SVEC"
CSSS      CALL output(S, 1, NVOSYM, 1, NVOSYM, NVOSYM, NVOSYM, 1)
      DO 18 J=1,NREDUC
      DO 18 I=1,NVOSYM
      ZVEC(I,J) = (S(I,J)+T(I,J))/TWO
      YVEC(I,J) = (S(I,J)-T(I,J))/TWO
   18 CONTINUE
CSSS      Write(6,*) "The ZVEC"
CSSS      CALL output(ZVEC, 1, NVOSYM, 1, NVOSYM, NVOSYM, NVOSYM, 1)
CSSS      Write(6,*) "The YVEC"
CSSS      CALL output(YVEC, 1, NVOSYM, 1, NVOSYM, NVOSYM, NVOSYM, 1)

CSSS      IF(IWRPA.GT.1) THEN
CSSS      write(6,*) ' Z vector'
CSSS      CALL OUTMXD(ZVEC,NVOSYM,NVOSYM,NVOSYM)    
CSSS      write(6,*) ' Y vector'
CSSS      CALL OUTMXD(YVEC,NVOSYM,NVOSYM,NVOSYM)    
CSSS      END IF
      RETURN
      END
