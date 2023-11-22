
      SUBROUTINE GSORTHOG(IRREPX,LENGTH,NDIM,VINPUT,TMP1,TMP2,
     &                    ILIST1,ILIST2,SPINAD,IOLDEST,MAXORD,
     &                    OVRLAP,RESID)
C
C FORMS GRAM-SCHMIDT ORTHOGONALIZATION OF INPUT VECTOR VINPUT 
C TO THOSE RESIDING ON LIST ILIST2.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SPINAD
      DIMENSION VINPUT(LENGTH),TMP1(LENGTH),TMP2(LENGTH)
      DIMENSION OVRLAP(NDIM)
      IGET(I)=1+MOD(IOLDEST+MAXORD-I,MAXORD+1)
C
      CALL DCOPY(LENGTH,VINPUT,1,TMP1,1)
      IF (SPINAD) CALL SPNTSING(IRREPX,TMP1,TMP2,1)
      resid=ddot(length,tmp1,1,vinput,1)
      if (resid.gt.0.d0) then
         resid=sqrt(resid)
      else
         if (resid.eq.0.d0) then
            print *, '@GSORTHOG: Error dividing by zero.'
         else
            print *, '@GSORTHOG: Error taking square root of negative.'
         end if
         call aces_exit(1)
      end if
      call dscal(length,1.0d0/resid,vinput,1)
      CALL DCOPY(LENGTH,VINPUT,1,TMP1,1)
      IF (SPINAD) CALL SPNTSING(IRREPX,TMP1,TMP2,1)

      DO 10 I=1,NDIM
       CALL GETLST(TMP2,IGET(I),1,1,ILIST1,ILIST2)
       OVRLAP(I)=-SDOT(LENGTH,TMP1,1,TMP2,1)
       CALL DAXPY(LENGTH,OVRLAP(I),TMP2,1,VINPUT,1)
10    CONTINUE
C
      CALL DCOPY(LENGTH,VINPUT,1,TMP1,1)
      IF (SPINAD) CALL SPNTSING(IRREPX,TMP1,TMP2,1)
      Z=DDOT(LENGTH,VINPUT,1,TMP1,1)
      RESID2=SQRT(Z)
      CALL DSCAL(LENGTH,1.0D0/RESID2,VINPUT,1)
C
      RETURN
      END
