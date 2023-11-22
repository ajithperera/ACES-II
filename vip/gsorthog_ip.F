      SUBROUTINE GSORTHOG_IP(IRREPX,LENGTH,NDIM,VINPUT,TMP1,TMP2,
     &                       ILIST1,ILIST2,SPINAD,IOLDEST,MAXORD,
     &                       RESID2, ICALC)
C
C FORMS GRAM-SCHMIDT ORTHOGONALIZATION OF INPUT VECTOR VINPUT 
C TO THOSE RESIDING ON LIST ILIST2.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SPINAD
      DIMENSION VINPUT(LENGTH),TMP1(LENGTH),TMP2(2*LENGTH)
      IGET(I)=1+MOD(IOLDEST+MAXORD-I,MAXORD+1)
C
      CALL SCOPY(LENGTH,VINPUT,1,TMP1,1)
      IF(SPINAD)THEN
       CALL SPNTSING(LENGTH,TMP1,TMP2,2*LENGTH, IRREPX, ICALC)
      ENDIF 
      resid=sqrt(sdot(length,tmp1,1,vinput,1))
      call sscal(length,1.0d0/resid,vinput,1)
      CALL SCOPY(LENGTH,VINPUT,1,TMP1,1)
      IF(SPINAD)THEN
       CALL SPNTSING(LENGTH,TMP1,TMP2,2*LENGTH, IRREPX, ICALC)
      ENDIF 
c      write(6,"(a)") "Gsorth-1st"
c      write(6,"(6(1x,F10.5))") (tmp1(i),i=1,length)
c      call checksum("s-old   :", tmp1,length)
c      write(6,*)

      DO 10 I=1,NDIM
       CALL GETLST(TMP2,IGET(I),1,1,ILIST1,ILIST2)
c        write(6,"(a)") "loop"
c        write(6,"(6(1x,F10.5))") (tmp2(j),j=1,length)
c        call checksum("s-old   :", tmp2,length)
       OVRLAP =-SDOT(LENGTH,TMP1,1,TMP2,1)
c        write(6,"(F10.5)") OVRLAP
       CALL SAXPY(LENGTH,OVRLAP,TMP2,1,VINPUT,1)
10    CONTINUE
c      write(6,"(a)") "after"
c      write(6,"(6(1x,F10.5))") (vinput(i),i=1,length)
c      call checksum("vinput  :", vinput,length)
C
      CALL SCOPY(LENGTH,VINPUT,1,TMP1,1)
      IF(SPINAD)THEN
       CALL SPNTSING(LENGTH,TMP1,TMP2,2*LENGTH, IRREPX, ICALC)
      ENDIF 
      Z=SDOT(LENGTH,VINPUT,1,TMP1,1)
      RESID2=SQRT(Z)
C        IF (RESID2 .LE. 0.1) THEN
C           WRITE(6,*) 'NEAR DEGENERACY IN GSORTHOG', RESID2
C        ENDIF
      CALL SSCAL(LENGTH,1.0D0/RESID2,VINPUT,1)
C
      RETURN
      END
