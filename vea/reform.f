      SUBROUTINE REFORM(R, P, SCR,NDIM,NDIMVEC,IRREPX,MAXCOR,
     $   ISIDE,IUHF,ICALC)
C
C THIS ROUTINE CALCULATES THE COMPLETE R MATRIX FROM THE
C EXPANSION VECTORS AND MATRIX VECTOR PRODUCTS ON DISK.
C ALSO THE OVERLAPMATRIX P OF THE PROJECTED VECTORS IS CALCULATED
C
CEND
      IMPLICIT INTEGER (A-Z)
      PARAMETER(MAXORD=100)
      DOUBLE PRECISION SCR,R,SDOT,X,P
      DIMENSION SCR(MAXCOR), R(10000),P(10000)
      COMMON/LISTDAV/LISTC, LISTHC, LISTH0
      COMMON/EXTINF/ITER,IOLDEST
C
      IGET(I)=1+MOD(IOLDEST+MAXORD-I,MAXORD+1)
      INDXF(I,J,N)=I+(J-1)*N
C
C  FIRST DETERMINE OVERLAP FROM VECTORS PROJECTED ON EXCITATION PATTERN (ON COLUMN 5 LISTH0)
C
      NSIZEC = NDIMVEC
      I000=1
      I010=I000+NSIZEC
      I020=I010+NSIZEC
      I030=I020+NSIZEC
      CALL GETLST(SCR(I020),5,1,1,1,LISTH0)
      DO 5 I=1, NDIM
      CALL GETLST(SCR(I010),IGET(I),1,1,ISIDE,LISTC)
      CALL VECPRD(SCR(I020),SCR(I010),SCR(I000),NSIZEC)
      IF(IUHF.EQ.0)THEN
       CALL SPNTSING(NSIZEC,SCR(I000),SCR(I030),MAXCOR-I030+1,IRREPX,
     $      ICALC)
      ENDIF
      DO 7 J=1,NDIM
         CALL GETLST(SCR(I010),IGET(J),1,1,ISIDE,LISTC)
         P(INDXF(I,J,MAXORD)) = SDOT(NSIZEC,SCR(I000),1,SCR(I010),1)
         P(INDXF(J,I,MAXORD)) = P(INDXF(I,J,MAXORD))
 7    CONTINUE
 5    CONTINUE
C
C REFORM R MATRIX
C
      IOFF1=1
      IOFF2=IOFF1+NDIMVEC
      do 10 i=1,ndim
       call getlst(scr(ioff1),iget(i),1,1,iside,listc)
       if(iuhf.eq.0)then
        call spntsing(ndimvec,scr(ioff1),scr(ioff2),maxcor-ioff2+1,
     $       irrepx, icalc)
       endif
       do 11 j=1,ndim
        call getlst(scr(ioff2),iget(j),1,1,iside,listhc)
        x=sdot(ndimvec,scr(ioff1),1,scr(ioff2),1)
        ipos=indxf(i,j,maxord)
        r(ipos)=x
11     continue
10    continue
c
      return
      end
