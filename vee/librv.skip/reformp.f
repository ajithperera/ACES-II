      SUBROUTINE REFORMP(P, SCR,NDIM,NDIMVEC,IRREPX,MAXCOR,
     $   ISIDE,IUHF)
C
C THIS ROUTINE CALCULATES THE OVERLAPMATRIX P OF THE PROJECTED VECTORS
C
CEND
      IMPLICIT INTEGER (A-Z)
      PARAMETER(MAXORD=100)
      DOUBLE PRECISION SCR,SDOT,P
      DIMENSION SCR(MAXCOR), P(10000)
      COMMON/EXTINF/ITER,IOLDEST
      COMMON/LISTPROJ/LISTH0, ICOLPR1, ICOLPR2
C
      IGET(I)=1+MOD(IOLDEST+MAXORD-I,MAXORD+1)
      INDXF(I,J,N)=I+(J-1)*N
C
      LISTC = 470
C
C DETERMINE OVERLAP FROM VECTORS PROJECTED ON EXCITATION PATTERN (ON ICOLPR1 LISTH0)
C
      NSIZEC = NDIMVEC
      I000=1
      I010=I000+NSIZEC
      I020=I010+NSIZEC
      I030=I020+NSIZEC
      CALL GETLST(SCR(I020),ICOLPR1,1,1,1,LISTH0)
      DO 5 I=1, NDIM
      CALL GETLST(SCR(I010),IGET(I),1,1,ISIDE,LISTC)
      CALL VECPRD(SCR(I020),SCR(I010),SCR(I000),NSIZEC)
      IF(IUHF.EQ.0)THEN
        CALL SPNTSING(IRREPX,SCR(I000),SCR(I030),MAXCOR-I030+1)
      ENDIF
      DO 7 J=1,NDIM
         CALL GETLST(SCR(I010),IGET(J),1,1,ISIDE,LISTC)
         P(INDXF(I,J,MAXORD)) = SDOT(NSIZEC,SCR(I000),1,SCR(I010),1)
         P(INDXF(J,I,MAXORD)) = P(INDXF(I,J,MAXORD))
 7    CONTINUE
 5    CONTINUE
C
      RETURN
      END

