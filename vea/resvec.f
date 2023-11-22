      SUBROUTINE RESVEC(SCR, MAXCOR, BUF, MAXORD, ISIDE, IBUFLC)
C
C CALCULATE RESIDUAL : [H C(new) - E C(new)]
C THE RESIDUAL IS PUT IN SCR, THE FIRST NSIZEC ELEMENTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION BUF(MAXORD*(2*MAXORD+3)+5), SCR(MAXCOR)
C
      COMMON/LISTDAV/LISTC, LISTHC, LISTH0
      COMMON/EXTINF/NDIMR,IOLDEST
      COMMON/EXTINF2/ROOT
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
C
      I000=1
      I010=I000+NSIZEC
      I020=I010+NSIZEC
      ITOP=I020+NSIZEC
      CALL FORMS(NSIZEC,NDIMR,SCR(I000),SCR(I020),BUF(IBUFLC),
     &           ISIDE,LISTHC,IOLDEST,MAXORD)
      CALL FORMS(NSIZEC,NDIMR,SCR(I010),SCR(I020),BUF(IBUFLC),
     &           ISIDE,LISTC,IOLDEST,MAXORD)
C
C WRITE CURRENT EIGENVECTOR TO LISTH0
C
C      Write(6, "(6(F10.5))") (scr(i000+i-1),i=1,10)
C      Write(6,*)
C      Write(6, "(6(F10.5))") (scr(i010+i-1),i=1,10)

      CALL PUTLST(SCR(I010), 1+ISIDE, 1, 1, 1, LISTH0)
C
      CALL SAXPY(NSIZEC,-ROOT,SCR(I010),1,SCR(I000),1)
C      Write(6,*) " From RESVEC"
C      Write(6, "(6(F10.5))") (scr(i000+i-1),i=1,nsizec)
C
      RETURN
      END
