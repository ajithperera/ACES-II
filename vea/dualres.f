      SUBROUTINE DUALRES(SCR, MAXCOR, IRREPX, ISIDE, IUHF, ISPIN)
C
C  NEW EXPANSION VECTOR IS DETERMINED BASED ON THE DUAL SPACE
C  ERROR VECTOR.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION SCR(MAXCOR)
C
      COMMON/LISTDAV/LISTC, LISTHC, LISTH0
      COMMON/EXTINF/NDIMR,IOLDEST
      COMMON/EXTINF2/ROOT
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)
C
         I000 = 1
         I010 = I000 + NSIZEC
         I020 = I010 + NSIZEC
         CALL GETLST(SCR(I000), 1+ISIDE, 1, 1, 1, LISTH0)
         CALL PUTS(SCR(I000),NSIZEC,ISPIN,IUHF,LS1IN,LS2IN)
         CALL EADIR(SCR(I010), (MAXCOR-I010+1)*IINTFP,IUHF,
     $      3-ISIDE, ISPIN)
         CALL LOADS(SCR(I010),NSIZEC,ISPIN,IUHF,LS1OUT,LS2OUT,.FALSE.,
     $      SCR(I020), MAXCOR-I020+1)
         CALL SAXPY(NSIZEC,-ROOT,SCR(I000),1,SCR(I010),1)
         CALL SCOPY(NSIZEC,SCR(I010),1,SCR(I000),1)
C
         RETURN
         END