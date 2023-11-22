C         
      SUBROUTINE READDG(DOO,DVV,DVO,DOV,LENOO,LENVV,LENVO,ISPIN)
C
C Reads ground state density from disk.  
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DOO(*),DVV(*),DVO(*),DOV(*)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/LISTDENS/LDENS
      DATA ONE,ONEM,HALF,HALFM/1.0D0,-1.0D0,0.5D0,-0.5D0/
C
C OCC-OCC
C
      CALL GETLST (DOO,1,1,1,ISPIN,ldens)
CSSS      call checksum("DOO  :", Doo,lenoo)
C
C VRT-VRT
C
      CALL GETLST (DVV,1,1,1,2+ISPIN,ldens)
      call checksum("DVV  :", Doo,lenvv)
C
C VRT-OCC
C
      CALL GETLST (DVO,1,1,1,4+ISPIN,ldens)
      CALL GETLST (DOV,1,1,1,ISPIN,90)
      CALL SAXPY  (LENVO,ONE,DOV,1,DVO,1)
CSSS      Print*, "The vo dens"
CSSS      print*,  (0.5D0*dvo(i), i=1, lenvo)
CSSS     call checksum("DvO  :", Doo,lenvo)
C OCC-VRT
C
      CALL GETLST (DOV,1,1,1,ISPIN,190)
CSSS      call checksum("DOV  :", Dov,lenvo)
CSSS      write(6,*)
CSSS   print*,  (dov(i), i=1, lenvo)

      RETURN
      END
