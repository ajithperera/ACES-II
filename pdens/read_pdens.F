C         
      SUBROUTINE READ_PDENS(DOO,DVV,DVO,DOV,LENOO,LENVV,LENVO,
     &                      ISPIN)
C
C Reads ground state density from disk.  
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DOO(*),DVV(*),DVO(*),DOV(*)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      DATA ONE,ONEM,HALF,HALFM/1.0D0,-1.0D0,0.5D0,-0.5D0/
C
C OCC-OCC
C
      CALL GETLST (DOO,1,1,1,ISPIN,160)
C
C VRT-VRT
C
      CALL GETLST (DVV,1,1,1,2+ISPIN,160)
C
C VRT-OCC
C
      CALL GETLST (DVO,1,1,1,4+ISPIN,160)
C
C OCC-VRT
C
      CALL GETLST (DOV,1,1,1,6+ISPIN,160)
#ifdef _DEBUG_LVM2
      write(6,*) 
      print*, "@-READ_PDENS, occ-occ density"
      print*,  (doo(i), i=1, lenoo)
      write(6,*)
      print*, "@-READ_PDENS, vrt-vrt density"
      print*,  (dvv(i), i=1, lenvv)
      write(6,*)
      print*, "@-READ_PDENS, vrt-occ density"
      print*,  (dvo(i), i=1, lenvo)
      write(6,*)
      print*, "@-READ_PDENS, occ-vrt density"
      print*,  (dov(i), i=1, lenvo)
#endif
C
      RETURN
      END
