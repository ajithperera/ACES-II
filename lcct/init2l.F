C     Initializes lists needed for Lambda-based CCSD(T)
C
      SUBROUTINE INIT2L(IUHF)
      IMPLICIT INTEGER(A-Z)
      LOGICAL TRIPNI,TRIPNI1,TRIPIT,T3STOR
      COMMON/TRIPLES/TRIPNI,TRIPNI1,TRIPIT,T3STOR
      COMMON /SYM/ POP1(8),POP2(8),VRT1(8),VRT2(8),NT(2),
     &   NF1(2),NF2(2)
      COMMON /T2ILIS/ LIST2I1,LIST2I2,LIST2I3
#include "syminf.com" /* for nirrep */
C
C CREATE LISTS FOR DOUBLES IN GRADIENT CALCULATIONS
C
      IMODE=0
c      CALL INIPCK(1,13,14,146,IMODE,0,1)
c      CALL INIPCK(1, 1, 3,144,IMODE,0,1)
c      IF (IUHF.EQ.1) CALL INIPCK(1,2,4,145,IMODE,0,1)
C
C CREATE LISTS FOR SINGLES IN GRADIENT CALCULATIONS
C
c       DO 20 ISPIN=1,IUHF+1
c        CALL UPDMOI(1,NT(ISPIN),ISPIN,190,0,0)
c        CALL UPDMOI(1,NT(ISPIN),ISPIN+2,190,0,0)
c        call aces_list_memset(ispin,  190,0)
c        call aces_list_memset(ispin+2,190,0)
c20     CONTINUE
C
c      IF (IFLAGS(87).EQ.3.AND.(TRIPIT.OR.TRIPNI)) RETURN
C
C CREATE LISTS FOR TRIPLES IN GRADIENT CALCULATIONS
C
C  ACTUALLY LISTS G(AB;CI) AND G(IJ;KA)
C

      LISTOFF1=55
      LISTOFF2=39
      LISTOFF3=59
      IF (IUHF.EQ.1) THEN
         CALL INIPCK(1, 3,16,107+LISTOFF1,IMODE,0,1)
         CALL INIPCK(1, 4,17,108+LISTOFF1,IMODE,0,1)
         CALL INIPCK(1,14,11,109+LISTOFF1,IMODE,0,1)
      END IF
      CALL INIPCK(1,14,18,110+LISTOFF1,IMODE,0,1)
      IF (IUHF.EQ.1) THEN
         CALL INIPCK(1, 1, 9,127+LISTOFF2,IMODE,0,1)
         CALL INIPCK(1, 2,10,128+LISTOFF2,IMODE,0,1)
         CALL INIPCK(1,13,18,129+LISTOFF2,IMODE,0,1)
      END IF
      CALL INIPCK(1,13,11,130+LISTOFF2,IMODE,0,1)
c     YAU - initialize all gamma lists to zero
      if (iuhf.eq.0) then
         do iGrp = 1, nirrep
            call aces_list_memset(iGrp,110+LISTOFF1,0)
         end do
         do iGrp = 1, nirrep
            call aces_list_memset(iGrp,130+LISTOFF2,0)
         end do
      else
         do iFam = 107+LISTOFF1, 110+LISTOFF1
            do iGrp = 1, nirrep
               call aces_list_memset(iGrp,iFam,0)
            end do
         end do
         do iFam = 127+LISTOFF2, 130+LISTOFF2
            do iGrp = 1, nirrep
               call aces_list_memset(iGrp,iFam,0)
            end do
         end do
      end if
C
C CREATE LISTS FOR TRIPLE CONTRIBUTION TO THE DENSITY MATRICES
C
      CALL UPDMOI(1,NF1(1),1,160,IMODE,0)
      IF (IUHF.EQ.1) CALL UPDMOI(1,NF1(2),2,160,IMODE,0)
      CALL UPDMOI(1,NF2(1),3,160,IMODE,0)
      IF (IUHF.EQ.1) CALL UPDMOI(1,NF2(2),4,160,IMODE,0)
      CALL UPDMOI(1,NT(1), 5,160,IMODE,0)
      IF (IUHF.EQ.1) CALL UPDMOI(1,NT(2), 6,160,IMODE,0)
      do iGrp = 1, 6, 2-iuhf
         call aces_list_memset(iGrp,160,0)
      end do
C
C CREATE LISTS FOR DELTA T(IJ,AB) CONTRIBUTION TO LAMBDA
C (QCISD AND CCSD ONLY)
C
      IF (IUHF.EQ.1) THEN
         CALL INIPCK(1,1,3,LIST2I1,IMODE,0,1)
         CALL INIPCK(1,2,4,LIST2I2,IMODE,0,1)
         do iFam = LIST2I1, LIST2I2
            do iGrp = 1, nirrep
               call aces_list_memset(iGrp,iFam,0)
            end do
         end do
         CALL INIPCK(1,1,3,201,IMODE,0,1)
         CALL INIPCK(1,2,4,202,IMODE,0,1)
         do iFam = 201, 202
            do iGrp = 1, nirrep
               call aces_list_memset(iGrp,iFam,0)
            end do
         end do
      END IF
      CALL INIPCK(1,13,14,LIST2I3,IMODE,0,1)
      do iGrp = 1, nirrep
         call aces_list_memset(iGrp,LIST2I3,0)
      end do
      CALL INIPCK(1,13,14,203,IMODE,0,1)
      do iGrp = 1, nirrep
         call aces_list_memset(iGrp,203,0)
      end do
      CALL UPDMOI(1,NT(1),1,LIST2I3+1,IMODE,0)
      IF (IUHF.EQ.1) THEN
         CALL UPDMOI(1,NT(2),2,LIST2I3+1,IMODE,0)
      END IF
      call aces_list_memset(1,LIST2I3+1,0)
      IF (IUHF.EQ.1) THEN
         call aces_list_memset(2,LIST2I3+1,0)
      END IF
      IF (IUHF.EQ.1) THEN
         CALL INIPCK(1,1,3,114+LISTOFF3,IMODE,0,1)
         CALL INIPCK(1,2,4,115+LISTOFF3,IMODE,0,1)
         do iFam = 114+LISTOFF3, 115+LISTOFF3
            do iGrp = 1, nirrep
               call aces_list_memset(iGrp,iFam,0)
            end do
         end do
      END IF
      CALL INIPCK(1,13,14,116+LISTOFF3,IMODE,0,1)
      do iGrp = 1, nirrep
         call aces_list_memset(iGrp,116+LISTOFF2,0)
      end do

      RETURN
      END

