C
C   THIS ROUTINE INITIALIZES LISTS ON THE 
C   GAMLAM FILE FOR UCC AND XCC TRIPLES GRADIENT CALCULATIONS
C
C CODED FEB/92 PS
C
      SUBROUTINE INIT3(IUHF)
      IMPLICIT INTEGER(A-Z)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      LOGICAL TRIPNI,TRIPNI1,TRIPIT,T3STOR,NONHF
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON/TRIPLES/TRIPNI,TRIPNI1,TRIPIT,T3STOR
      COMMON/NHFREF/NONHF
      COMMON /SYM/ POP1(8),POP2(8),VRT1(8),VRT2(8),NT(2),
     &             NF1(2),NF2(2)
#include "syminf.com" /* for nirrep */
C
      IMODE=0
c&      IF(.NOT.XCC) THEN
        LISTOFF2=0
        LISTOFF1=0
       IF(IUHF.EQ.1) THEN
        CALL INIPCK(1,3,16,107+LISTOFF1,IMODE,0,1)
        CALL INIPCK(1,4,17,108+LISTOFF1,IMODE,0,1)
        CALL INIPCK(1,14,11,109+LISTOFF1,IMODE,0,1)
       ENDIF
       CALL INIPCK(1,14,18,110+LISTOFF1,IMODE,0,1)
       IF(IUHF.EQ.1)THEN
        CALL INIPCK(1,1,9,127+LISTOFF2,IMODE,0,1)
        CALL INIPCK(1,2,10,128+LISTOFF2,IMODE,0,1)
        CALL INIPCK(1,13,18,129+LISTOFF2,IMODE,0,1)
       ENDIF
       CALL INIPCK(1,13,11,130+LISTOFF2,IMODE,0,1)
cYAU - initialize all gamma lists to zero
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
c&      ENDIF
C
C CREATE LISTS FOR TRIPLE CONTRIBUTION TO THE DENSITY MATRICES
C
                     CALL UPDMOI(1,NF1(1),1,160,IMODE,0)
      IF (IUHF.EQ.1) CALL UPDMOI(1,NF1(2),2,160,IMODE,0)
                     CALL UPDMOI(1,NF2(1),3,160,IMODE,0)
      IF (IUHF.EQ.1) CALL UPDMOI(1,NF2(2),4,160,IMODE,0) 
                     CALL UPDMOI(1, NT(1),5,160,IMODE,0)
      IF (IUHF.EQ.1) CALL UPDMOI(1, NT(2),6,160,IMODE,0)
      do iGrp = 1, 6, 2-iuhf
         call aces_list_memset(iGrp,160,0)
      end do

      RETURN
      END

