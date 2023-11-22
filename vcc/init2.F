C
C   THIS ROUTINE INITIALIZES THE T2 LISTS ON THE 
C   GAMLAM FILE WITHIN ROHF-MBPT(3), MBPT(4),
C   QCISD, AND CCSD GRADIENT CALCULATIONS
C
C CODED AUGUST/90, JUNE/91 JG
C
      SUBROUTINE INIT2(IUHF,PCCD)
      IMPLICIT INTEGER(A-Z)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      LOGICAL TRIPNI,TRIPNI1,TRIPIT,T3STOR,NONHF,YESNO
      LOGICAL PCCD
      CHARACTER*80 FNAME
      COMMON/FLAGS/IFLAGS(100)
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON/TRIPLES/TRIPNI,TRIPNI1,TRIPIT,T3STOR
      COMMON/NHFREF/NONHF
      COMMON /SYM/ POP1(8),POP2(8),VRT1(8),VRT2(8),NT(2),
     &             NF1(2),NF2(2)
#include "syminf.com" /* for nirrep */
C
C-----------------------------------------------------------------------
C     See if GAMLAM is there. If so we create list 146 with IMODE=0.
C     This is important in triples EOM and some kinds of triples deriv-
C     atives. JDW 3/26/96.
C-----------------------------------------------------------------------
C
c      CALL GFNAME('GAMLAM  ',FNAME,ILENGTH)
c      INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=YESNO)
c      IF(YESNO)THEN
c        IMODE = 0
c      ELSE
c        IMODE=1
c      ENDIF
C
C CREATE LISTS FOR DOUBLES IN GRADIENT CALCULATIONS
C
      IMODE=0
      CALL INIPCK(1,13,14,146,IMODE,0,1)
      CALL INIPCK(1, 1, 3,144,IMODE,0,1)
      IF (IUHF.EQ.1) CALL INIPCK(1,2,4,145,IMODE,0,1)
C
C CREATE LISTS FOR SINGLES IN GRADIENT CALCULATIONS
C
      IF(M4SDQ.OR.QCISD.OR.CCSD.OR.PCCD.OR.(NONHF.AND. 
     &                            (MBPT3) .OR.CCD))THEN
       DO 20 ISPIN=1,IUHF+1
        CALL UPDMOI(1,NT(ISPIN),ISPIN,190,0,0)
        CALL UPDMOI(1,NT(ISPIN),ISPIN+2,190,0,0)
        call aces_list_memset(ispin,  190,0)
        call aces_list_memset(ispin+2,190,0)
20     CONTINUE
      ENDIF

      IF (IFLAGS(87).EQ.3.AND.(TRIPIT.OR.TRIPNI)) RETURN
C
C CREATE LISTS FOR TRIPLES IN GRADIENT CALCULATIONS
C
C  ACTUALLY LISTS G(AB;CI) AND G(IJ;KA)
C
      IF (TRIPNI) THEN

         IF (CCSD) THEN
            LISTOFF1=55
            LISTOFF2=39
         ELSE
            LISTOFF2=0
            LISTOFF1=0
         END IF
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
C The list numbers, 101-104 is uese eleswhere with different dimensions
C (see cphf and vdint MES). 04/2008, Ajith Perera
C
         IF (CCSD.OR.QCISD) THEN
            IF (IUHF.EQ.1) THEN
               CALL INIPCK(1,1,3,1,IMODE,0,1)
               CALL INIPCK(1,2,4,2,IMODE,0,1)
               do iFam = 1, 2
               do iGrp = 1, nirrep
                  call aces_list_memset(iGrp,iFam,0)
               end do
               end do
            END IF
            CALL INIPCK(1,13,14,3,IMODE,0,1)
            do iGrp = 1, nirrep
               call aces_list_memset(iGrp,3,0)
            end do
            IF (TRIPNI1) THEN
               CALL UPDMOI(1,NT(1),1,4,IMODE,0)
               IF (IUHF.EQ.1) THEN
                  CALL UPDMOI(1,NT(2),2,4,IMODE,0)
               END IF
               call aces_list_memset(1,4,0)
               IF (IUHF.EQ.1) THEN
                  call aces_list_memset(2,4,0)
               END IF
            END IF
            IF (IUHF.EQ.1) THEN
               CALL INIPCK(1,1,3,114,IMODE,0,1)
               CALL INIPCK(1,2,4,115,IMODE,0,1)
               do iFam = 114, 115
               do iGrp = 1, nirrep
                  call aces_list_memset(iGrp,iFam,0)
               end do
               end do
            END IF
            CALL INIPCK(1,13,14,116,IMODE,0,1)
            do iGrp = 1, nirrep
               call aces_list_memset(iGrp,116,0)
            end do
         END IF

c     END IF (TRIPNI)
      END IF

      RETURN
      END

