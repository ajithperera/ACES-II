
C CREATES LISTS USED IN EOM CALCULATIONS

#include "disttype.h"

      SUBROUTINE MAKLST2(LISTZ1,LISTZ1OFF,LISTZ2,IUHF)
      IMPLICIT INTEGER (A-Z)
      integer iarg1
      DIMENSION MAXOO(2),MAXVV(2),MAXVO(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /EXTRAP/ MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON /SYM/    POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /FLAGS/ IFLAGS(100)

      iarg1=0
      call aces_io_remove(54,'DERGAM')
C
C     CALCULATE MAXIMUM SIZES OF VV, OO AND VO VECTORS
C
      MAXVV(1)=0
      MAXOO(1)=0
      MAXVO(1)=0
      DO IRREP = 1, NIRREP
         MAXVV(1)=MAX(MAXVV(1),IRPDPD(IRREP,19))
         MAXOO(1)=MAX(MAXOO(1),IRPDPD(IRREP,21))
         MAXVO(1)=MAX(MAXVO(1),IRPDPD(IRREP,9))
      END DO
      if (iuhf.eq.1) then
         MAXVV(2)=0
         MAXOO(2)=0
         MAXVO(2)=0
         DO IRREP = 1, NIRREP
            MAXVV(2)=MAX(MAXVV(2),IRPDPD(IRREP,20))
            MAXOO(2)=MAX(MAXOO(2),IRPDPD(IRREP,22))
            MAXVO(2)=MAX(MAXVO(2),IRPDPD(IRREP,10))
         END DO
      end if
C
C     CREATE SINGLES VECTOR LISTS AND DIPOLE MOMENT LISTS
C
      DO 11 ISPIN=1,1+IUHF
         CALL UPDMOI(1,MAXVO(ISPIN),ISPIN,490,0,0)
         CALL UPDMOI(1,MAXVO(ISPIN),ISPIN+LISTZ1OFF,LISTZ1,0,0)
         CALL UPDMOI(1,MAXVO(ISPIN),ISPIN,493,0,0)
         CALL UPDMOI(1,MAXVO(ISPIN),ISPIN+2,490,0,0)
         CALL UPDMOI(1,MAXVO(ISPIN),9,447+ISPIN,0,0)
         call aces_list_memset(ispin,490,0)
         call aces_list_memset(ispin+listz1off,LISTZ1,0)
         call aces_list_memset(ispin,493,0)
         call aces_list_memset(ispin+2,490,0)
         call aces_list_memset(9,447+ISPIN,0)
C
         DO 17 IRREP=1,NIRREP
            LENOO=IRPDPD(IRREP,20+ISPIN)
            LENVV=IRPDPD(IRREP,18+ISPIN)
            LENVO=IRPDPD(IRREP,8+ISPIN)
            NUMRECS=6
            CALL UPDMOI(NUMRECS,LENOO,IRREP,475+ISPIN,0,0)
            CALL UPDMOI(NUMRECS,LENVV,IRREP,477+ISPIN,0,0)
            CALL UPDMOI(NUMRECS,LENVO,IRREP,479+ISPIN,0,0)
 17      CONTINUE
 11   CONTINUE

      IENTER=0
      IOFF=0
C
C     CREATE AREA FOR Q(AB) AND Q(IJ) THREE-BODY INTERMEDIATES
C
      DO 110 ISPIN=1,1+IUHF
         MAXAB=0
         MAXIJ=0
         MAXIA=0
         DO 120 IRREP=1,NIRREP
            MAXAB=MAX(MAXAB,IRPDPD(IRREP,18+ISPIN))
            MAXIJ=MAX(MAXIJ,IRPDPD(IRREP,20+ISPIN))
            MAXIA=MAX(MAXIA,IRPDPD(IRREP,8+ISPIN))
 120     CONTINUE
         CALL UPDMOI(1,MAXVO(ISPIN),ISPIN,490,0,0)
         CALL UPDMOI(1,MAXAB,ISPIN,492,0,0)
         CALL UPDMOI(1,MAXIJ,ISPIN,491,0,0)
         CALL UPDMOI(1,MAXVO(ISPIN),ISPIN+2,490,0,0)
         CALL UPDMOI(1,MAXAB,ISPIN+2,492,0,0)
         CALL UPDMOI(1,MAXIJ,ISPIN+2,491,0,0)
 110  CONTINUE
C
C     CALCULATE MAXIMUM SIZE OF DAVIDSON LISTS
C
      MAXLEN=0
      DO 25 IRREPX=1,NIRREP
         LEN=0
         DO 13 ISPIN=1,1+IUHF
            LEN=LEN+IRPDPD(IRREPX,8+ISPIN)
 13      CONTINUE
         LEN=LEN+IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
         IF(IUHF.NE.0)THEN
            LEN=LEN+IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
            LEN=LEN+IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
         ENDIF
         MAXLEN=MAX(MAXLEN,LEN)
 25   CONTINUE
      NUMLST=2
      MAXEXP0 = MAXEXP
      NREDUCE0 = NREDUCE
      DO 27 I=1,NUMLST
         CALL UPDMOI(MAXEXP0,MAXLEN,I,470,0,0)
         CALL UPDMOI(MAXEXP0,MAXLEN,I,471,0,0)
         CALL UPDMOI(6+NREDUCE0,MAXLEN,I,472,0,0)
 27   CONTINUE
C
C     NOW MAKE DENOMINATOR AND T2 LISTS
C
      DO 30 ISPIN=3,3-2*IUHF,-1
C
C DENOMINATOR AND T2 LISTS
C
         TTYPEL=ISYTYP(1,43+ISPIN)
         TTYPER=ISYTYP(2,43+ISPIN)
         CALL INIPCK2(1,TTYPEL,TTYPER,460+ISPIN,iarg1,0,1)
         CALL INIPCK2(1,TTYPEL,TTYPER,LISTZ2-1+ISPIN,iarg1,0,1)
         CALL INIPCK2(1,TTYPEL,TTYPER,447+ISPIN,iarg1,0,1)
         CALL INIPCK2(1,TTYPEL,TTYPER,443+ISPIN,iarg1,0,1)
C
C     FOR AO-BASED ALGORITHMS, NEED AO T2 LISTS
C
         IF(IFLAGS(93).EQ.2)THEN
            CALL INIPCK2(1,TTYPEL,TTYPER,280+ISPIN,iarg1,0,1)
            TTYPEL=15
            CALL INIPCK2(1,TTYPEL,TTYPER,213+ISPIN,iarg1,0,1)
            CALL INIPCK2(1,TTYPEL,TTYPER,463+ISPIN,iarg1,0,1)
         ENDIF
 30   CONTINUE

C
C     MAKE RESORTED R2 AND L2 LISTS
C
      CALL INIPCK2(1,h_DistType_AI,h_DistType_ai,426,iarg1,0,1)
      CALL INIPCK2(1,h_DistType_Ai,h_DistType_Ai,428,iarg1,0,1)
      CALL INIPCK2(1,h_DistType_AI,h_DistType_AI,434,iarg1,0,1)
      CALL INIPCK2(1,h_DistType_AI,h_DistType_ai,437,iarg1,0,1)
      CALL INIPCK2(1,h_DistType_Ai,h_DistType_aI,439,iarg1,0,1)
      CALL INIPCK2(1,h_DistType_AI,h_DistType_AI,440,iarg1,0,1)
      CALL INIPCK2(1,h_DistType_ai,h_DistType_ai,441,iarg1,0,1)
      CALL INIPCK2(1,h_DistType_AI,h_DistType_ai,442,iarg1,0,1)
      CALL INIPCK2(1,h_DistType_Ai,h_DistType_aI,443,iarg1,0,1)
      CALL INIPCK2(1,h_DistType_AI,h_DistType_AI,454,iarg1,0,1)
      CALL INIPCK2(1,h_DistType_AI,h_DistType_ai,457,iarg1,0,1)
      CALL INIPCK2(1,h_DistType_Ai,h_DistType_aI,459,iarg1,0,1)
      do iGrp = 1, nirrep
         call aces_list_memset(iGrp,426,0)
         call aces_list_memset(iGrp,428,0)
         call aces_list_memset(iGrp,434,0)
         call aces_list_memset(iGrp,437,0)
         call aces_list_memset(iGrp,439,0)
         call aces_list_memset(iGrp,440,0)
         call aces_list_memset(iGrp,441,0)
         call aces_list_memset(iGrp,442,0)
         call aces_list_memset(iGrp,443,0)
         call aces_list_memset(iGrp,454,0)
         call aces_list_memset(iGrp,457,0)
         call aces_list_memset(iGrp,459,0)
      end do
      IF (IUHF.NE.0) THEN
         CALL INIPCK2(1,h_DistType_AI,h_DistType_AI,424,iarg1,0,1)
         CALL INIPCK2(1,h_DistType_ai,h_DistType_ai,425,iarg1,0,1)
         CALL INIPCK2(1,h_DistType_ai,h_DistType_AI,427,iarg1,0,1)
         CALL INIPCK2(1,h_DistType_aI,h_DistType_aI,429,iarg1,0,1)
         CALL INIPCK2(1,h_DistType_ai,h_DistType_ai,435,iarg1,0,1)
         CALL INIPCK2(1,h_DistType_ai,h_DistType_AI,436,iarg1,0,1)
         CALL INIPCK2(1,h_DistType_aI,h_DistType_Ai,438,iarg1,0,1)
         CALL INIPCK2(1,h_DistType_ai,h_DistType_ai,455,iarg1,0,1)
         CALL INIPCK2(1,h_DistType_ai,h_DistType_AI,456,iarg1,0,1)
         CALL INIPCK2(1,h_DistType_aI,h_DistType_Ai,458,iarg1,0,1)
         do iGrp = 1, nirrep
            call aces_list_memset(iGrp,424,0)
            call aces_list_memset(iGrp,425,0)
            call aces_list_memset(iGrp,427,0)
            call aces_list_memset(iGrp,429,0)
            call aces_list_memset(iGrp,435,0)
            call aces_list_memset(iGrp,436,0)
            call aces_list_memset(iGrp,438,0)
            call aces_list_memset(iGrp,455,0)
            call aces_list_memset(iGrp,456,0)
            call aces_list_memset(iGrp,458,0)
         end do
      END IF

      RETURN
      END
