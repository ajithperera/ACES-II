
C CREATES LISTS USED IN EOM CALCULATIONS

#include "disttype.h"
C
      SUBROUTINE MAKE_LISTS4_PERTRB_DENS(ICORE, MAXCOR, IUHF)
C
      IMPLICIT INTEGER (A-Z)
      LOGICAL EOMCC, NODAVID, CIS, ESPROP
C
      DIMENSION ICORE(MAXCOR),MAXOO(2),MAXVV(2),MAXVO(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/    POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
C
#include "jodaflags.com"
C
      CALL ACES_IO_REMOVE(54,'DERGAM')
      INEWFIL= 0
      IZILCH = 0
      EOMCC  = .FALSE.
      NODAVID= .FALSE.
          CIS= .FALSE.
       ESPROP= .FALSE.        
C
C CALCULATE MAXIMUM SIZES OF VV, OO AND VO VECTORS
C
      MAXVV(1)=0
      MAXOO(1)=0
      MAXVO(1)=0

      DO IRREP = 1, NIRREP
         MAXVV(1)=MAX(MAXVV(1),IRPDPD(IRREP,19))
         MAXOO(1)=MAX(MAXOO(1),IRPDPD(IRREP,21))
         MAXVO(1)=MAX(MAXVO(1),IRPDPD(IRREP,9))
      END DO
C
      IF (IUHF .EQ. 1) THEN
         MAXVV(2)=0
         MAXOO(2)=0
         MAXVO(2)=0
         DO IRREP = 1, NIRREP
            MAXVV(2)=MAX(MAXVV(2),IRPDPD(IRREP,20))
            MAXOO(2)=MAX(MAXOO(2),IRPDPD(IRREP,22))
            MAXVO(2)=MAX(MAXVO(2),IRPDPD(IRREP,10))
         END DO
      ENDIF
C
C CREATE SINGLES VECTOR LISTS AND DIPOLE MOMENT LISTS
C
      DO ISPIN=1,1+IUHF
C
        CALL UPDMOI(1,MAXVO(ISPIN),ISPIN,490,INEWFIL,0)
        INEWFIL=0
CSSS        CALL UPDMOI(1,MAXVO(ISPIN),ISPIN,493,INEWFIL,0)
        CALL UPDMOI(1,MAXVO(ISPIN),ISPIN+2,490,INEWFIL,0)
CSSS        CALL UPDMOI(1,MAXVO(ISPIN),9,447+ISPIN,INEWFIL,0)
        call aces_list_memset(ispin,  490,0)
CSSS        call aces_list_memset(ispin,  493,0)
        call aces_list_memset(ispin+2,490,0)
CSSS        call aces_list_memset(9,      447+ISPIN,0)
C
#ifdef _KEEP
        DO IRREP=1,NIRREP
           LENOO=IRPDPD(IRREP,20+ISPIN) 
           LENVV=IRPDPD(IRREP,18+ISPIN) 
           LENVO=IRPDPD(IRREP,8+ISPIN) 
           IF (EOMCC)THEN
             NUMRECS=6
           ELSE
             NUMRECS=3
           ENDIF
           CALL UPDMOI(NUMRECS,LENOO,IRREP,475+ISPIN,INEWFIL,0)
           INEWFIL=0
           CALL UPDMOI(NUMRECS,LENVV,IRREP,477+ISPIN,INEWFIL,0)
           CALL UPDMOI(NUMRECS,LENVO,IRREP,479+ISPIN,INEWFIL,0)
       ENDDO
#endif
C
      ENDDO
C
#ifdef _KEEP

      IENTER=0
        IOFF=0
      IF(CIS)THEN
       call aces_io_remove(51,'GAMLAM')
      ENDIF
C
      CALL UPDMOI(3,NFMI(1),1,160,IENTER,IOFF)
      IENTER=0
        IOFF=0
      CALL UPDMOI(3,NFEA(1),3,160,IENTER,IOFF)
      CALL UPDMOI(3,NT(1),5,160,IENTER,IOFF)
C
      IF (IUHF.NE.0)THEN
         CALL UPDMOI(3,NFMI(2),2,160,IENTER,IOFF)
         CALL UPDMOI(3,NFEA(2),4,160,IENTER,IOFF)
         CALL UPDMOI(3,NT(2),6,160,IENTER,IOFF)
      ENDIF
#endif
C
C CREATE AREA FOR Q(AB) AND Q(IJ) THREE-BODY INTERMEDIATES
C
#ifdef _KEEP
      DO ISPIN=1,1+IUHF
         MAXAB=0
         MAXIJ=0
         MAXIA=0
C
         DO IRREP=1,NIRREP
            MAXAB=MAX(MAXAB,IRPDPD(IRREP,18+ISPIN))
            MAXIJ=MAX(MAXIJ,IRPDPD(IRREP,20+ISPIN))
            MAXIA=MAX(MAXIA,IRPDPD(IRREP,8+ISPIN))
         ENDDO
       
         CALL UPDMOI(1,MAXVO(ISPIN),ISPIN,490,INEWFIL,0)
         INEWFIL=0
         CALL UPDMOI(1,MAXAB,ISPIN,492,INEWFIL,0)
         CALL UPDMOI(1,MAXIJ,ISPIN,491,INEWFIL,0)
         CALL UPDMOI(1,MAXVO(ISPIN),ISPIN+2,490,INEWFIL,0)
         CALL UPDMOI(1,MAXAB,ISPIN+2,492,INEWFIL,0)
         CALL UPDMOI(1,MAXIJ,ISPIN+2,491,INEWFIL,0)
      ENDDO
#endif
C
C CALCULATE MAXIMUM SIZE OF DAVIDSON LISTS
C
#ifdef _NOSKIP
      MAXLEN=0
      DO IRREPX=1,NIRREP
         LEN=0
       DO ISPIN=1,1+IUHF
          LEN=LEN+IRPDPD(IRREPX,8+ISPIN)
       ENDDO
C
       LEN=LEN+IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
C
       IF (IUHF.NE.0)THEN
          LEN=LEN+IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
          LEN=LEN+IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
       ENDIF
       MAXLEN=MAX(MAXLEN,LEN)
      ENDDO
C
      NUMLST=2
C
      IF (NODAVID) THEN
         MAXEXP0 = 0
        NREDUCE0 = 0
      ELSE
         MAXEXP0 = MAXEXP
        NREDUCE0 = NREDUCE
      ENDIF
C
      DO I=1,NUMLST
         CALL UPDMOI(MAXEXP0,MAXLEN,I,470,0,0)
         CALL UPDMOI(MAXEXP0,MAXLEN,I,471,0,0)
         CALL UPDMOI(6+NREDUCE0,MAXLEN,I,472,0,0)
      ENDDO
#endif
C
C NOW MAKE DENOMINATOR AND T2 LISTS   
C       
      DO ISPIN=3,3-2*IUHF,-1
C
C DENOMINATOR AND T2 LISTS
C
         TTYPEL=ISYTYP(1,43+ISPIN)
         TTYPER=ISYTYP(2,43+ISPIN)
         CALL INIPCK2(1,TTYPEL,TTYPER,460+ISPIN,IZILCH,IZILCH,1)
         CALL INIPCK2(1,TTYPEL,TTYPER,447+ISPIN,IZILCH,IZILCH,1)
         CALL INIPCK2(1,TTYPEL,TTYPER,443+ISPIN,IZILCH,IZILCH,1)
C
C FOR AO-BASED ALGORITHMS, NEED AO T2 LISTS
C
         IF (IFLAGS(93).EQ.2) THEN
            CALL INIPCK2(1,TTYPEL,TTYPER,280+ISPIN,IZILCH,IZILCH,1)
            TTYPEL=15
            CALL INIPCK2(1,TTYPEL,TTYPER,213+ISPIN,IZILCH,IZILCH,1)
            CALL INIPCK2(1,TTYPEL,TTYPER,463+ISPIN,IZILCH,IZILCH,1)
         ENDIF
C
      ENDDO
C
C MAKE RESORTED R2 AND L2 LISTS
C

      CALL INIPCK2(1,h_DistType_AI,h_DistType_ai,426,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_Ai,h_DistType_Ai,428,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_AI,h_DistType_AI,434,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_AI,h_DistType_ai,437,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_Ai,h_DistType_aI,439,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_AI,h_DistType_AI,440,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_ai,h_DistType_ai,441,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_AI,h_DistType_ai,442,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_Ai,h_DistType_aI,443,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_AI,h_DistType_AI,454,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_AI,h_DistType_ai,457,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_Ai,h_DistType_aI,459,IZILCH,0,1)

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
      enddo
C
      IF (IUHF.NE.0) THEN
         CALL INIPCK2(1,h_DistType_AI,h_DistType_AI,424,IZILCH,0,1)
         CALL INIPCK2(1,h_DistType_ai,h_DistType_ai,425,IZILCH,0,1)
         CALL INIPCK2(1,h_DistType_ai,h_DistType_AI,427,IZILCH,0,1)
         CALL INIPCK2(1,h_DistType_aI,h_DistType_aI,429,IZILCH,0,1)
         CALL INIPCK2(1,h_DistType_ai,h_DistType_ai,435,IZILCH,0,1)
         CALL INIPCK2(1,h_DistType_ai,h_DistType_AI,436,IZILCH,0,1)
         CALL INIPCK2(1,h_DistType_aI,h_DistType_Ai,438,IZILCH,0,1)
         CALL INIPCK2(1,h_DistType_ai,h_DistType_ai,455,IZILCH,0,1)
         CALL INIPCK2(1,h_DistType_ai,h_DistType_AI,456,IZILCH,0,1)
         CALL INIPCK2(1,h_DistType_aI,h_DistType_Ai,458,IZILCH,0,1)

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
C
      IF (IFLAGS(91).GT.1) THEN
         CALL INIPCK2(1,h_DistType_AltB,h_DistType_IltJ,114,IZILCH,0,1)
         CALL INIPCK2(1,h_DistType_altb,h_DistType_iltj,115,IZILCH,0,1)
         CALL INIPCK2(1,h_DistType_Ab,  h_DistType_Ij,  116,IZILCH,0,1)
      END IF
C
      RETURN
      END
