
      SUBROUTINE MAKLST(ICORE,MAXCOR,IUHF)
C
C CREATES LISTS USED IN EOM CALCULATIONS
C
      IMPLICIT INTEGER (A-Z)
      LOGICAL CIS,EOMCC,FULDIAG,INCORE,READGUES,ESPROP
      LOGICAL CCSD, MBPT, PARTEOM, NODAVID
      DIMENSION ICORE(MAXCOR),MAXOO(2),MAXVV(2),MAXVO(2)
      COMMON/METH/CIS,EOMCC,FULDIAG,INCORE,READGUES
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /EXTRAP/ MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON /SYM/    POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /PROPGRAD/ ESPROP,IOPTROOT,IOPTSYM
      COMMON /FLAGS/ IFLAGS(100)
      COMMON/EOMINFO/CCSD, MBPT, PARTEOM, NODAVID
      COMMON/LISTDENS/LDENS
C
c      INEWFIL=1
      call aces_io_remove(54,'DERGAM')
      INEWFIL=0
C
C CALCULATE MAXIMUM SIZES OF VV, OO AND VO VECTORS
C
      DO 10 ISPIN=1,1+IUHF
         MAXVV(ISPIN)=0
         MAXOO(ISPIN)=0
         MAXVO(ISPIN)=0
         DO 20 IRREP=1,NIRREP
            MAXVV(ISPIN)=MAX(MAXVV(ISPIN),IRPDPD(IRREP,18+ISPIN))
            MAXOO(ISPIN)=MAX(MAXOO(ISPIN),IRPDPD(IRREP,20+ISPIN))
            MAXVO(ISPIN)=MAX(MAXVO(ISPIN),IRPDPD(IRREP,8+ISPIN))
 20      CONTINUE
 10   CONTINUE
C
C CREATE SINGLES VECTOR LISTS AND DIPOLE MOMENT LISTS
C
      DO 11 ISPIN=1,1+IUHF
         CALL UPDMOI(1,MAXVO(ISPIN),ISPIN,490,INEWFIL,0)
         INEWFIL=0
         CALL UPDMOI(1,MAXVO(ISPIN),ISPIN,493,INEWFIL,0)
         CALL UPDMOI(1,MAXVO(ISPIN),ISPIN+2,490,INEWFIL,0)
         CALL UPDMOI(1,MAXVO(ISPIN),9,447+ISPIN,INEWFIL,0)
         call aces_list_memset(ispin,490,0)
         call aces_list_memset(ispin,493,0)
         call aces_list_memset(ispin+2,490,0)
         call aces_list_memset(9,447+ISPIN,0)
 11   CONTINUE
C
      IENTER=0
      IOFF=0
      CALL UPDMOI(3,NFMI(1),1,ldens,IENTER,IOFF)
      IENTER=0
      IOFF=0
      CALL UPDMOI(3,NFEA(1),3,ldens,IENTER,IOFF)
      CALL UPDMOI(3,NT(1),5,ldens,IENTER,IOFF)
      IF(IUHF.NE.0)THEN
         CALL UPDMOI(3,NFMI(2),2,ldens,IENTER,IOFF)
         CALL UPDMOI(3,NFEA(2),4,ldens,IENTER,IOFF)
         CALL UPDMOI(3,NT(2),6,ldens,IENTER,IOFF)
      ENDIF
C
C CREATE AREA FOR Q(AB) AND Q(IJ) THREE-BODY INTERMEDIATES
C
      DO 110 ISPIN=1,1+IUHF
         MAXAB=0
         MAXIJ=0
         DO 120 IRREP=1,NIRREP
            MAXAB=MAX(MAXAB,IRPDPD(IRREP,18+ISPIN))
            MAXIJ=MAX(MAXIJ,IRPDPD(IRREP,20+ISPIN))
 120     CONTINUE
         CALL UPDMOI(1,MAXAB,ISPIN,492,INEWFIL,0)
         CALL UPDMOI(1,MAXIJ,ISPIN,491,INEWFIL,0)
C esprop
         CALL UPDMOI(1,MAXAB,ISPIN+2,492,INEWFIL,0)
         CALL UPDMOI(1,MAXIJ,ISPIN+2,491,INEWFIL,0)
 110  CONTINUE
C
C CALCULATE MAXIMUM SIZE OF DAVIDSON LISTS
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
      MAXEXP=IFLAGS(97)
      NREDUCE = 0
      DO 27 I=1,NUMLST
        IF (.NOT. NODAVID) THEN
          CALL UPDMOI(MAXEXP,MAXLEN,I,470,0,0)
          CALL UPDMOI(MAXEXP,MAXLEN,I,471,0,0)
        ENDIF
         CALL UPDMOI(3+NREDUCE,MAXLEN,I,472,0,0)
 27   CONTINUE
C
C NOW MAKE DENOMINATOR AND T2 LISTS
C
      DO 30 ISPIN=3,3-2*IUHF,-1
C
C DENOMINATOR AND T2 LISTS
C
         TTYPEL=ISYTYP(1,43+ISPIN)
         TTYPER=ISYTYP(2,43+ISPIN)
         CALL INIPCK2(1,TTYPEL,TTYPER,460+ISPIN,0,0,1)
         CALL INIPCK2(1,TTYPEL,TTYPER,447+ISPIN,0,0,1)
         CALL INIPCK2(1,TTYPEL,TTYPER,443+ISPIN,0,0,1)
         CALL INIPCK2(1,TTYPEL,TTYPER,453+ISPIN,0,0,1)
C
C FOR AO-BASED ALGORITHMS, NEED AO T2 LISTS
C
CSSS         IF(IFLAGS(93).EQ.2)THEN
CSSS            TTYPEL=15
CSSS            CALL INIPCK2(1,TTYPEL,TTYPER,413+ISPIN,0,0,1)
CSSS            CALL INIPCK2(1,TTYPEL,TTYPER,463+ISPIN,0,0,1)
CSSS         ENDIF

C FOR AO-BASED ALGORITHMS, NEED AO T2 LISTS
C
       IF(IFLAGS(93).EQ.2)THEN
           CALL INIPCK2(1,TTYPEL,TTYPER,280+ISPIN,IZILCH,IZILCH,1)
           TTYPEL=15
           CALL INIPCK2(1,TTYPEL,TTYPER,213+ISPIN,IZILCH,IZILCH,1)
           CALL INIPCK2(1,TTYPEL,TTYPER,463+ISPIN,IZILCH,IZILCH,1)
       ENDIF
C
 30   CONTINUE
C
C MAKE RESORTED T2 LIST
C
      if (IUHF.eq.0) then
         DO ILIST=37,43,2
            TTYPEL=ISYTYP(1,ILIST)
            TTYPER=ISYTYP(2,ILIST)
            CALL INIPCK2(1,TTYPEL,TTYPER,400+ILIST,0,0,1)
            do iGrp = 1, nirrep
               call aces_list_memset(iGrp,400+ILIST,0)
            end do
         END DO
            CALL INIPCK2(1,ISYTYP(1,42),ISYTYP(2,42),400+42,0,0,1)
            do iGrp = 1, nirrep
               call aces_list_memset(iGrp,400+42,0)
            end do
      else
         DO ILIST=34,43
            TTYPEL=ISYTYP(1,ILIST)
            TTYPER=ISYTYP(2,ILIST)
            CALL INIPCK2(1,TTYPEL,TTYPER,400+ILIST,0,0,1)
            do iGrp = 1, nirrep
               call aces_list_memset(iGrp,400+ILIST,0)
            end do
         END DO
      end if
C
      RETURN
      END