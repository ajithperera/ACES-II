      SUBROUTINE AUXIOI(PRESNT,INIT,IRREPX,NLIST,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISTSZ,POP,VRT,DIRPRD
      LOGICAL PRESNT,INIT
c     PARAMETER(NLIST=12)
      DIMENSION BUF(4096)
C
      COMMON /AUXIO/ DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     1               NRECS,LUAUX
      COMMON /T3IOOF/ IJKPOS(8,8,8,2),IJKLEN(36,8,4),IJKOFF(36,8,4),
     1                NCOMB(4)
C
      COMMON /SYM/    POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /FLAGS/ IFLAGS(100)
C
#include "icdacc.com"
C
C     This routine sets up a variety of arrays for auxiliary (triples
C     and other) I/O (common blocks AUXIO and T3IOOF). It also opens
C     the file AUXMOI, and, if requested, initializes it.
C
C   --- In this particular version, all lists are of triples type ---
C   ---     There is a different version for use with vcc5q       ---
C
C     PRESNT - Specifies if AUXIOI already exists or not. If PRESNT is
C              true, we open it as old.
C     INIT   - Specifies if we must initialize AUXIOI or not. If INIT
C              is true, we initialize it.
C     IRREPX - Overall symmetry of lists. This is used to calculate
C              starting words of parts of lists.
C     NLIST  - Number of triples lists to be created. Note this variable
C              should be specified as though as UHF calculation is
C              being performed. Logic in this routine reduces disk space.
C              Also, there is a reduction automatically handled in
C              CCSD based noniterative calculations. For normal energy
C              calculations (e.g. CCSDT), NLIST is 8. For EOM-CCSDT-1
C              property calculations it is 12. For EOM-CCSDT property
C              calculations it will be 16 (for two sets of T and two sets
C              of lambda amplitudes.
C
C     N.B. There is some logical interplay between PRESNT and INIT. Thus,
C          if PRESNT is false, it is likely that INIT should be true,
C          since we may end up trying to read from a file which has not
C          been written.
C
C     LNPHYR : Length (in double precision words) of physical records.
C     LUAUX  : Unit number for T3 file (or other auxiliary file).
C
      LNPHYR = 4096
      LUAUX  = 91
C
      CALL IZERO(DISTSZ,8 * 100)
      CALL IZERO(NDISTS,8 * 100)
      CALL IZERO(INIWRD,8 * 100)
C
C     Calculate distribution sizes and numbers of distributions (blocks of
C     4).
C
      NBLOCK = NLIST / 4
      DO 10 IBLOCK=1,NBLOCK
      IPOS = (IBLOCK - 1) * 4 + 1
      CALL DIMT3(DISTSZ(1,IPOS),NDISTS(1,IPOS))
   10 CONTINUE
C
C     Initial words of symmetry sublists. Note that this logic is
C     supposed to let us handle lists which are not totally symmetric.
C     We cannot currently have amplitudes of different symmetries
C     around. Think about that for possible future things.
C
      IF(IUHF.EQ.0)THEN
C
      DO   30 ILIST =1,NLIST
      DO   20 IRREPR=1,NIRREP
C
      IRREPL = DIRPRD(IRREPR,IRREPX)
C
      IF(ILIST.EQ.1.AND.IRREPR.EQ.1)THEN
      INIWRD(IRREPR,ILIST) = 1
      ENDIF
C
      IF(               IRREPR.GT.1)THEN
        if(ilist.eq.3.or.ilist.eq.4.or.ilist.eq.7.or.ilist.eq.8)then
          nsize = 0
        else
          NSIZE = DISTSZ(DIRPRD(IRREPR-1,IRREPX),ILIST) * 
     1            NDISTS(       IRREPR-1        ,ILIST)
        endif
      INIWRD(IRREPR,ILIST) = INIWRD(IRREPR-1,ILIST) + NSIZE
      ENDIF
C
      IF(ILIST.GT.1.AND.IRREPR.EQ.1)THEN
        if(ilist.eq.4.or.ilist.eq.5.or.ilist.eq.8)then
          nsize = 0
        else
          NSIZE = DISTSZ(DIRPRD(NIRREP  ,IRREPX),ILIST-1) *
     1            NDISTS(       NIRREP          ,ILIST-1)
        endif
      INIWRD(IRREPR,ILIST) = INIWRD(NIRREP,ILIST-1) + NSIZE
      ENDIF
C
   20 CONTINUE
   30 CONTINUE
C
      ELSE
C
      DO  130 ILIST =1,NLIST
      DO  120 IRREPR=1,NIRREP
C
      IRREPL = DIRPRD(IRREPR,IRREPX)
C
      IF(ILIST.EQ.1.AND.IRREPR.EQ.1)THEN
      INIWRD(IRREPR,ILIST) = 1
      ENDIF
C
      IF(               IRREPR.GT.1)THEN
      NSIZE = DISTSZ(DIRPRD(IRREPR-1,IRREPX),ILIST) * 
     1        NDISTS(       IRREPR-1        ,ILIST)
      INIWRD(IRREPR,ILIST) = INIWRD(IRREPR-1,ILIST) + NSIZE
      ENDIF
C
      IF(ILIST.GT.1.AND.IRREPR.EQ.1)THEN
      NSIZE = DISTSZ(DIRPRD(NIRREP  ,IRREPX),ILIST-1) *
     1        NDISTS(       NIRREP          ,ILIST-1)
      INIWRD(IRREPR,ILIST) = INIWRD(NIRREP,ILIST-1) + NSIZE
      ENDIF
C
  120 CONTINUE
  130 CONTINUE
C
      ENDIF
C
C     For certain noniterative triples calculations make addresses of
C     5-8 the same as 1-4.
C
      IF(IFLAGS(2).EQ.12.OR.IFLAGS(2).EQ.26.OR.IFLAGS(2).EQ.27.OR.
     &   IFLAGS(2).EQ.28.OR.IFLAGS(2).EQ.31)THEN
        CALL ICOPY(NIRREP,INIWRD(1,1),1,INIWRD(1,5),1)
        CALL ICOPY(NIRREP,INIWRD(1,2),1,INIWRD(1,6),1)
        CALL ICOPY(NIRREP,INIWRD(1,3),1,INIWRD(1,7),1)
        CALL ICOPY(NIRREP,INIWRD(1,4),1,INIWRD(1,8),1)
      ENDIF
C
      NWDTOT = INIWRD(NIRREP,NLIST) + 
     1         DISTSZ(DIRPRD(NIRREP,IRREPX),NLIST) *
     1         NDISTS(       NIRREP        ,NLIST)
      WRITE(6,1010) NWDTOT
 1010 FORMAT(' @AUXIOI-I, Total number of words on AUXMOI ',I15)
C
C     Open the file(s).
C
      IF(PRESNT)THEN
c Nevin 8-31-95 modified to account for direct access in words or bytes
c      OPEN(UNIT=LUAUX,STATUS='UNKNOWN',FILE='AUXMOI',ACCESS='DIRECT',
c    1      FORM='UNFORMATTED',RECL=LNPHYR * 8)
       OPEN(UNIT=LUAUX,STATUS='UNKNOWN',FILE='AUXMOI',ACCESS='DIRECT',
     1      FORM='UNFORMATTED',RECL=LNPHYR * iintfp * idaccm)
      ELSE
       OPEN(UNIT=LUAUX,STATUS='NEW',FILE='AUXMOI',ACCESS='DIRECT',
     1      FORM='UNFORMATTED',RECL=LNPHYR * iintfp * idaccm)
      ENDIF
C
      IF(INIT)THEN
C
C     Initialize the lists.
C
      NRECS = NWDTOT/LNPHYR + 1
      CALL ZERO(BUF,LNPHYR)
      DO   50 IREC=1,NRECS
      WRITE(LUAUX,REC=IREC) BUF
   50 CONTINUE
C
      ELSE
C
C     Do not initialize, just calculate number of records.
C
      NRECS = NWDTOT/LNPHYR + 1
C
      ENDIF
C
C     Calculate pointers, lengths, and offsets for T3 i/o.
C
      CALL OFFT3
      RETURN
      END
