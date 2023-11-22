      SUBROUTINE AUXIOI(PRESNT,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISTSZ,POP,VRT,DIRPRD
      LOGICAL PRESNT
      PARAMETER(NLIST=8)
      DIMENSION BUF(4096)
C
      COMMON /AUXIO/ DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     1               NRECS,LUAUX
      COMMON /T3IOOF/ IJKPOS(8,8,8,2),IJKLEN(36,8,4),IJKOFF(36,8,4),
     1                NCOMB(4)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      COMMON /SYM/    POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
 
#include "icdacc.com"

      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /FLAGS/  IFLAGS(100)
C
C     If PRESNT is .TRUE., then AUXMOI already exists and we open it
C     as 'OLD' and do not initialize the contents. If PRESNT is .FALSE.,
C     AUXMOI is not there and we initialize it.
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
C     Calculate distribution sizes and numbers of distributions.
C
      CALL DIMT3(DISTSZ(1,1),NDISTS(1,1))
      CALL DIMT3(DISTSZ(1,5),NDISTS(1,5))
C
C     Initial words of symmetry sublists. Note that this logic assumes
C     that the lists are overall total symmetric. This will change if
C     we do EOM or Fock space with triples.
C
      IF(IUHF.EQ.0)THEN
C
      DO   30 ILIST =1,NLIST
      DO   20 IRREP =1,NIRREP
C
      IF(ILIST.EQ.1.AND.IRREP.EQ.1)THEN
      INIWRD(IRREP,ILIST) = 1
      ENDIF
C
      IF(               IRREP.GT.1)THEN
        if(ilist.eq.3.or.ilist.eq.4.or.ilist.eq.7.or.ilist.eq.8)then
          nsize = 0
        else
          NSIZE = DISTSZ(IRREP-1,ILIST) * NDISTS(IRREP-1,ILIST)
        endif
      INIWRD(IRREP,ILIST) = INIWRD(IRREP-1,ILIST) + NSIZE
      ENDIF
C
      IF(ILIST.GT.1.AND.IRREP.EQ.1)THEN
        if(ilist.eq.4.or.ilist.eq.5.or.ilist.eq.8)then
          nsize = 0
        else
          NSIZE = DISTSZ(NIRREP,ILIST-1) * NDISTS(NIRREP,ILIST-1)
        endif
      INIWRD(IRREP,ILIST) = INIWRD(NIRREP,ILIST-1) + NSIZE
      ENDIF
C
   20 CONTINUE
   30 CONTINUE
C
      ELSE
C
      DO  130 ILIST =1,NLIST
      DO  120 IRREP =1,NIRREP
C
      IF(ILIST.EQ.1.AND.IRREP.EQ.1)THEN
      INIWRD(IRREP,ILIST) = 1
      ENDIF
C
      IF(               IRREP.GT.1)THEN
      NSIZE = DISTSZ(IRREP-1,ILIST) * NDISTS(IRREP-1,ILIST)
      INIWRD(IRREP,ILIST) = INIWRD(IRREP-1,ILIST) + NSIZE
      ENDIF
C
      IF(ILIST.GT.1.AND.IRREP.EQ.1)THEN
      NSIZE = DISTSZ(NIRREP,ILIST-1) * NDISTS(NIRREP,ILIST-1)
      INIWRD(IRREP,ILIST) = INIWRD(NIRREP,ILIST-1) + NSIZE
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
      IF(IFLAGS(2).EQ.12.OR.(IFLAGS(2).GE.26.AND.IFLAGS(2).LE.31))THEN
      CALL ICOPY(NIRREP,INIWRD(1,1),1,INIWRD(1,5),1)
      CALL ICOPY(NIRREP,INIWRD(1,2),1,INIWRD(1,6),1)
      CALL ICOPY(NIRREP,INIWRD(1,3),1,INIWRD(1,7),1)
      CALL ICOPY(NIRREP,INIWRD(1,4),1,INIWRD(1,8),1)
      ENDIF
C
      NWDTOT = INIWRD(NIRREP,NLIST) + 
     1         DISTSZ(NIRREP,NLIST) * NDISTS(NIRREP,NLIST)
      WRITE(6,1010) NWDTOT
 1010 FORMAT(' @AUXIOI-I, Total number of words on AUXMOI ',I15)
C
C     Open the file(s).
C
      IF(PRESNT)THEN
      OPEN(UNIT=LUAUX,STATUS='OLD',FILE='AUXMOI',ACCESS='DIRECT',
     1     FORM='UNFORMATTED',RECL=LNPHYR * idaccm * iintfp)
      NRECS = NWDTOT/LNPHYR + 1
      ELSE
c YAU - old
c      OPEN(UNIT=LUAUX,STATUS='NEW',FILE='AUXMOI',ACCESS='DIRECT',
c     1     FORM='UNFORMATTED',RECL=LNPHYR * 8)
c YAU - new
      OPEN(UNIT=LUAUX,STATUS='NEW',FILE='AUXMOI',ACCESS='DIRECT',
     1     FORM='UNFORMATTED',RECL=LNPHYR * idaccm * iintfp)
c YAU - end
C
C     Initialize the file(s).
C
      NRECS = NWDTOT/LNPHYR + 1
      CALL ZERO(BUF,LNPHYR)
      DO   40 IREC=1,NRECS
      WRITE(LUAUX,REC=IREC) BUF
   40 CONTINUE
C
      ENDIF
C
C     Calculate pointers, lengths, and offsets for T3 i/o.
C
      CALL OFFT3
C
      RETURN
      END
