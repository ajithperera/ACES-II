      SUBROUTINE GETLIST(Z,INIDIS,NDIS,ICACHE,IRREPR,ILIST,IRREPX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISTSZ,DIRPRD
      LOGICAL INCORE
#include "trps.par" /* for NUMT3 */
      DIMENSION Z(1),BUF(4096)
C
      COMMON /AUXIO/ DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     1               NUMRECS,LUAUX
      COMMON /T3FILE/ ALLT3(NUMT3)
      COMMON /T3CORE/ INCORE
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
C
      IRREPL = DIRPRD(IRREPX,IRREPR)
C
C     Number of words to be read.
C
      NWORDS = NDIS * DISTSZ(IRREPL,ILIST)
c      write(6,*) ' nwords ',nwords
C
      IF(NWORDS.LE.0) RETURN
C
C     First word on the file to be read.
C
      NFIRST = INIWRD(IRREPR,ILIST) + (INIDIS-1)*DISTSZ(IRREPL,ILIST)
C
C     Try to read from ALLT3. All we need to know is NWORDS and NFIRST
C     (apparently).
C
      IF(INCORE)THEN
      CALL DCOPY(NWORDS,ALLT3(NFIRST),1,Z,1)
      RETURN
      ENDIF
C
C     First record to be read.
C
      INIREC = NFIRST / LNPHYR
C      
      NCOMP = LNPHYR * INIREC
      IF(NFIRST.GT.NCOMP)THEN
      INIREC = INIREC + 1
      ENDIF
C
C     First word in INIREC which is needed (preceding ones not used
C     at this time).
C
      NUSE1 = NFIRST - (INIREC-1) * LNPHYR
C
C     Number of remaining words in INIREC.
C
      NREM = LNPHYR - NUSE1 + 1
C
C     Are all data on INIREC ? If so read them and return.
C
c      write(6,*) inirec,nuse1,nwords,nrem,lnphyr
      IF(NREM.GE.NWORDS)THEN
      READ(LUAUX,REC=INIREC) BUF
C
      DO   10 IWORD=1,NWORDS
      Z(IWORD) = BUF(NUSE1 + IWORD - 1)
   10 CONTINUE
C
      RETURN
      ENDIF
C
C     If not read the first record and put everything beyond NUSE1 in Z.
C
      READ(LUAUX,REC=INIREC) BUF
C
      DO   20 IWORD=1,NREM
      Z(IWORD) = BUF(NUSE1 + IWORD - 1)
   20 CONTINUE
C
C     Now deal with remaining words to be read.
C
      NLEFT = NWORDS - NREM
      NRECS = NLEFT / LNPHYR
      NCOMP = LNPHYR * NRECS
      IF(NCOMP.LT.NLEFT)THEN
      NRECS = NRECS + 1
C
C     We have NRECS-1 records of LNPHYR words plus one with less than
C     LNPHYR words to be read.
C
      IF(NRECS.GT.1)THEN
C
      DO   40 IREC =1,NRECS-1
C     READ(LUAUX,REC=INIREC + IREC -1) BUF
      READ(LUAUX,REC=INIREC + IREC   ) BUF
      DO   30 IWORD=1,LNPHYR
      Z(NREM + (IREC-1)*LNPHYR + IWORD) = BUF(IWORD)
   30 CONTINUE
   40 CONTINUE
C
      ENDIF
C
C     Deal with contribution from last record.
C
      NLAST = NLEFT - (NRECS-1) * LNPHYR
C     READ(LUAUX,REC=INIREC + NRECS - 1) BUF
      READ(LUAUX,REC=INIREC + NRECS    ) BUF
C
      DO   50 IWORD=1,NLAST
      Z(NREM + (NRECS-1)*LNPHYR + IWORD) = BUF(IWORD)
   50 CONTINUE
C
      ELSE
C
C     Distributions end at the end of a physical record (ie NLEFT is an
C     integer multiple of LNPHYR).
C
      DO   70 IREC =1,NRECS
C     READ(LUAUX,REC=INIREC + IREC -1) BUF
      READ(LUAUX,REC=INIREC + IREC   ) BUF
      DO   60 IWORD=1,LNPHYR
      Z(NREM + (IREC-1)*LNPHYR + IWORD) = BUF(IWORD)
   60 CONTINUE
   70 CONTINUE
      ENDIF
C
      RETURN
      END
