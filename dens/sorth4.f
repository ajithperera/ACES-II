      SUBROUTINE SORTH4(ICORE,MAXCOR,IUHF)
C
C   THIS ROUTINE DRIVES THE SORT OF THE H4(IA,JB) INTERMEDIATES.
C   THE NEW ORDERING IS REQUIRED FOR THE EVALUATION OF THE
C   G5 AND G6 INTERMEDIATES
C
C   ALL H4 ARE RESORTED IN THE FOLLOWING WAY :
C
C    H(IA,JB) :   A,I ; B,J ---->     A,J ; B,I
C
CEND
C
C  CODED SEPTEMBER/90  JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL CHANGE
      INTEGER DIRPRD,DISSYT
      CHARACTER*4 SPTYP,SPAAAA(2),SPABAB(2),SPABBA(2)
      DIMENSION ICORE(MAXCOR)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/FLAGS/IFLAGS(100)
C
      DATA SPAAAA /'AAAA','BBBB'/,SPABAB /'ABAB','BABA'/
      DATA SPABBA /'AABB','BBAA'/
C
C CHECK WHICH LISTS HAVE BEEN USED FOR H INTERMEDIATES
C 
      IOFFLIST=0
      IF(IFLAGS(3).EQ.2) IOFFLIST=200
C    
C  LOOP OVER ISPINS TO PERFORM SORTS, WRITE ALL H4-AMPLITUDES BACK TO
C  THE SAME LIST
C
      DO 1000 ISPIN=1,IUHF+1
C
C  DO FIRST AAAA AND BBBB
C
      LISTH=53+ISPIN+IOFFLIST
      NSIZE=ISYMSZ(ISYTYP(1,LISTH),ISYTYP(2,LISTH))
      I001=1
      I002=I001+IINTFP*NSIZE
      I003=I002+IINTFP*NSIZE
      I004=I003+2*NOCCO(ISPIN)*NVRTO(ISPIN)*IINTFP
      IF(MAXCOR.LT.I004) THEN
       CALL INSMEM('SORTH4',I004,MAXCOR)
      ENDIF
      CALL GETALL(ICORE(I001),NSIZE,1,LISTH)
      SPTYP=SPAAAA(ISPIN)
      CALL SSTRNG(ICORE(I001),ICORE(I002),NSIZE,NSIZE,ICORE(I003),
     &            SPTYP)
      CALL PUTALL(ICORE(I002),NSIZE,1,LISTH)
C
C  ABAB CASES 
C
      LISTH=57+ISPIN+IOFFLIST
      NSIZE=ISYMSZ(ISYTYP(1,LISTH),ISYTYP(2,LISTH))
      I001=1
      I002=I001+IINTFP*NSIZE
      I003=I002+IINTFP*NSIZE
      I004=I003+2*IINTFP*NVRTO(ISPIN)*NOCCO(3-ISPIN)
      IF(MAXCOR.LT.I004) THEN
       CALL INSMEM('SORTH4',I004,MAXCOR)
      ENDIF
      CALL GETALL(ICORE(I001),NSIZE,1,LISTH)
      SPTYP=SPABAB(ISPIN)
      CALL SSTRNG(ICORE(I001),ICORE(I002),NSIZE,NSIZE,ICORE(I003),
     &            SPTYP)
      CALL PUTALL(ICORE(I002),NSIZE,1,LISTH)
     
C
C  DO ABBA CASES
C
      LISTH=55+ISPIN+IOFFLIST
      NSIZE=ISYMSZ(ISYTYP(1,LISTH),ISYTYP(2,LISTH))
      I001=1
      I002=I001+IINTFP*NSIZE
      I003=I002+IINTFP*NSIZE
      I004=I003+IINTFP*(NVRTO(ISPIN)+NVRTO(3-ISPIN))*
     &                 (NOCCO(ISPIN)+NOCCO(3-ISPIN))
      IF(MAXCOR.LT.I004) THEN
       CALL INSMEM('SORTH4',I004,MAXCOR)
      ENDIF
      CALL GETALL(ICORE(I001),NSIZE,1,LISTH)
      SPTYP=SPABBA(ISPIN)
      CALL SSTRNG(ICORE(I001),ICORE(I002),NSIZE,NSIZE,ICORE(I003),
     &            SPTYP)
      LISTN=20+ISPIN
      CHANGE=.TRUE.
      CALL NEWTYP(LISTH,ISYTYP(1,LISTN),ISYTYP(2,LISTN),CHANGE)
      IOFF=I002
      DO 110 IRREP=1,NIRREP
       NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTN))
       DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTN))
       CALL PUTLST(ICORE(IOFF),1,NUMSYT,1,IRREP,LISTH)
       IOFF=IOFF+IINTFP*NUMSYT*DISSYT
110   CONTINUE 
1000  CONTINUE
      RETURN
      END
