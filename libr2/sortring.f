      SUBROUTINE SORTRING(ICORE,MAXCOR,IUHF,IRREPX,IMODE,
     &                    LISTAAAA,LISTBBBB,LISTABAB,LISTBABA,
     &                    LISTBAAB,LISTABBA,SORTAAAA)
C
C   THIS ROUTINE DRIVES THE SORT OF THE H(IA,JB) INTERMEDIATES.
C
C   ALL H ARE RESORTED IN THE FOLLOWING WAY :
C
C    H(IA,JB) :   A,I ; B,J ---->     A,J ; B,I  [IMODE=1]
C    H(IA,JB) :   A,J ; B,I ---->     A,I ; B,J  [IMODE=2]
C
C NOTE THAT FOR RHF, ONLY LISTBAAB AND LIST ABAB ARE RESORTED [LIKE
C LISTS 56 AND 58, RESPECTIVELY].  THE AAAA TERMS ARE RESORTED ONLY
C IF SORTAAA IS SET TO .TRUE.
C
CEND
C
C  CODED SEPTEMBER/93  JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL CHANGE,SORTAAAA
      INTEGER DIRPRD,DISSYT
      CHARACTER*4 SPTYP,SPAAAA(2),SPABAB(2),SPABBA(2),SPABBA2(2)
      DIMENSION ICORE(MAXCOR)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON/MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA SPAAAA /'AAAA','BBBB'/,SPABAB /'ABAB','BABA'/
      DATA SPABBA /'AABB','BBAA'/,SPABBA2/'ABBA','BAAB'/ 
C    
C  LOOP OVER ISPINS TO PERFORM SORTS, WRITE ALL H4-AMPLITUDES BACK TO
C  THE SAME LIST
C
C
      DO 1000 ISPIN=1,IUHF+1
C
C  DO FIRST AAAA AND BBBB (UHF ONLY)
C
      IF(SORTAAAA)THEN
C
       LISTH=LISTAAAA-1+ISPIN
       NSIZE=IDSYMSZ(IRREPX,ISYTYP(1,LISTH),ISYTYP(2,LISTH))
       I001=1
       I002=I001+IINTFP*NSIZE
       I003=I002+IINTFP*NSIZE
       I004=I003+2*NOCCO(ISPIN)*NVRTO(ISPIN)*IINTFP
       IF(MAXCOR.LT.I004) THEN
        CALL INSMEM('SORTH',I004,MAXCOR)
       ENDIF
       CALL GETALL(ICORE(I001),NSIZE,IRREPX,LISTH)
       SPTYP=SPAAAA(ISPIN)
       CALL GSSTRNG(ICORE(I001),ICORE(I002),NSIZE,NSIZE,ICORE(I003),
     &             SPTYP,IRREPX)
       CALL PUTALL(ICORE(I002),NSIZE,IRREPX,LISTH)
C
      ENDIF
C
C  ABAB CASES 
C
      IF(ISPIN.EQ.1)THEN
       LISTH=LISTBAAB
      ELSE
       LISTH=LISTABBA
      ENDIF
      NSIZE=IDSYMSZ(IRREPX,ISYTYP(1,LISTH),ISYTYP(2,LISTH))
      I001=1
      I002=I001+IINTFP*NSIZE
      I003=I002+IINTFP*NSIZE
      I004=I003+2*IINTFP*NVRTO(ISPIN)*NOCCO(3-ISPIN)
      IF(MAXCOR.LT.I004) THEN
       CALL INSMEM('SORTH',I004,MAXCOR)
      ENDIF
      CALL GETALL(ICORE(I001),NSIZE,IRREPX,LISTH)
      SPTYP=SPABAB(ISPIN)
      CALL GSSTRNG(ICORE(I001),ICORE(I002),NSIZE,NSIZE,ICORE(I003),
     &            SPTYP,IRREPX)
      CALL PUTALL(ICORE(I002),NSIZE,IRREPX,LISTH)
C
C  DO ABBA CASES
C
      IF(ISPIN.EQ.1)THEN
       LISTH=LISTABAB
      ELSE
       LISTH=LISTBABA
      ENDIF
      NSIZE=IDSYMSZ(IRREPX,ISYTYP(1,LISTH),ISYTYP(2,LISTH))
      I001=1
      I002=I001+IINTFP*NSIZE
      I003=I002+IINTFP*NSIZE
      I004=I003+IINTFP*(NVRTO(ISPIN)+NVRTO(3-ISPIN))*
     &                 (NOCCO(ISPIN)+NOCCO(3-ISPIN))
      IF(MAXCOR.LT.I004) THEN
       CALL INSMEM('SORTH',I004,MAXCOR)
      ENDIF
      CALL GETALL(ICORE(I001),NSIZE,IRREPX,LISTH)
      IF(IMODE.EQ.1)THEN
       SPTYP=SPABBA(ISPIN)
       CALL GSSTRNG(ICORE(I001),ICORE(I002),NSIZE,NSIZE,ICORE(I003),
     &             SPTYP,IRREPX)
       LISTN=20+ISPIN
      ELSE
       SPTYP=SPABBA2(ISPIN)
       CALL GSSTRNG(ICORE(I001),ICORE(I002),NSIZE,NSIZE,ICORE(I003),
     &             SPTYP,IRREPX)
       LISTN=19-ISPIN
      ENDIF
   
      CHANGE=.TRUE.
      CALL NEWTYP2(IRREPX,LISTH,ISYTYP(1,LISTN),ISYTYP(2,LISTN),CHANGE)
      IOFF=I002
      DO 110 IRREPR=1,NIRREP
       IRREPL=DIRPRD(IRREPX,IRREPR)
       NUMSYT=IRPDPD(IRREPR,ISYTYP(2,LISTN))
       DISSYT=IRPDPD(IRREPL,ISYTYP(1,LISTN))
       CALL PUTLST(ICORE(IOFF),1,NUMSYT,1,IRREPR,LISTH)
       IOFF=IOFF+IINTFP*NUMSYT*DISSYT
110   CONTINUE 
1000  CONTINUE
      RETURN
      END
