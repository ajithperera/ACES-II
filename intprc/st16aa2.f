      SUBROUTINE ST16AA2(W,BUF,IBUF,NUMIRW,ISYM,IPW,IPDIS,IPDSZ,
     &                   MXCOR,IRREPA,ILNBUF,ISPIN)
C
C SET UP PROCESSING OF ABCD INTEGRALS DIRECTLY FROM HF2 FILE
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION W,BUF
      LOGICAL GRAD
      DIMENSION W(MXCOR),NUMIRW(1),IPW(8),IPDIS(8),IPDSZ(8)
      DIMENSION ISYM(1),BUF(ILNBUF),IBUF(ILNBUF),IRREPA(*)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
C
      NSIZE=MXCOR/IINTFP
C
C CALCULATE TOTAL SIZE OF ABCD INTEGRALS
C
      ITOTSIZ=ISYMSZ(ISYTYP(1,230+ISPIN),ISYTYP(2,230+ISPIN))
C
C INITIALIZE BOTTOM ADDRESS FOR PASSES
C
      NPASS=0
      ILOWADR=1
      IRRWRT=1
      IDIS1=1
      DISLEFT=IRPDPD(IRRWRT,ISYTYP(2,230+ISPIN))-IDIS1+1
1     CONTINUE
C
C CALCULATE TOP ADDRESS FOR THIS PASS
C
      ITOPADR=MIN(ILOWADR+NSIZE-1,ITOTSIZ)
      NWORDS=ITOPADR-ILOWADR+1
C
      CALL RDABCDAA(W,NSIZE,ITOPADR,ILOWADR,BUF,IBUF,IRREPA,
     &              IPDSZ,IPDIS,ISYM,IPW,NUMIRW,ILNBUF,ISPIN)
      IOFFW=1
      NPASS=NPASS+1
C
C NOW SOME COMPLICATED LOGIC.  WRITE OUT WHATEVER WE HAVE IN CORE.
C
      NLIST=230+ISPIN
2     NDISLEFT=IRPDPD(IRRWRT,ISYTYP(2,NLIST))-IDIS1+1
      NWDSLEFT=NDISLEFT*IRPDPD(IRRWRT,ISYTYP(1,NLIST))
      IF(NWDSLEFT.LE.NWORDS)THEN
       CALL PUTLST(W(IOFFW),IDIS1,NDISLEFT,1,IRRWRT,NLIST)
       NWORDS=NWORDS-NDISLEFT*IRPDPD(IRRWRT,ISYTYP(1,NLIST))
       IOFFW=IOFFW+NDISLEFT*IRPDPD(IRRWRT,ISYTYP(1,NLIST))
       ILOWADR=ILOWADR+NDISLEFT*IRPDPD(IRRWRT,ISYTYP(1,NLIST))
       IRRWRT=IRRWRT+1
       IDIS1=1
       IDONE=0
       IF(IRRWRT.LE.NIRREP)GOTO 2
      ELSE
       NDISWRIT=NWORDS/IRPDPD(IRRWRT,ISYTYP(1,NLIST))
       CALL PUTLST(W(IOFFW),IDIS1,NDISWRIT,1,IRRWRT,NLIST)
       IDIS1=IDIS1+NDISWRIT
       ILOWADR=ILOWADR+NDISWRIT*IRPDPD(IRRWRT,ISYTYP(1,NLIST))
       IDONE=1
      ENDIF
C
      IF(ILOWADR.LE.ITOTSIZ)GOTO 1
C
      CLOSE(UNIT=25,STATUS='DELETE')
C
      WRITE(6,1000)NPASS
1000  FORMAT(T3,'@ST16AA2-I, ABCD integral processing required ',
     &       I5,' passes.')
C
      RETURN
      END