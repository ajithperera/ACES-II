      SUBROUTINE ULIST(Z,BUF,LENGTH,IDIS,IOFF,
     &                 ICACHE,IRREPL,IRREPR,LIST,IUP)
C
C THIS SUBROUTINE DUMPS THE ARRAY Z ON LIST `LIST'. THE PARAMETERS
C ARE:
C
C   Z....... ARRAY TO DUMP TO DISK
C   BUF..... SCRATCH ARRAY OF LENGTH DISSIZ
C   LENGTH.. LENGTH OF ARRAY Z
C   IDIS.... THE NUMBER OF DISTRIBUTION TO WHICH Z IS DUMPED
C   IOFF.... OFFSET IN THE DISTRIBUTION TO WHICH Z IS DUMPED
C   ICACHE.. CACHE NUMBER
C   IRREPL.. IRREP OF THE LEFT SIDE
C   IRREPR.. IRREP OF THE RIGHT SIDE
C   LIST.... LIST NUMBER
C   IUP..... 0 = INITIALIZE, 1 = COPY, 2 = ADD Z TO LIST
C
CEND
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DISSIZ
      DIMENSION Z(LENGTH),BUF(1000)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
C 
      DATA ONE/1.D0/
C
      IF(IUP.EQ.0) THEN
       DISSIZ=IRPDPD(IRREPL,ISYTYP(1,LIST))
       CALL IZERO(BUF,IINTFP*DISSIZ)
c YAU : old
c      CALL ICOPY(LENGTH*IINTFP,Z,1,BUF(IOFF+1),1)
c YAU : new
       CALL DCOPY(LENGTH,Z,1,BUF(IOFF+1),1)
c YAU : end
      ELSE IF(IUP.EQ.1) THEN
       CALL GETLST(BUF,IDIS,1,ICACHE,IRREPR,LIST)
c YAU : old
c      CALL ICOPY(LENGTH*IINTFP,Z,1,BUF(IOFF+1),1)
c YAU : new
       CALL DCOPY(LENGTH,Z,1,BUF(IOFF+1),1)
c YAU : end
      ELSE IF(IUP.EQ.2) THEN
       CALL GETLST(BUF,IDIS,1,ICACHE,IRREPR,LIST)
       CALL SAXPY(LENGTH,ONE,Z,1,BUF(IOFF+1),1) 
      ENDIF 
       DISSIZ=IRPDPD(IRREPL,ISYTYP(1,LIST))
      CALL PUTLST(BUF,IDIS,1,ICACHE,IRREPR,LIST)
      RETURN
      END
