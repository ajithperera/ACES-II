      SUBROUTINE W5OAB2(ICORE,MAXCOR,IUHF,T1,TARGET,
     &                  IRREP,LISTW,ISPINT,WDSZ)
C
C THIS ROUTINE PERFORMS THE CONTRACTION
C                                        
C             Z(Ef,Am) = W(Ef,Ag) * T(g,m)
C
C USING AN ALGORITHM WHICH DOES NOT ASSUME THAT W IS HELD IN
C CORE, BUT RATHER THAT N**3 ELEMENTS CAN FIT IN CORE, WHERE
C N IS THE NUMBER OF BASIS FUNCTIONS.  THIS ROUTINE PERFORMS 
C THE CONTRACTION FOR ONLY ONE IRREP SO IT MUST BE CALLED FROM
C INSIDE THE MAIN IRREP LOOP OF THE CORRESPONDING CONTRACTION IN
C W5AB1.
C
C THE TARGET ARRAY IS PASSED IN AS Z(Ef,Am).
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION TARGET,T1,ONE
      DIMENSION ICORE(MAXCOR),T1(1),TARGET(1),IOFFT1(8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      DATA ONE /1.0/
C
C COMPUTE T1 OFFSETS
C
      IOFF=1
      DO 5 IRREPT=1,NIRREP
       IOFFT1(IRREPT)=IOFF
       IOFF=IOFF+POP(IRREPT,2)*VRT(IRREPT,2)
5     CONTINUE
C
C THE W INTERMEDIATES ARE ON DISK IN THE MATRIX W(Ef,Ag)
C
C
C FIRST LOOP OVER IRREDUCIBLE REPRESENTATIONS OF g
C
      NFIRST0=1
      IOFFTAR=1
      DO 10 IRREPG=1,NIRREP
       IRREPA=DIRPRD(IRREPG,IRREP)
       IRREPM=IRREPG
C
       NUMA=VRT(IRREPA,3-ISPINT)
       NUMG=VRT(IRREPG,ISPINT)
       NUMM=POP(IRREPM,ISPINT)
       NUMEFA=WDSZ*NUMA
C
C DETERMINE HOW MANY W(EfA) VECTORS CAN FIT IN CORE.
C
       IF(NUMEFA.EQ.0)GOTO 10
       NINCOR=MAXCOR/(NUMEFA*IINTFP)
       IF(NINCOR.EQ.0)THEN
        WRITE(6,1000)
1000    FORMAT(T3,'@W5OAB2-F, Sorry.  Not enough memory to perform',
     &            ' W(Ef,Ag)*T(g,m) contraction.')
        CALL INSMEM('W5OAB2',NUMEFA*IINTFP,MAXCOR)
       ENDIF
C
C GOOD.  WE CAN ACTUALLY DO THIS CALCULATION.  NOW LOAD UP
C SOME W(EfA,g) AND GO TO WORK.
C
       NLEFT=NUMG
       NFIRST=NFIRST0
       IOFFT=IOFFT1(IRREPM)
1      NBATCH=MIN(NLEFT,NINCOR)
       NREAD =NBATCH*NUMA
       CALL GETLST(ICORE,NFIRST,NREAD,1,IRREP,LISTW)

CSSS       call checksum("WABCD",ICORE,NREAD,S)
       IF(MIN(NUMEFA,NUMM,NBATCH).NE.0)THEN
        CALL XGEMM('N','N',NUMEFA,NUMM,NBATCH,ONE,ICORE,NUMEFA,
     &             T1(IOFFT),NUMG,ONE,TARGET(IOFFTAR),NUMEFA)
       ENDIF
       NFIRST=NFIRST+NREAD
       NLEFT=NLEFT-NBATCH
       IOFFT=IOFFT+NBATCH
       IF(NLEFT.NE.0)GOTO 1
       IOFFTAR=IOFFTAR+NUMEFA*NUMM
       NFIRST0=NFIRST0+NUMA*NUMG
10    CONTINUE
      RETURN
      END
