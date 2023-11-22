C
C PICKS UP INTEGRALS OFF OF SORT FILE AND PUTS THEM OUT ON MOINTS.
C  FOR FULL SORTS OF <IJ|KL> TYPE QUANTITIES (ALL I,J FOR A GIVEN K,L).
C
      SUBROUTINE CABCD1(W,W1,BUF,IBUF,ICHAIN,DISSIZ,NDBCK,
     &                  NBKINT,NBUCK,IRECL,NUMIRW,ISYM,
     &                  IPW,IPDIS,NWDIS,ISPIN,NMO,IRREPA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSIZ,BUFSIZ,DIRPRD,VRT,POP
      INTEGER AND
      CHARACTER*80 FNAME
      DIMENSION W(DISSIZ*NDBCK),W1(1),BUF(NBKINT),IBUF(NBKINT)
      DIMENSION ICHAIN(NBUCK),ISYM(1),NUMIRW(2*NIRREP),
     &          IPW(8),IPDIS(8) ,IRREPA(NMO,2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON/SYM/POP(8,2),VRT(8,2),NJUNK(6)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      IUPKI(INT)=AND(INT,IALONE)
      IUPKJ(INT)=AND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=AND(ISHFT(INT,-2*IBITWD),IALONE)
      IUPKL(INT)=AND(ISHFT(INT,-3*IBITWD),IALONE)
      INDX1(I,J,N)=I+(J-1)*N
      INDX2(I,J)=I+((J-1)*(J-2))/2
      JDIS=0
      IRREPW=1
      NOCC=NOCCO(ISPIN)
      NVRT=NVRTO(ISPIN)
C
      CALL GFNAME('SRTFIL  ',FNAME,ILENGTH)
      OPEN(UNIT=20,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &ACCESS='DIRECT',RECL=IRECL)
C
      DO 10 IDIS=1,NBUCK
       ISTART=(IDIS-1)*NDBCK+1
       DO 1 IRREP=1,NIRREP
       IRREPST=IRREP
       IF(IPDIS(IRREP).GE.ISTART) THEN
        GO TO 2
       ENDIF
1      CONTINUE
       IRREPST=IRREPST+1
2      CONTINUE
       IRREPST=IRREPST-1
       IOFFST1=IPW(IRREPST)
       NNST=NUMIRW(IRREPST+NIRREP)
       IOFFST=IOFFST1+NNST*(ISTART-1-IPDIS(IRREPST))
       CALL ZERO(W,DISSIZ*NDBCK)
5      READ(20,REC=ICHAIN(IDIS))BUF,IBUF,NUT,ICHAN
       ICHAIN(IDIS)=ICHAN
       DO 30 INT=1,NUT
        I=IUPKI(IBUF(INT))
        J=IUPKJ(IBUF(INT))
        K=IUPKK(IBUF(INT))
        L=IUPKL(IBUF(INT))
        ISQ1=INDX1(I,J,NVRT)
        ISQ2=INDX2(K,L)
        IRREP12=DIRPRD(IRREPA(I+NOCC,ISPIN),IRREPA(J+NOCC,ISPIN))
        IOFF=IPW(IRREP12)
        NN=NUMIRW(IRREP12+NIRREP)
        IADR=IOFF+ISYM(ISQ1+NWDIS)
     &           +NN*(ISYM(ISQ2)-1-IPDIS(IRREP12))-IOFFST
        W(IADR)=BUF(INT)
30     CONTINUE
       IF(ICHAIN(IDIS).NE.0)GOTO 5
C
C  IF WE REACH THIS POINT, THE ENTIRE BUCKET IS IN CORE. WRITE
C  THEM TO A DIRECT ACCESS FILE 
C
       IOFFW=1
       DO 40 ID=1,NDBCK
39      JDIS=JDIS+1
        IF(JDIS.GT.NUMIRW(IRREPW)) THEN
        IRREPW=IRREPW+1
        IF(IRREPW.GT.NIRREP) GO TO 35
         JDIS=0
         GO TO 39
        ENDIF
        CALL ASSYM(IRREPW,VRT(1,ISPIN),1,1,W1,W(IOFFW))
        CALL PUTLST(W1,JDIS,1,2,IRREPW,230+ISPIN)
        IOFFW=IOFFW+NUMIRW(IRREPW+NIRREP)
40     CONTINUE
35     CONTINUE
10    CONTINUE
      NUMREC=0
      DO 50 I=1,NIRREP
       NUMREC=NUMREC+NUMIRW(IRREP+NIRREP)
50    CONTINUE
      WRITE(LUOUT,300) NUMREC
300   FORMAT(T3,'@CABCD0-I, Wrote ',I6,' logical records of ',
     &          'symmetry blocked integrals.')
      CLOSE(UNIT=20,STATUS='DELETE')
      RETURN
      END