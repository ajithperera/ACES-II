
      SUBROUTINE SG1AAAA(BUCK,BUF,BUF1,IBUCK,NINBCK,ICHAIN,NBUCK,
     &                   NDBCK,NBKINT,NREC,IRECL,ILNBUF,
     &                   BUF2,BUF3,ISPIN,IRREPDO)
C
C OUT OF CORE SORT OF THE MO GAMMAS FOR A SPECIFIC IRREP.
C
CEND
      IMPLICIT INTEGER (A-Z)
      PARAMETER (LUSRT=20)
      CHARACTER*80 FNAME
      DOUBLE PRECISION BUF(ILNBUF),X,BUF1(*)
      DOUBLE PRECISION BUF2(*),BUF3(*),ONE,TWO
      DOUBLE PRECISION BUCK(NBKINT,NBUCK)
      DOUBLE PRECISION FABIJ,FAIBJ,FABCD,FIJKL,FABCI,FIJKA
      INTEGER IBUCK(NBKINT,NBUCK),NINBCK(NBUCK),ICHAIN(NBUCK)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      LOGICAL MBPT4
      LOGICAL ROHF
      LOGICAL TDA,EOM
      LOGICAL DABCD,GABCD
C   
      COMMON /EXCITE/ TDA,EOM
      COMMON /STATSYM/IRREPX
      COMMON /METH/ MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /MOPOPS/ MOPOP(8)
      COMMON /FACTORS/ FABIJ,FABCD,FIJKL,FAIBJ,FABCI,FIJKA
      COMMON /REF/ ROHF
      COMMON /IABCD/ IOFFABCD
      COMMON /ABCD/ GABCD,DABCD
C
      DATA ONE,TWO /1.D0,2.D0/
C
      INDXF(I,J,N)=I+(J-1)*N
      INDXT(I,J)  =I+(J*(J-1))/2
      NNP1O2(I)   =(I*(I+1))/2
      IPACK(I,J)  =OR(I,ISHFT(J,IBITWD))
C
      NREC=1
      CALL ZERO(BUCK,NBKINT*NBUCK)
      CALL IZERO(IBUCK,NBKINT*NBUCK)
      CALL IZERO(ICHAIN,NBUCK)
      CALL IZERO(NINBCK,NBUCK)
      CALL GFNAME('SRTFIL  ',FNAME,ILENGTH)
      OPEN(UNIT=20,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &     ACCESS='DIRECT',RECL=IRECL)
C
      IF(TDA)GOTO 3000
C
C FIRST PROCESS THE G(AB,IJ) LIST
C EACH ELEMENT IS PUT INTO FOUR PLACES
C
C   G(AI,BJ), -G(BI,AJ), -G(AJ,BI), G(BJ,AI)
C
C   FACTOR 0.5 HERE, IN ORDER TO BE CONSISTENT WITH SECOND ORDER
C
      CALL SYMPOS('T',POP(1,ISPIN),POP(1,ISPIN),IRREPDO,1,ILOGREC)
      CALL SYMPOS('T',VRT(1,ISPIN),VRT(1,ISPIN),IRREPDO,1,ITHRU1)
      LISTABIJ=113+ISPIN
       NUMOCC=POP(IRREPDO,ISPIN)
       NUMVRT=VRT(IRREPDO,ISPIN)
       NUMAO=MOPOP(IRREPDO)
       MAXL=NNP1O2(NUMAO)
       DO 101 INDXJ=2,NUMOCC
        INDJ=INDXJ
        DO 102 INDXI=1,INDXJ-1
         ILOGREC=ILOGREC+1
         INDI=INDXI
         CALL GETLST(BUF,ILOGREC,1,1,1,LISTABIJ)
         ITHRU=ITHRU1
         DO 103 INDXB=2,NUMVRT 
          INDB=INDXB+NUMOCC
*VOCL LOOP,NOVREC
CDIR$ IVDEP
          DO 104 INDXA=1,INDXB-1
           ITHRU=ITHRU+1
           X=BUF(ITHRU)*FABIJ
           INDA=INDXA+NUMOCC
           GINDXL=INDXT(INDI,INDA)
           GINDXR=INDXT(INDJ,INDB)
           IBKET=1+(GINDXR-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXL,GINDXR)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
           IBKET=1+(GINDXL-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXR,GINDXL)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
           GINDXL=INDXT(INDJ,INDA)
           GINDXR=INDXT(INDI,INDB)
           IBKET=1+(GINDXR-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=-X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXL,GINDXR)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
           IBKET=1+(GINDXL-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=-X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXR,GINDXL)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
104       CONTINUE
103      CONTINUE
102     CONTINUE
101    CONTINUE
C
C  FOR SECOND ORDER NOTHING ELSE TO DO
C
      IF(MBPT2)GOTO 2000
C
3000  CONTINUE
C
C NOW PROCESS THE G(AI,BJ) LIST.  VALUES ARE PUT
C  IN TWO PLACES IN SQUARE CANONICAL LIST [G(AB,IJ) AND G(IJ,AB)]
C  FACTOR OF 0.5 AGAIN IN ORDER TO BE CONSISTENT
C
      CALL SYMPOS('F',VRT(1,ISPIN),POP(1,ISPIN),IRREPDO,1,ILOGREC)
      CALL SYMPOS('F',VRT(1,ISPIN),POP(1,ISPIN),IRREPDO,1,IOFFIR)
      LISTAIBJ=122+ISPIN
       NUMOCC=POP(IRREPDO,ISPIN)
       NUMVRT=VRT(IRREPDO,ISPIN)
       NUMAO=MOPOP(IRREPDO)
       MAXL=NNP1O2(NUMAO)
       DO 201 INDXJ=1,NUMOCC
        INDJ=INDXJ
        DO 202 INDXB=1,NUMVRT
         ILOGREC=ILOGREC+1
         INDB=INDXB+NUMOCC
         CALL GETLST(BUF,ILOGREC,1,1,1,LISTAIBJ)
         DO 203 INDXI=1,INDXJ
          INDI=INDXI
          GINDXR=INDXT(INDI,INDJ)
          DO 204 INDXA=1,INDXB
           IPOS=INDXF(INDXA,INDXI,NUMVRT)+IOFFIR
           X=BUF(IPOS)*FAIBJ
           INDA=INDXA+NUMOCC
           GINDXL=INDXT(INDA,INDB)
           IBKET=1+(GINDXR-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXL,GINDXR)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
           IBKET=1+(GINDXL-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXR,GINDXL)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
204       CONTINUE
          DO 205 INDXA=INDXB,NUMVRT
           IPOS=INDXF(INDXA,INDXI,NUMVRT)+IOFFIR
           X=BUF(IPOS)*FAIBJ
           INDA=INDXA+NUMOCC
           GINDXL=INDXT(INDB,INDA)
           IBKET=1+(GINDXR-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXL,GINDXR)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
           IBKET=1+(GINDXL-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXR,GINDXL)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
205       CONTINUE
203      CONTINUE
         DO 206 INDXI=1,NUMOCC
          INDI=INDXI
          GINDXL=INDXT(INDI,INDB)
          DO 207 INDXA=1,NUMVRT
           IPOS=INDXF(INDXA,INDXI,NUMVRT)+IOFFIR
           X=BUF(IPOS)*FAIBJ
           INDA=INDXA+NUMOCC
           GINDXR=INDXT(INDJ,INDA)
           IBKET=1+(GINDXR-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=-X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXL,GINDXR)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
207       CONTINUE
206      CONTINUE
202     CONTINUE
201    CONTINUE
C
       IF(TDA)GOTO 2000
C
C NOW PROCESS THE G(IJ,KL) LIST
C EACH ELEMENT CONTRIBUTES IN THE FOLLOWING WAY : 
C
C  G(IK,JL), G(JL,IK), - G(JK,IL), - G(IL,JK)
C
C  FACTOR : 0.5
C
C A TOTAL OF FOUR CONTRIBUTIONS
C
      CALL SYMPOS('T',POP(1,ISPIN),POP(1,ISPIN),IRREPDO,1,ILOGREC)
      CALL SYMPOS('T',POP(1,ISPIN),POP(1,ISPIN),IRREPDO,1,IOFFIR)
      LISTIJKL=110+ISPIN
       NUMOCC=POP(IRREPDO,ISPIN)
       NUMAO=MOPOP(IRREPDO)
       MAXL=NNP1O2(NUMAO)
       DO 301 INDXL=2,NUMOCC
        INDL=INDXL
        DO 302 INDXK=1,INDXL-1
         ILOGREC=ILOGREC+1
         INDK=INDXK
         CALL GETLST(BUF,ILOGREC,1,1,1,LISTIJKL)
         CALL EXPND3(BUF(IOFFIR+1),BUF1,NUMOCC)
         DO 303 INDXJ=1,INDXL
          INDJ=INDXJ
          GINDXR=INDXT(INDJ,INDL)
          DO 304 INDXI=1,INDXK
           IPOS=INDXF(INDXI,INDXJ,NUMOCC)
           X=BUF1(IPOS)*FIJKL
           INDI=INDXI
           GINDXL=INDXT(INDI,INDK)
           IBKET=1+(GINDXR-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXL,GINDXR)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
304       CONTINUE
          DO 305 INDXI=INDXK,NUMOCC
           IPOS=INDXF(INDXI,INDXJ,NUMOCC)
           X=BUF1(IPOS)*FIJKL
           INDI=INDXI
           GINDXL=INDXT(INDK,INDI)
           IBKET=1+(GINDXR-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXL,GINDXR)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
305       CONTINUE
303      CONTINUE
         DO 306 INDXJ=1,INDXK
          INDJ=INDXJ
          GINDXR=INDXT(INDXJ,INDXK)
          DO 307 INDXI=1,INDXL    
           IPOS=INDXF(INDXI,INDXJ,NUMOCC)
           X=BUF1(IPOS)*FIJKL
           INDI=INDXI
           GINDXL=INDXT(INDI,INDL)
           IBKET=1+(GINDXR-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=-X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXL,GINDXR)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
307       CONTINUE
          DO 308 INDXI=INDXL,NUMOCC
           IPOS=INDXF(INDXI,INDXJ,NUMOCC)
           X=BUF1(IPOS)*FIJKL
           INDI=INDXI
           GINDXL=INDXT(INDL,INDI)   	
           IBKET=1+(GINDXR-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=-X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXL,GINDXR)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
308       CONTINUE
306      CONTINUE
302     CONTINUE
301    CONTINUE
C
C NOW PROCESS THE G(AB,CD) LIST
C
C AGAIN THERE ARE FOUR CONTRIBUTIONS
C
C  G(AC,BD), G(BD,AC), - G(AD,BC), -G(BC,AD)
C
C  FACTOR : 0.5 
C
      NPASS=1
      IF(EOM.AND.DABCD) NPASS=2
C
C
      DO 410 IPASS=1,NPASS
C
       IF(DABCD) THEN
C
        IF(IPASS.EQ.1) THEN
C
         NUMDSZ=IRPDPD(1,ISYTYP(2,143+ISPIN))
         DISDSZ=IRPDPD(1,ISYTYP(1,143+ISPIN))
C
         I0TA=1
         CALL GETLST(BUF2,1,NUMDSZ,1,1,153+ISPIN)
         CALL GETLST(BUF3,1,NUMDSZ,1,1,143+ISPIN)
         IF(CCSD) THEN
           CALL GETLST(BUF(I0TA),1,1,1,ISPIN,157)
           CALL FTAU(BUF2,BUF(I0TA),BUF(I0TA),DISDSZ,
     &               NUMDSZ,POP(1,ISPIN),POP(1,ISPIN),
     &               VRT(1,ISPIN),VRT(1,ISPIN),1,ISPIN,ONE)
         ENDIF
         MBPT4=M4DQ.OR.M4SDQ.OR.M4SDTQ
         IF(MBPT4) THEN
          CALL SSCAL(NUMDSZ*DISDSZ,TWO,BUF2,1)
          CALL SAXPY(NUMDSZ*DISDSZ,ONE,BUF1,1,BUf2,1)
         ENDIF
C
        ELSE
C
         NUMDSZ=IRPDPD(IRREPX,ISYTYP(2,143+ISPIN))
         DISDSZ=IRPDPD(1,ISYTYP(1,143+ISPIN))
C
         I0TA=1
         I0RA=I0TA+NT(ISPIN)
C
         CALL GETLST(BUF2,1,NUMDSZ,1,IRREPX,460+ISPIN)
         CALL GETLST(BUF3,1,NUMDSZ,1,IRREPX,443+ISPIN)
C
         CALL GETLST(BUF(I0TA),1,1,1,ISPIN,157)
         CALL GETLST(BUF(I0RA),1,1,1,2+ISPIN,490)
C
         CALL DTAU(1,IRREPX,1,IRREPX,BUF2,
     &             BUF(I0TA),BUF(I0TA),BUF(I0RA),BUF(I0RA),
     &             ISPIN,ONE)
C
        ENDIF
C
       ENDIF
C
      CALL SYMPOS('T',VRT(1,ISPIN),VRT(1,ISPIN),IRREPDO,1,ILOGREC)
      CALL SYMPOS('T',VRT(1,ISPIN),VRT(1,ISPIN),IRREPDO,1,IOFFIR)
      LISTABCD=IOFFABCD+130+ISPIN
       NUMVRT=VRT(IRREPDO,ISPIN)
       OCC=POP(IRREPDO,ISPIN)
       NUMAO=MOPOP(IRREPDO)
       MAXL=NNP1O2(NUMAO)
       DO 401 INDXD=2,NUMVRT
        INDD=INDXD+OCC
        DO 402 INDXC=1,INDXD-1
         ILOGREC=ILOGREC+1
         INDC=INDXC+OCC
C
         IF(.NOT.DABCD) THEN
C
          CALL GETLST(BUF,ILOGREC,1,1,1,LISTABCD)
C
         ELSE
C
          CALL DIRG2(BUF,BUF2,BUF3,BUF2,BUF3,NUMDSZ,DISDSZ,
     &               ILOGREC,1,1)
C
         ENDIF
         CALL EXPND3(BUF(IOFFIR+1),BUF1,NUMVRT)
         DO 403 INDXB=1,INDXD
          INDB=INDXB+OCC
          GINDXR=INDXT(INDB,INDD)
          DO 404 INDXA=1,INDXC
           IPOS=INDXF(INDXA,INDXB,NUMVRT)
           X=BUF1(IPOS)*FABCD
           INDA=INDXA+OCC
           GINDXL=INDXT(INDA,INDC)
           IBKET=1+(GINDXR-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXL,GINDXR)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
404       CONTINUE
          DO 405 INDXA=INDXC,NUMVRT
           IPOS=INDXF(INDXA,INDXB,NUMVRT)
           X=BUF1(IPOS)*FABCD
           INDA=INDXA+OCC
           GINDXL=INDXT(INDC,INDA)
           IBKET=1+(GINDXR-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXL,GINDXR)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
405       CONTINUE
403      CONTINUE
         DO 406 INDXB=1,INDXC
          INDB=INDXB+OCC
          GINDXR=INDXT(INDB,INDC)
          DO 407 INDXA=1,INDXD
           IPOS=INDXF(INDXA,INDXB,NUMVRT)
           X=BUF1(IPOS)*FABCD
           INDA=INDXA+OCC
           GINDXL=INDXT(INDA,INDD)
           IBKET=1+(GINDXR-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=-X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXL,GINDXR)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
407       CONTINUE
          DO 408 INDXA=INDXD,NUMVRT
           IPOS=INDXF(INDXA,INDXB,NUMVRT)
           X=BUF1(IPOS)*FABCD
           INDA=INDXA+OCC
           GINDXL=INDXT(INDD,INDA)
           IBKET=1+(GINDXR-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=-X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXL,GINDXR)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
408       CONTINUE
406      CONTINUE
402     CONTINUE
401    CONTINUE
410   CONTINUE
C
C NOW PROCESS THE G(ABCI) LIST
C
C  FOR MBPT(3) AND CCD SKIP THIS PART AND GO TO 2000
C  HOWEVER, FOR ROHF-MBPT(3) THIS PART HAS TO CONSIDERED
C
C FOUR CONTRIBUTIONS : G(BI,AC) AND -G(AI,BC) AND G(AC,BI) AND -G(BC,AI)
C
C  FACTOR 0.5 D0
C
      IF((MBPT3.AND.(.NOT.ROHF)).OR.CCD)GOTO 2000
      CALL SYMPOS('F',VRT(1,ISPIN),POP(1,ISPIN),IRREPDO,1,ILOGREC)
      CALL SYMPOS('T',VRT(1,ISPIN),VRT(1,ISPIN),IRREPDO,1,IOFFIR)
      LISTABCI=126+ISPIN
       NUMOCC=POP(IRREPDO,ISPIN)
       NUMVRT=VRT(IRREPDO,ISPIN)
       NUMAO=MOPOP(IRREPDO)
       MAXL=NNP1O2(NUMAO)
       DO 501 INDXI=1,NUMOCC
        INDI=INDXI
        DO 502 INDXC=1,NUMVRT
         ILOGREC=ILOGREC+1
         INDC=INDXC+NUMOCC
         CALL GETLST(BUF,ILOGREC,1,1,1,LISTABCI)
         CALL EXPND3(BUF(IOFFIR+1),BUF1,NUMVRT)
         DO 503 INDXB=1,NUMVRT
          INDB=INDXB+NUMOCC  
          GINDXR=INDXT(INDI,INDB)
          DO 504 INDXA=1,INDXC
           IPOS=INDXF(INDXA,INDXB,NUMVRT)
           X=BUF1(IPOS)*FABCI
           INDA=INDXA+NUMOCC
           GINDXL=INDXT(INDA,INDC)
           IBKET=1+(GINDXR-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXL,GINDXR)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
           IBKET=1+(GINDXL-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXR,GINDXL)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
504       CONTINUE
          DO 505 INDXA=INDXC,NUMVRT
           IPOS=INDXF(INDXA,INDXB,NUMVRT)
           X=BUF1(IPOS)*FABCI
           INDA=INDXA+NUMOCC
           GINDXL=INDXT(INDC,INDA)
           IBKET=1+(GINDXR-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXL,GINDXR)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
           IBKET=1+(GINDXL-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXR,GINDXL)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
505       CONTINUE
503      CONTINUE
502     CONTINUE
501    CONTINUE
C
C NOW PROCESS THE G(IJKA) LIST
C
C AGAIN FOUR CONTRIBUTIONE G(IK,AJ) AND -G(JK,AI) AND 
C G(AJ,IK) AND -G(AI,JK)
C
      CALL SYMPOS('F',POP(1,ISPIN),VRT(1,ISPIN),IRREPDO,1,ILOGREC)
      CALL SYMPOS('T',POP(1,ISPIN),POP(1,ISPIN),IRREPDO,1,IOFFIR)
      LISTIJKA=106+ISPIN
       NUMVRT=VRT(IRREPDO,ISPIN)
       NUMOCC=POP(IRREPDO,ISPIN)
       NUMAO=MOPOP(IRREPDO)
       MAXL=NNP1O2(NUMAO)
       DO 601 INDXA=1,NUMVRT
        DO 602 INDXK=1,NUMOCC
         ILOGREC=ILOGREC+1
         INDA=INDXA+NUMOCC
         INDK=INDXK
         CALL GETLST(BUF,ILOGREC,1,1,1,LISTIJKA)
         CALL EXPND3(BUF(IOFFIR+1),BUF1,NUMOCC)
         DO 603 INDXJ=1,NUMOCC
          INDJ=INDXJ
          GINDXR=INDXT(INDJ,INDA)
          DO 604 INDXI=1,INDXK
           IPOS=INDXF(INDXI,INDXJ,NUMOCC)
           X=BUF1(IPOS)*FIJKA
           INDI=INDXI
           GINDXL=INDXT(INDI,INDK)
           IBKET=1+(GINDXR-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXL,GINDXR)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
           IBKET=1+(GINDXL-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXR,GINDXL)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
604       CONTINUE
          DO 605 INDXI=INDXK,NUMOCC
           IPOS=INDXF(INDXI,INDXJ,NUMOCC)
           X=BUF1(IPOS)*FIJKA
           INDI=INDXI
           GINDXL=INDXT(INDK,INDI)
           IBKET=1+(GINDXR-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXL,GINDXR)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
           IBKET=1+(GINDXL-1)/NDBCK
           NINBCK(IBKET)=NINBCK(IBKET)+1
           BUCK(NINBCK(IBKET),IBKET)=X
           IBUCK(NINBCK(IBKET),IBKET)=IPACK(GINDXR,GINDXL)
           IF(NINBCK(IBKET).EQ.NBKINT)THEN
            CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &                 NINBCK(IBKET),NBKINT,NREC)
           ENDIF
605       CONTINUE
603      CONTINUE
602     CONTINUE
601    CONTINUE
C
C FLUSH REMAINING BUFFERS
C
2000  CONTINUE
C
      DO 5000 IBKET=1,NBUCK
       IF(NINBCK(IBKET).NE.0)THEN
        CALL PLUNK(LUSRT,BUCK(1,IBKET),IBUCK(1,IBKET),ICHAIN(IBKET),
     &             NINBCK(IBKET),NBKINT,NREC)
       ENDIF
5000  CONTINUE
C
      RETURN
      END
