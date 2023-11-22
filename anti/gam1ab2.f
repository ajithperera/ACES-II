      SUBROUTINE GAM1AB2(GAMMA,BUF,BUF1,BUF2,SCR1,SCR2,
     &                  LENGAM,IRREPDO,LIST1,LIST2,IUHF)
C
C THIS ROUTINE FORMS THE ABAB SYMMETRY TYPE AO GAMMA LISTS
C  FOR SPIN CASE AB AND WRITES THEM TO THE BACKTRANSFORMATION
C  MOINTS FILE FOR A PARTICULAR IRREP. THIS IS THE INCORE VERSION!
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION GAMMA(LENGAM),BUF(*),X,ONE,TWO
      DOUBLE PRECISION BUF1(*),BUF2(*),SCR1(*),SCR2(*)
      DOUBLE PRECISION FABIJ,FAIBJ,FABCD,FIJKL,FABCI,FIJKA
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      LOGICAL MBPT4
      LOGICAL ROHF,MOEQAO
      LOGICAL GABCD,DABCD
      LOGICAL TDA,EOM
C
      COMMON /METH/ MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /MOPOPS/ MOPOP(8)
      COMMON /FACTORS/ FABIJ,FABCD,FIJKL,FAIBJ,FABCI,FIJKA
      COMMON /REF/ ROHF
      COMMON /SIZES/ MOEQAO
      COMMON /IABCD/ IOFFABCD
      COMMON /ABCD/ GABCD,DABCD
      COMMON /EXCITE/ TDA,EOM
      COMMON /STATSYM/IRREPX
C
      DATA ONE,TWO /1.D0,2.D0/
C
      INDXF(I,J,N)=I+(J-1)*N
      INDXT(I,J)  =I+(J*(J-1))/2
      NNP1O2(I)   =(I*(I+1))/2
C
      IF(IUHF.EQ.0) THEN
C
C NOW PROCESS THE G(Ab,Ci) LIST
C
      CALL SYMPOS('F',VRT(1,1),POP(1,2),IRREPDO,1,ILOGREC)
      CALL SYMPOS('F',VRT(1,1),VRT(1,2),IRREPDO,1,IOFFIR)
       NUMI=POP(IRREPDO,2)
       NUMC=VRT(IRREPDO,1)
       NUMB=VRT(IRREPDO,2)
       NUMA=VRT(IRREPDO,1)
       NUMAO=MOPOP(IRREPDO)
       OCCA=POP(IRREPDO,1)
       OCCB=POP(IRREPDO,2)
       MAXL=NNP1O2(NUMAO)
       DO 501 INDXI=1,NUMI
        DO 502 INDXC=1,NUMC
         ILOGREC=ILOGREC+1
         INDI=INDXI
         INDC=INDXC+OCCA
         CALL GETLST(BUF,ILOGREC,1,1,1,130)
         DO 503 INDXB=1,NUMB
          INDB=INDXB+OCCB
          GINDXR=INDXT(INDI,INDB)
          DO 504 INDXA=1,INDXC
           IPOS=INDXF(INDXA,INDXB,NUMA)+IOFFIR
           X=BUF(IPOS)*FABCI
           INDA=INDXA+OCCA
           GINDXL=INDXT(INDA,INDC)

           gamma(indxf(gindxl,gindxr,maxl))=
     &        gamma(indxf(gindxl,gindxr,maxl))+x

           gamma(indxf(gindxr,gindxl,maxl))=
     &        gamma(indxf(gindxr,gindxl,maxl))+x

504       CONTINUE

          DO 505 INDXA=INDXC,NUMA
           IPOS=INDXF(INDXA,INDXB,NUMA)+IOFFIR
           X=BUF(IPOS)*FABCI
           INDA=INDXA+OCCA
           GINDXL=INDXT(INDC,INDA)

           gamma(indxf(gindxl,gindxr,maxl))=
     &        gamma(indxf(gindxl,gindxr,maxl))+x

           gamma(indxf(gindxr,gindxl,maxl))=
     &        gamma(indxf(gindxr,gindxl,maxl))+x

505       CONTINUE
503      CONTINUE
502     CONTINUE
501    CONTINUE
C
      ELSE
C
C NOW PROCESS THE G(Ab,Ci) LIST
C
C  UHF VERSION : THE ELEMENTS ARE PLACED IN G(AC,bi)
C
      CALL SYMPOS('F',VRT(1,1),POP(1,2),IRREPDO,1,ILOGREC)
      CALL SYMPOS('F',VRT(1,1),VRT(1,2),IRREPDO,1,IOFFIR)
       NUMI=POP(IRREPDO,2)
       NUMC=VRT(IRREPDO,1)
       NUMB=VRT(IRREPDO,2)
       NUMA=VRT(IRREPDO,1)
       NUMAO=MOPOP(IRREPDO)
       OCCA=POP(IRREPDO,1)
       OCCB=POP(IRREPDO,2)
       MAXL=NNP1O2(NUMAO)
       DO 1501 INDXI=1,NUMI
        DO 1502 INDXC=1,NUMC
         ILOGREC=ILOGREC+1
         INDI=INDXI
         INDC=INDXC+OCCA
         CALL GETLST(BUF,ILOGREC,1,1,1,130)
         DO 1503 INDXB=1,NUMB
          INDB=INDXB+OCCB
          GINDXR=INDXT(INDI,INDB)
          DO 1504 INDXA=1,INDXC
           IPOS=INDXF(INDXA,INDXB,NUMA)+IOFFIR
           X=BUF(IPOS)*FABCI
           INDA=INDXA+OCCA
           GINDXL=INDXT(INDA,INDC)

           gamma(indxf(gindxl,gindxr,maxl))=
     &        gamma(indxf(gindxl,gindxr,maxl))+x

1504       CONTINUE
          DO 1505 INDXA=INDXC,NUMA
           IPOS=INDXF(INDXA,INDXB,NUMA)+IOFFIR
           X=BUF(IPOS)*FABCI
           INDA=INDXA+OCCA
           GINDXL=INDXT(INDC,INDA)

           gamma(indxf(gindxl,gindxr,maxl))=
     &        gamma(indxf(gindxl,gindxr,maxl))+x

1505       CONTINUE
1503      CONTINUE
1502     CONTINUE
1501    CONTINUE
C
C NOW PROCESS THE G(Ab,Ic) LIST
C
      CALL SYMPOS('F',VRT(1,2),POP(1,1),IRREPDO,1,ILOGREC)
      CALL SYMPOS('F',VRT(1,1),VRT(1,2),IRREPDO,1,IOFFIR)
       NUMI=POP(IRREPDO,1)
       NUMC=VRT(IRREPDO,2)
       NUMB=VRT(IRREPDO,2)
       NUMA=VRT(IRREPDO,1)
       NUMAO=MOPOP(IRREPDO)
       OCCA=POP(IRREPDO,1)
       OCCB=POP(IRREPDO,2)
       MAXL=NNP1O2(NUMAO)
       DO 2501 INDXC=1,NUMC
        INDC=INDXC+OCCB
        DO 2502 INDXI=1,NUMI
         ILOGREC=ILOGREC+1
         INDI=INDXI
         CALL GETLST(BUF,ILOGREC,1,1,1,129)
         DO 2503 INDXB=1,INDXC
          INDB=INDXB+OCCB
          GINDXR=INDXT(INDB,INDC)
          DO 2504 INDXA=1,NUMA
           IPOS=INDXF(INDXA,INDXB,NUMA)+IOFFIR
           X=BUF(IPOS)*FABCI
           INDA=INDXA+OCCA
           GINDXL=INDXT(INDI,INDA)

           gamma(indxf(gindxl,gindxr,maxl))=
     &        gamma(indxf(gindxl,gindxr,maxl))+x

2504       CONTINUE
2503      CONTINUE
          DO 2505 INDXB=INDXC,NUMB
           INDB=INDXB+OCCB
           GINDXR=INDXT(INDC,INDB)
           DO 2506 INDXA=1,NUMA
           IPOS=INDXF(INDXA,INDXB,NUMA)+IOFFIR
           X=BUF(IPOS)*FABCI
           INDA=INDXA+OCCA
           GINDXL=INDXT(INDI,INDA)

           gamma(indxf(gindxl,gindxr,maxl))=
     &        gamma(indxf(gindxl,gindxr,maxl))+x

2506       CONTINUE
2505      CONTINUE
2502     CONTINUE
2501    CONTINUE
C
       ENDIF
C
       IF(IUHF.EQ.0) THEN
C
C NOW PROCESS THE G(Ij,Ka) LIST
C
      CALL SYMPOS('F',VRT(1,2),POP(1,1),IRREPDO,1,ILOGREC)
      CALL SYMPOS('F',POP(1,1),POP(1,2),IRREPDO,1,IOFFIR)
       NUMA=VRT(IRREPDO,2)
       NUMK=POP(IRREPDO,1)
       NUMJ=POP(IRREPDO,2)
       NUMI=POP(IRREPDO,1)
       NUMAO=MOPOP(IRREPDO)
       OCCA=POP(IRREPDO,1)
       OCCB=POP(IRREPDO,2)
       MAXL=NNP1O2(NUMAO)
       DO 601 INDXA=1,NUMA
        DO 602 INDXK=1,NUMK
         ILOGREC=ILOGREC+1
         INDA=INDXA+OCCB
         INDK=INDXK
         CALL GETLST(BUF,ILOGREC,1,1,1,110)
         DO 603 INDXJ=1,NUMJ
          INDJ=INDXJ
          GINDXR=INDXT(INDJ,INDA)
          DO 604 INDXI=1,INDXK
           IPOS=INDXF(INDXI,INDXJ,NUMI)+IOFFIR
           X=BUF(IPOS)*FIJKA
           INDI=INDXI
           GINDXL=INDXT(INDI,INDK)

           gamma(indxf(gindxl,gindxr,maxl))=
     &        gamma(indxf(gindxl,gindxr,maxl))+x


           gamma(indxf(gindxr,gindxl,maxl))=
     &        gamma(indxf(gindxr,gindxl,maxl))+x

604       CONTINUE
          DO 605 INDXI=INDXK,NUMI
           IPOS=INDXF(INDXI,INDXJ,NUMI)+IOFFIR
           X=BUF(IPOS)*FIJKA
           INDI=INDXI
           GINDXL=INDXT(INDK,INDI)

           gamma(indxf(gindxl,gindxr,maxl))=
     &        gamma(indxf(gindxl,gindxr,maxl))+x


           gamma(indxf(gindxr,gindxl,maxl))=
     &        gamma(indxf(gindxr,gindxl,maxl))+x

605       CONTINUE
603      CONTINUE
602     CONTINUE
601    CONTINUE
C
      ELSE
C
C NOW PROCESS THE G(Ij,Ka) LIST
C
C  UHF VERSION
C
      CALL SYMPOS('F',VRT(1,2),POP(1,1),IRREPDO,1,ILOGREC)
      CALL SYMPOS('F',POP(1,1),POP(1,2),IRREPDO,1,IOFFIR)
       NUMA=VRT(IRREPDO,2)
       NUMK=POP(IRREPDO,1)
       NUMJ=POP(IRREPDO,2)
       NUMI=POP(IRREPDO,1)
       NUMAO=MOPOP(IRREPDO)
       OCCA=POP(IRREPDO,1)
       OCCB=POP(IRREPDO,2)
       MAXL=NNP1O2(NUMAO)
       DO 1601 INDXA=1,NUMA
        DO 1602 INDXK=1,NUMK
         ILOGREC=ILOGREC+1
         INDA=INDXA+OCCB
         INDK=INDXK
         CALL GETLST(BUF,ILOGREC,1,1,1,110)
         DO 1603 INDXJ=1,NUMJ
          INDJ=INDXJ
          GINDXR=INDXT(INDJ,INDA)
          DO 1604 INDXI=1,INDXK
           IPOS=INDXF(INDXI,INDXJ,NUMI)+IOFFIR
           X=BUF(IPOS)*FIJKA
           INDI=INDXI
           GINDXL=INDXT(INDI,INDK)

           gamma(indxf(gindxl,gindxr,maxl))=
     &        gamma(indxf(gindxl,gindxr,maxl))+x


1604       CONTINUE
          DO 1605 INDXI=INDXK,NUMI
           IPOS=INDXF(INDXI,INDXJ,NUMI)+IOFFIR
           X=BUF(IPOS)*FIJKA
           INDI=INDXI
           GINDXL=INDXT(INDK,INDI)

           gamma(indxf(gindxl,gindxr,maxl))=
     &        gamma(indxf(gindxl,gindxr,maxl))+x


1605       CONTINUE
1603      CONTINUE
1602     CONTINUE
1601    CONTINUE
C
C NOW PROCESS THE G(Ij,Ak) LIST
C
      CALL SYMPOS('F',VRT(1,1),POP(1,2),IRREPDO,1,ILOGREC)
      CALL SYMPOS('F',POP(1,1),POP(1,2),IRREPDO,1,IOFFIR)
       NUMA=VRT(IRREPDO,1)
       NUMK=POP(IRREPDO,2)
       NUMJ=POP(IRREPDO,2)
       NUMI=POP(IRREPDO,1)
       NUMAO=MOPOP(IRREPDO)
       OCCA=POP(IRREPDO,1)
       OCCB=POP(IRREPDO,2)
       MAXL=NNP1O2(NUMAO)
       DO 2601 INDXK=1,NUMK
        DO 2602 INDXA=1,NUMA
         ILOGREC=ILOGREC+1
         INDA=INDXA+OCCA
         INDK=INDXK
         CALL GETLST(BUF,ILOGREC,1,1,1,109)
         DO 2603 INDXJ=1,INDXK
          INDJ=INDXJ
          GINDXR=INDXT(INDJ,INDK)
          DO 2604 INDXI=1,NUMI
           IPOS=INDXF(INDXI,INDXJ,NUMI)+IOFFIR
           X=BUF(IPOS)*FIJKA
           INDI=INDXI
           GINDXL=INDXT(INDI,INDA)

           gamma(indxf(gindxl,gindxr,maxl))=
     &        gamma(indxf(gindxl,gindxr,maxl))+x


2604       CONTINUE
2603      CONTINUE
          DO 2606 INDXJ=INDXK,NUMJ
           INDJ=INDXJ
           GINDXR=INDXT(INDK,INDJ)
          DO 2605 INDXI=1,NUMI
           IPOS=INDXF(INDXI,INDXJ,NUMI)+IOFFIR
           X=BUF(IPOS)*FIJKA
           INDI=INDXI
           GINDXL=INDXT(INDI,INDA)

           gamma(indxf(gindxl,gindxr,maxl))=
     &        gamma(indxf(gindxl,gindxr,maxl))+x


2605       CONTINUE
2606      CONTINUE
2602     CONTINUE
2601    CONTINUE
C
       ENDIF
C
      RETURN
      END
