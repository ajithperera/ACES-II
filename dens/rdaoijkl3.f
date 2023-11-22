      SUBROUTINE RDAOIJKL3(G6,XIA,BUF,IBUF,IAOSYM,IMAP,POP,
     &                     IlNBUF,LUINT,IUHF,NAO,IRREPX)
C
C THIS ROUTINE LOADS THE AO INTEGRALS FROM THE IJKL FILE (ALL FOUR
C INDICES HAVE DIFFERENT SYMMETRIES) AND CONTRACTS THEM WITH THE G6
C AMPLITUDES.  SOME COMPLICATED STUFF!
C
CEND
C
C JANUARY 94/JG,  KARLSRUHE AND AUSTIN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL DONE
      INTEGER AND,OR,DIRPRD,POP
      DIMENSION BUF(ILNBUF),IBUF(ILNBUF)
      DIMENSION XIA(*),G6(*),IAOSYM(*),IMAP(*)
      DIMENSION IOFFR(8,8),IOFFL(8,8),IOFFG6(8)
      DIMENSION IOFFMO(8),POP(8),IOFFX(8),ISIZE(8)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     &             IRPDPDAO(8),IRPDPDAOS(8),ISTART(8,8),ISTARTMO(8,3)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/G6OOC/IPOPS(8,8),IPOPE(8,8),NDONE(8,8),NSTARTM,
     &             ISTARTI,NSIZG,IOFFSET,DONE
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW 
C
      IUPKI(INT)=AND(INT,IALONE)
      IUPKJ(INT)=AND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=AND(ISHFT(INT,-2*IBITWD),IALONE)
      IUPKL(INT)=AND(ISHFT(INT,-3*IBITWD),IALONE)
      IPACK(I,J,K,L)=OR(OR(OR(I,ISHFT(J,IBITWD)),ISHFT(K,2*IBITWD)),
     &                  ISHFT(L,3*IBITWD))
      INDX(I,J)=J+(I*(I-1))/2
      INDX3(I,J)=I+(J*(J-1))/2
      INDX2(I,J,N)=I+(J-1)*N
C
      DATA HALF /0.5D0/
C
C GET OFFSETS FOR XIA ARRAY (SYMMETRY IRREPX)
C
      IOFFX1=0
      DO 9 IRREPR=1,NIRREP
       IRREPL=DIRPRD(IRREPR,IRREPX)
       IOFFX(IRREPR)=IOFFX1
       IOFFX1=IOFFX1+POP(IRREPR)*IAOPOP(IRREPL)
9     CONTINUE
C
C FILL IOFFR,IOFFL,AND IOFFG6
C
      DO 12 IRREP=1,NIRREP
       IOFFRR=0
       IOFFLL=0
       DO 11 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREP)
        IOFFR(IRREPR,IRREP)=IOFFRR
        IOFFL(IRREPR,IRREP)=IOFFLL
        IOFFRR=IOFFRR+POP(IRREPR)*IAOPOP(IRREPL)
        IOFFLL=IOFFLL+IAOPOP(IRREPR)*IAOPOP(IRREPL)
11     CONTINUE
       ISIZE(IRREP)=IOFFRR
12    CONTINUE
      IOFFG=-IOFFSET
      DO 13 IRREPR=1,NIRREP
       IRREPL=DIRPRD(IRREPX,IRREPR)
       IOFFG6(IRREPR)=IOFFG
       IOFFG=IOFFG+ISIZE(IRREPR)*IRPDPDAO(IRREPL)
13    CONTINUE
C
C FILL IMAP ARRAY 
C
      IOFF=0
      DO 14 IRREP=1,NIRREP
       ITHRU=0
       DO 15 IAO=1,IAOPOP(IRREP)
        ITHRU=ITHRU+1
        IOFF=IOFF+1
        IMAP(IOFF)=ITHRU
15     CONTINUE
14    CONTINUE
C
C START READING THE AO 2E INTEGRALS
C
      TOTAL=0.0D0
      CALL TIMER(1)
      NAOBUF=0
      NUMINT=0
      CALL LOCATE(LUINT,'TWOELSUP')
C
C REENTRY POINT FOR LOOP OVER ALL INTEGRAL BUFFERS
C
1     READ(LUINT)BUF,IBUF,NUT
      NUT8=8*NUT
      NAOBUF=NAOBUF+1
C
C PROCESS THE AO 2E INTEGRALS
C
      DO 10 INT=1,NUT
C
C X IS THE INTEGRAL VALUE
C
       X=BUF(INT)
C
C INDI,INDJ,INDK,INDL ARE THE INDICES OF (IJ|KL) = <IK|JL>
C     
       INDI=IUPKI(IBUF(INT))
       INDJ=IUPKJ(IBUF(INT))
       INDK=IUPKK(IBUF(INT))
       INDL=IUPKL(IBUF(INT))
C
       IJ=INDX(MAX(INDJ,INDI),MIN(INDJ,INDI))
       KL=INDX(MAX(INDL,INDK),MIN(INDL,INDK))
C
C GET THE PREFACTORS FOR THE INTEGRALS
C
       IF(INDK.EQ.INDL) THEN
        X=X*HALF
       ENDIF
       IF(INDI.EQ.INDJ) THEN
        X=X*HALF
       ENDIF
       IF(IJ.EQ.KL) THEN
        X=X*HALF
       ENDIF
C
C SYMMETRIES OF AO INTEGRALS
C
       ISYMI=IAOSYM(INDI)
       ISYMJ=IAOSYM(INDJ)
       ISYMK=IAOSYM(INDK)
       ISYML=IAOSYM(INDL)
       ISYMIX=DIRPRD(ISYMI,IRREPX)
       ISYMJX=DIRPRD(ISYMJ,IRREPX)
       ISYMKX=DIRPRD(ISYMK,IRREPX)
       ISYMLX=DIRPRD(ISYML,IRREPX)
C
C NUMBER OF AOS IN THE CORRESPONDING IRREPS
C
       NAOI=IAOPOP(ISYMI)
       NAOJ=IAOPOP(ISYMJ)
       NAOK=IAOPOP(ISYMK)
       NAOL=IAOPOP(ISYML)
C
C OFFSETS FOR XIA ARRAY IN MIXED REPRESENTATION
C
       IOFFXI=IOFFX(ISYMIX)
       IOFFXJ=IOFFX(ISYMJX)
       IOFFXK=IOFFX(ISYMKX)
       IOFFXL=IOFFX(ISYMLX)
C
       I=IMAP(INDI)
       J=IMAP(INDJ)
       K=IMAP(INDK)
       L=IMAP(INDL)
C
       ISYMIK=DIRPRD(IAOSYM(INDI),IAOSYM(INDK))
       ISYMIL=DIRPRD(IAOSYM(INDI),IAOSYM(INDL))
       ISYMJK=DIRPRD(IAOSYM(INDJ),IAOSYM(INDK))
       ISYMJL=DIRPRD(IAOSYM(INDJ),IAOSYM(INDL))
       ISYMIKX=DIRPRD(ISYMIK,IRREPX)
       ISYMILX=DIRPRD(ISYMIL,IRREPX)
       ISYMJKX=DIRPRD(ISYMJK,IRREPX)
       ISYMJLX=DIRPRD(ISYMJL,IRREPX)
C
       IOFFG6IK=IOFFG6(ISYMIKX)
       IOFFG6IL=IOFFG6(ISYMILX)
       IOFFG6JK=IOFFG6(ISYMJKX)
       IOFFG6JL=IOFFG6(ISYMJLX)
C
       NSIZIK=IRPDPDAO(ISYMIK)
       NSIZIL=IRPDPDAO(ISYMIL)
       NSIZJK=IRPDPDAO(ISYMJK)
       NSIZJL=IRPDPDAO(ISYMJL)

       IOFFRIK=IOFFR(ISYMIX,ISYMIKX)
       IOFFRJK=IOFFR(ISYMJX,ISYMJKX)
       IOFFRIL=IOFFR(ISYMIX,ISYMILX)
       IOFFRJL=IOFFR(ISYMJX,ISYMJLX)
       IOFFRKI=IOFFR(ISYMKX,ISYMIKX)
       IOFFRLI=IOFFR(ISYMLX,ISYMILX)
       IOFFRKJ=IOFFR(ISYMKX,ISYMJKX)
       IOFFRLJ=IOFFR(ISYMLX,ISYMJLX)
C
       IOFFLIK=IOFFL(ISYMI,ISYMIK)
       IOFFLJK=IOFFL(ISYMJ,ISYMJK)
       IOFFLIL=IOFFL(ISYMI,ISYMIL)
       IOFFLJL=IOFFL(ISYMJ,ISYMJL)
       IOFFLKI=IOFFL(ISYMK,ISYMIK)
       IOFFLLI=IOFFL(ISYML,ISYMIL)
       IOFFLKJ=IOFFL(ISYMK,ISYMJK)
       IOFFLLJ=IOFFL(ISYML,ISYMJL)
C
C FIRST CONTRIBUTION: X(IOCCI,I) = G(L,J,K,IOCCI) * I(L,J,K,I)
C
       DO 101 IOCCI=IPOPS(ISYMIX,ISYMIKX),IPOPE(ISYMIX,ISYMIKX)
        IADR1=IOFFG6IK+NSIZJL*(IOFFRIK+(IOCCI-1)*NAOK+K-1)
     &         +IOFFLJL+L+(J-1)*NAOL
        XIA(I+(IOCCI-1)*NAOI+IOFFXI)=XIA(I+(IOCCI-1)*NAOI+IOFFXI)+
     &    X*G6(IADR1)
101    CONTINUE
C
C SECOND CONTRIBUTION: X(IOCCI,I) = G(K,J,L,IOCCI) * I(K,J,L,I)
C
       DO 102 IOCCI=IPOPS(ISYMIX,ISYMILX),IPOPE(ISYMIX,ISYMILX)
        IADR2=IOFFG6IL+NSIZJK*(IOFFRIL+(IOCCI-1)*NAOL+L-1)
     &         +IOFFLJK+K+(J-1)*NAOK
        XIA(I+(IOCCI-1)*NAOI+IOFFXI)=XIA(I+(IOCCI-1)*NAOI+IOFFXI)+
     &    X*G6(IADR2)
102    CONTINUE
C
C THIRD CONTRIBUTION: X(IOCCJ,J) = G(L,I,K,IOCCJ) * I(L,I,K,J)
C
       DO 103 IOCCJ=IPOPS(ISYMJX,ISYMJKX),IPOPE(ISYMJX,ISYMJKX)
        IADR3=IOFFG6JK+NSIZIL*(IOFFRJK+(IOCCJ-1)*NAOK+K-1)
     &         +IOFFLIL+L+(I-1)*NAOL
        XIA(J+(IOCCJ-1)*NAOJ+IOFFXJ)=XIA(J+(IOCCJ-1)*NAOJ+IOFFXJ)+
     &    X*G6(IADR3)
103    CONTINUE
C
C FOURTH CONTRIBUTION: X(IOCCJ,J) = G(K,I,L,IOCCJ) * I(K,I,L,J)
C
       DO 104 IOCCJ=IPOPS(ISYMJX,ISYMJLX),IPOPE(ISYMJX,ISYMJLX)
        IADR4=IOFFG6JL+NSIZIK*(IOFFRJL+(IOCCJ-1)*NAOL+L-1)
     &         +IOFFLIK+K+(I-1)*NAOK
        XIA(J+(IOCCJ-1)*NAOJ+IOFFXJ)=XIA(J+(IOCCJ-1)*NAOJ+IOFFXJ)+
     &    X*G6(IADR4)
104    CONTINUE
C
C FIFTH CONTRIBUTION: X(IOCCK,K) = G(J,L,I,IOCCK) * I(J,L,I,K)
C
       DO 105 IOCCK=IPOPS(ISYMKX,ISYMIKX),IPOPE(ISYMKX,ISYMIKX)
        IADR5=IOFFG6IK+NSIZJL*(IOFFRKI+(IOCCK-1)*NAOI+I-1)
     &         +IOFFLLJ+J+(L-1)*NAOJ
        XIA(K+(IOCCK-1)*NAOK+IOFFXK)=XIA(K+(IOCCK-1)*NAOK+IOFFXK)+
     &    X*G6(IADR5)
105    CONTINUE
C
C SIXTH CONTRIBUTION: X(IOCCL,L) = G(J,K,I,IOCCL) * I(J,K,I,L)
C
       DO 106 IOCCL=IPOPS(ISYMLX,ISYMILX),IPOPE(ISYMLX,ISYMILX)
        IADR6=IOFFG6IL+NSIZJK*(IOFFRLI+(IOCCL-1)*NAOI+I-1)
     &         +IOFFLKJ+J+(K-1)*NAOJ
        XIA(L+(IOCCL-1)*NAOL+IOFFXL)=XIA(L+(IOCCL-1)*NAOL+IOFFXL)+
     &    X*G6(IADR6)
106    CONTINUE
C
C SEVENTH CONTRIBUTION: X(IOCCK,K) = G(I,L,J,IOCCK) * I(I,L,J,K)
C
       DO 107 IOCCK=IPOPS(ISYMKX,ISYMJKX),IPOPE(ISYMKX,ISYMJKX)
        IADR7=IOFFG6JK+NSIZIL*(IOFFRKJ+(IOCCK-1)*NAOJ+J-1)
     &         +IOFFLLI+I+(L-1)*NAOI
        XIA(K+(IOCCK-1)*NAOK+IOFFXK)=XIA(K+(IOCCK-1)*NAOK+IOFFXK)+
     &    X*G6(IADR7)
107    CONTINUE
C
C EIGHTH CONTRIBUTION: X(IOCCL,L) = G(I,K,J,IOCCL) * I(I,K,J,L)
C
       DO 108 IOCCL=IPOPS(ISYMLX,ISYMJLX),IPOPE(ISYMLX,ISYMJLX)
        IADR8=IOFFG6JL+NSIZIK*(IOFFRLJ+(IOCCL-1)*NAOJ+J-1)
     &         +IOFFLKI+I+(K-1)*NAOI
        XIA(L+(IOCCL-1)*NAOL+IOFFXL)=XIA(L+(IOCCL-1)*NAOL+IOFFXL)+
     &    X*G6(IADR8)
108    CONTINUE
10    CONTINUE
C
      NUMINT=NUMINT+MAX(NUT,0)
C
      IF(NUT.NE.-1)GOTO 1
C
C ALL DONE, PRINT OUT TIMING INFORMATION
C
      CALL TIMER(1)
C
c      write(*,*) 'CPU time for the contraction : ',TIMENEW
c      write(*,*) 'Total CPU time : ', TIMENEW
C
c      WRITE(6,*)' processed ',numint,' ao basis integrals ',
c     &          'from ',naobuf-1,' buffers.'
C
c
      RETURN
      END