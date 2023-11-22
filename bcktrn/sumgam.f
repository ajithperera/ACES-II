      SUBROUTINE SUMGAM
C
C CHECKSUMS THE AO GAMMAS
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION Z1,Z2,Z3,Z4
      DIMENSION IDID(8)
      COMMON // ICORE(1) 
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /AOPOPS/ AOPOP(8),MOPOP(8),NAO,NAOSPH,NMO
      COMMON /SZAOGM/ AOGMSZ(2,4,100)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /ISTART/ I0,ICRSIZ 
C
      NNP1O2(I)=(I*(I+1))/2
C
C INITIALIZE NEW MOINTS FILE.
C
C
C CREATE AAAA LISTS.
C
        izilch=0
        if(izilch.eq.0)return
      Z1=0.0
      DO 10 IRREP=1,NIRREP
       NUMDIS=NNP1O2(AOPOP(IRREP))
       DISSIZ=NUMDIS
       CALL GETLST(ICORE(I0),1,NUMDIS,1,1,IRREP)
       if(irrep.eq.1)then
        write(13,'((I5,f20.10))')(j+1,icore(i0+j),j=0,1520,1)
       endif
       CALL ZNRM(ICORE(I0),NUMDIS*DISSIZ,Z1)
      write(6,'(F15.9)')snrm2(numdis*dissiz,icore(i0),1)
10    CONTINUE
      WRITE(6,100)Z1
100   FORMAT(T3,'@SUMGAM-I, AAAA norm of AO gamma: ',F20.10,'.')
C
C CREATE AABB LISTS.
C
      Z2=0.0
      ITHRU=0
      DO 20 IRREPA=2,NIRREP
       DO 21 IRREPB=1,IRREPA-1
        ITHRU=ITHRU+1
        DISSIZ=NNP1O2(AOPOP(IRREPA))
        NUMDIS=NNP1O2(AOPOP(IRREPB))
        CALL GETLST(ICORE(I0),1,NUMDIS,1,2,ITHRU)
      write(6,'(2I5,F15.9)')irrepa,irrepb,snrm2(numdis*dissiz,
     &      icore(i0),1)
        CALL ZNRM(ICORE(I0),NUMDIS*DISSIZ,Z2)
c         write(6,*)z2
21     CONTINUE
20    CONTINUE
      WRITE(6,200)Z2
200   FORMAT(T3,'@SUMGAM-I, AABB norm of AO gamma: ',F20.10,'.')
C
C CREATE ABAB LISTS.
C
      ITHRU=0
      Z3=0.0
      DO 30 IRREP=2,NIRREP
       DO 31 IRREPA=1,NIRREP
        IRREPB=DIRPRD(IRREPA,IRREP)
        IF(IRREPB.LT.IRREPA)GOTO 31
        ITHRU=ITHRU+1
        NUMDIS=AOPOP(IRREPA)*AOPOP(IRREPB)
        DISSIZ=NUMDIS
        CALL GETLST(ICORE(I0),1,NUMDIS,1,3,ITHRU)
      write(6,'(2I5,F15.9)')irrepa,irrepb,snrm2(numdis*dissiz,
     &      icore(i0),1)
c        CALL ZNRM(ICORE(I0),NUMDIS*DISSIZ,Z3)
c         write(6,*)z3
31     CONTINUE
30    CONTINUE
      WRITE(6,300)Z3
300   FORMAT(T3,'@SUMGAM-I, ABAB norm of AO gamma: ',F20.10,'.')
C
C CREATE ABCD LISTS.
C
      ITHRU=0
      Z4=0.0
      DO 40 IRREP=2,NIRREP
       DO 41 IRREPD=1,NIRREP
        IRREPC=DIRPRD(IRREPD,IRREP)
        IF(IRREPC.LT.IRREPD)GOTO 41
        IBOT=MAX(IRREPC,IRREPD)+1
        CALL IZERO(IDID,8)
        DO 42 ITMP=IBOT,NIRREP 
         IRREPA=DIRPRD(ITMP,IRREP)
         IRREPB=MIN(IRREPA,ITMP)      
         IRREPA=MAX(IRREPA,ITMP)
         IF(MAX(IDID(IRREPA),IDID(IRREPB)).NE.0)GOTO 42
         IDID(IRREPA)=1
         IDID(IRREPB)=1
         ITHRU=ITHRU+1
         DISSIZ=AOPOP(IRREPC)*AOPOP(IRREPD)
         NUMDIS=AOPOP(IRREPA)*AOPOP(IRREPB)
         CALL GETLST(ICORE(I0),1,NUMDIS,1,4,ITHRU)
      write(6,'(4I5,F15.9)')irrepa,irrepb,irrepc,irrepd,
     &      snrm2(numdis*dissiz,
     &      icore(i0),1)
c         CALL ZNRM(ICORE(I0),NUMDIS*DISSIZ,Z4)
c          write(6,*)' CONTRIBUTION: ',Z4
42      CONTINUE
41     CONTINUE
40    CONTINUE
      WRITE(6,400)Z4
400   FORMAT(T3,'@SUMGAM-I, ABCD norm of AO gamma: ',F20.10,'.')
C
      WRITE(6,1000)Z1+2.0*Z2+Z3+2.0*Z4
1000  FORMAT(T3,'@SUMGAM-I, Total norm of AO gamma: ',F20.10,'.')
      RETURN
      END
