      SUBROUTINE LARGEAMP(T2,MAXCOR,IUHF)
C
C LOCATE LARGEST ELEMENTS IN AMPLITUDE VECTORS AND PRINT THEM OUT.
C
      IMPLICIT INTEGER (A-Z)
      LOGICAL RHF
      CHARACTER*4 SPCASE(4)
      DOUBLE PRECISION T2(1)
      DIMENSION IRROFF(8),NUM(8),DSZ(8)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPX(255,2),DIRPRD(8,8)
      COMMON /FSSYMPOP/ FSDPDAN(8,22),FSDPDNA(8,22),FSDPDAA(8,22),
     &                  FSDPDIN(8,22),FSDPDNI(8,22),FSDPDII(8,22),
     &                  FSDPDAI(8,22),FSDPDIA(8,22)
      COMMON /FSSYM/ POPA(8,2),POPI(8,2),VRTA(8,2),VRTI(8,2),
     &               NTAN(2),NTNA(2),NTAA(2),
     &               NTIN(2),NTNI(2),NTII(2),
     &               NTAI(2),NTIA(2),
     &               NF1AN(2),NF1AA(2),NF1IN(2),NF1II(2),NF1AI(2),
     &               NF2AN(2),NF2AA(2),NF2IN(2),NF2II(2),NF2AI(2)
      DATA ZILCH  /0.0/
      DATA SPCASE /'AAAA','BBBB','BABA','ABAB'/
C
      RHF=IUHF.EQ.0
      IBOT=4-3*IUHF
      ntop=10
      WRITE(6,2000)
2000  FORMAT(T12,47('-'))
2001  FORMAT(T27,14('-'))
      WRITE(6,1000)
      WRITE(6,2000)
1000  FORMAT(T15,' Summary of largest excitation amplitudes ')
C
      DO 100 ISPIN=IBOT,4       
       WRITE(6,1001)SPCASE(ISPIN)
1001   FORMAT(T15,'            Spin case ',A4)
       WRITE(6,2001)
       IF(ISPIN.GT.2)THEN
        WRITE(6,2002)
2002    FORMAT(T21,'-',T31,'-')
       ELSEIF(ISPIN.EQ.2)THEN
        WRITE(6,2003)
2003    FORMAT(T16,'-',T21,'-',T26,'-',T31,'-')
       ENDIF
       IF(ISPIN.NE.3)WRITE(6,2004)
2004   FORMAT(T16,'i',T21,'j',T26,'k',T31,'a')
       IF(ISPIN.EQ.3)WRITE(6,2005)
2005   FORMAT(T16,'i',T21,'j',T26,'a',T31,'k')
C
C ABAB SPIN CASE
C
       IOFF=1
       DO 10 IRREP=1,NIRREP
        LISTT=95+ISPIN
        IF(ISPIN.NE.3)THEN
         NUM(IRREP)=FSDPDAN(IRREP,ISYTYP(2,LISTT))
         DSZ(IRREP)=IRPDPD (IRREP,ISYTYP(1,LISTT))
         CALL FSGET(T2(IOFF),1,NUM(IRREP),1,IRREP,LISTT,'NNAN')
        ELSE
         NUM(IRREP)=FSDPDNA(IRREP,ISYTYP(2,LISTT))
         DSZ(IRREP)=IRPDPD (IRREP,ISYTYP(1,LISTT))
         CALL FSGET(T2(IOFF),1,NUM(IRREP),1,IRREP,LISTT,'NNNA')
        ENDIF
        IRROFF(IRREP)=IOFF
        IOFF=IOFF+NUM(IRREP)*DSZ(IRREP)
10     CONTINUE
       NSIZE=IOFF-1
       DO 11 ILOCATE=1,MIN(NTOP,NSIZE)
        J=ISAMAX(NSIZE,T2,1)
        IF(ISPIN.LE.2)THEN
         call abs2pqrs(j,num,dsz,'PCK','FUL',
     &                 pop(1,ispin),pop(1,ispin),
     &                 popa(1,ispin),vrt(1,ispin),
     &                 indxi,indxj,indxk,indxa)
         write(6,'(T12,4I5,F20.10)')indxi,indxj,indxk,indxa,t2(j)
        ELSEIF(ISPIN.EQ.3)THEN
         call abs2pqrs(j,num,dsz,'FUL','FUL',
     &                 pop(1,1),pop(1,2),
     &                 vrt(1,1),popa(1,2),
     &                 indxi,indxj,indxk,indxa)
         write(6,'(T12,4I5,F20.10)')indxi,indxj,indxk,indxa,t2(j)
        ELSEIF(ISPIN.EQ.4)THEN
         call abs2pqrs(j,num,dsz,'FUL','FUL',
     &                 pop(1,1),pop(1,2),popa(1,1),vrt(1,2),
     &                 indxi,indxj,indxk,indxa)
         write(6,'(T12,4I5,F20.10)')indxi,indxj,indxk,indxa,t2(j)
        ENDIF
        T2(J)=ZILCH
11     CONTINUE
100   CONTINUE
      WRITE(6,2000)
      RETURN
      END
