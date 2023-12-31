C
C DRIVER FOR AO-BASED PARTICLE-PARTICLE LADDER CONTRACTION USING
C MULTIPASS ALGORITHM.
C
      SUBROUTINE AOLAD3(ICORE,MAXCOR,IUHF,TAU,IRREPX,LSTAO,LSTAOINC)
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION X,SNRM2,SDOT,ABLAD,ZILCH,ONE
      LOGICAL TAU,PRINT
      CHARACTER*2 SPLAB(3)
      CHARACTER*80 FNAME
C
      DIMENSION ICORE(MAXCOR)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     &             IRPDPDAO(8),IRPDPDAOS(8),ISTART(8,8),ISTARTMO(8,3)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/FLAGS/IFLAGS(100)
C
      DATA ZILCH,ONE /0.0D0,1.0D0/
      DATA SPLAB /'AA','BB','AB'/
C
      NNM1O2(I)=(I*(I-1))/2
C
      IONE=1
      LUINT=10
      NMO=NOCCO(1)+NVRTO(1)
      CALL GETREC(20,'JOBARC','NBASTOT ',IONE,NAO)
      CALL GETAOINF(IUHF,IRREPX)
      ILNBUF=600
      NAOSQ=NAO*NAO
      NAOTRI=(NAO*(NAO-1))/2
C
C CALCULATE MAXIMUM SIZE OF T2 IRREP FOR MEMORY ALLOCATION
C
      ISIZT2MAX=0
      DO 15 ISPIN=3,3-2*IUHF,-1
       DO 16 IRREPIJ=1,NIRREP
        IRREPAB=DIRPRD(IRREPIJ,IRREPX)
        ISIZT2MAX=MAX(ISIZT2MAX,IRPDPD(IRREPIJ,ISYTYP(2,43+ISPIN))
     &                     *IRPDPDAO(IRREPAB))
16     CONTINUE
15    CONTINUE       
C
      PRINT=IFLAGS(1).GE.10
      I000=1
      I010=I000+ILNBUF*IINTFP*8
      I020=I010+ILNBUF*8
      I030=I020+ILNBUF*8
      I040=I030+ILNBUF*8
      I050=I040+MAX(ILNBUF*8,NAO)
      I051=I050+NAO*NAO
      I060=I051+NAO*NAO
      I070=I060+NAO*NAO
      I080=MAX(I070+NAO*NAO,I000+ISIZT2MAX*2*IINTFP)
      I090=I080+8
      I100=I090+8
      I110=I100+8
      I111=I110+8
      I112=I111+NAOSQ
      I120=I112+NAOSQ+IINTFP
      NSIZE=MAXCOR-I120
C
      ioff=i111
      do 104 i=1,nao
       do 105 j=1,nao
        icore(ioff)=i+(j-1)*nao
        ioff=ioff+1
105    continue
104   continue 
      ioff=i112
      do 106 i=1,nao
       do 107 j=1,nao
        ido=max(i,j)
        jdo=min(i,j)
        icore(ioff)=jdo+(ido*(ido-1))/2
        ioff=ioff+1
107    continue
106   continue  
C      
C      
C MULTIPLE PASS LOGIC BEGINS HERE
C
      CALL MPINFO1(ICORE(I100),ICORE(I110),ITOTSIZ,ISPIN)
c
      IRRBEGIN=1
      ILOWADR=1
C
      MXCOR=MAXCOR-I120
C
      NPASS=1
      ILOWADR=1
      IRRWRT=1
      IDIS1=1
      DISLEFT=IRPDPDAO(IRRWRT)-IDIS1+1
      NSIZE=MXCOR
      IF(PRINT)THEN
       WRITE(6,6000)ITOTSIZ
       WRITE(6,6001)MXCOR
       WRITE(6,6002)1+(ITOTSIZ*IINTFP)/MXCOR
6000   FORMAT(T3,'@AOLAD3-I, Total number of integrals : ',I15,'.')
6001   FORMAT(T3,'@AOLAD3-I, Size of active core area  : ',I15,'.')
6002   FORMAT(T3,'@AOLAD3-I, Approximate number of passes : ',I10,'.')
6003   FORMAT(T3,'Beginning pass ',I10,' through AO integrals.') 
      ENDIF
1     CONTINUE
C
      J=1
      IOFF=I040
      DO 2 IRREP=1,NIRREP
       DO 3 I=1,IAOPOP(IRREP)
        ICORE(IOFF)=IRREP
        IOFF=IOFF+1
3      CONTINUE
2     CONTINUE
      CALL AOSYMVEC(ICORE(I050),NAO)
      CALL AOSYMVC2(ICORE(I051),ICORE(I040),NAO)
c      ITOPADR=MIN(ILOWADR+MXCOR-1,ITOTSIZ)
      IF(PRINT)WRITE(6,6003)NPASS
      CALL MPINFO2(ICORE(I090),ICORE(I080),ICORE(I070),
     &             IRRWRT,IDIS1,IDISSTRT,IRRSTRT,ICORE(I051),
     &             NAO,IRRBEGIN,IRREND,NSIZE,NWORDS)
      ITOPADR=ILOWADR+NWORDS-1
C
C LOAD UP AO INTEGRALS
C
      CALL GFNAME('IIII    ',FNAME,ILENGTH)
      OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),
     &     FORM='UNFORMATTED',STATUS='OLD')
      CALL ZERO(ICORE(I120),NWORDS+1)
      CALL AOLADMP(
     &              ICORE(I000),ICORE(I010),
     &              ICORE(I020),ICORE(I030),ICORE(I060),ICORE(I040),
     &              ICORE(I050),
     &              ICORE(I070),ICORE(I080),ICORE(I090),ICORE(I100),
     &              ICORE(I110),ICORE(I111),ICORE(I112),ICORE(I120),
     &              ILNBUF,LUINT,ISPIN,NAO,
     &              ISIZT2MIX,ILOWADR,IRRBEGIN,IRREND)
      CLOSE(UNIT=LUINT,STATUS='KEEP')
      IF(NIRREP.GT.1)THEN
       CALL GFNAME('IJIJ    ',FNAME,ILENGTH)
       OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),
     &      FORM='UNFORMATTED',STATUS='OLD')
       CALL AOLADMP(
     &              ICORE(I000),ICORE(I010),
     &              ICORE(I020),ICORE(I030),ICORE(I060),ICORE(I040),
     &              ICORE(I050),
     &              ICORE(I070),ICORE(I080),ICORE(I090),ICORE(I100),
     &              ICORE(I110),ICORE(I111),ICORE(I112),ICORE(I120),
     &              ILNBUF,LUINT,ISPIN,NAO,
     &              ISIZT2MIX,ILOWADR,IRRBEGIN,IRREND)
       CLOSE(UNIT=LUINT,STATUS='KEEP')
C
C SKIP NEXT TWO CALLS FOR IRREP=1 ONLY PROCESSING SINCE THESE INTEGRALS
C CANNOT CONTRIBUTE IN THESE CASES!
C
       IF(IRREND.GE.2)THEN
        CALL GFNAME('IIJJ    ',FNAME,ILENGTH)
        OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),
     &       FORM='UNFORMATTED',STATUS='OLD')
        CALL AOLADMP(
     &              ICORE(I000),ICORE(I010),
     &              ICORE(I020),ICORE(I030),ICORE(I060),ICORE(I040),
     &              ICORE(I050),
     &              ICORE(I070),ICORE(I080),ICORE(I090),ICORE(I100),
     &              ICORE(I110),ICORE(I111),ICORE(I112),ICORE(I120),
     &               ILNBUF,LUINT,ISPIN,NAO,
     &               ISIZT2MIX,ILOWADR,IRRBEGIN,IRREND)
        CLOSE(UNIT=LUINT,STATUS='KEEP')
        IF(NIRREP.GT.2.AND.IRREND.GE.2)THEN
         CALL GFNAME('IJKL    ',FNAME,ILENGTH)
         OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),
     &        FORM='UNFORMATTED',STATUS='OLD')
         CALL AOLADMP(
     &              ICORE(I000),ICORE(I010),
     &              ICORE(I020),ICORE(I030),ICORE(I060),ICORE(I040),
     &              ICORE(I050),
     &              ICORE(I070),ICORE(I080),ICORE(I090),ICORE(I100),
     &              ICORE(I110),ICORE(I111),ICORE(I112),ICORE(I120),
     &               ILNBUF,LUINT,ISPIN,NAO,
     &               ISIZT2MIX,ILOWADR,IRRBEGIN,IRREND)
         CLOSE(UNIT=LUINT,STATUS='KEEP')
        ENDIF
       ENDIF
      ENDIF
C
C FINISHED READING AO INTEGRALS.  CONTRACT THEM!
C
      I0Z=I000+ISIZT2MAX*IINTFP
      CALL LADCONT(ICORE(I000),ICORE(I0Z),ICORE(I120),
     &             ICORE(I090),ICORE(I080),ICORE(I100),IUHF,NPASS,
     &             IRREPX,LSTAOINC,LSTAO)
C
      IF(ITOPADR.NE.ITOTSIZ)THEN
       NPASS=NPASS+1
       ILOWADR=ITOPADR+1
       IDIS1=IDISSTRT
       IRRWRT=IRRSTRT
       GOTO 1
      ELSE
c       WRITE(6,1005)NPASS
c1005   FORMAT(T3,'@AOLADMP-I, AO integral processing required ',
c     &        I5,' passes.')
      ENDIF
C
      RETURN
      END
