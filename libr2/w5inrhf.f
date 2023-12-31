      SUBROUTINE W5INRHF(ICORE,MAXCOR,IUHF,TERM1,TERM2,IOFFLIST)
C
C THIS ROUTINE CALCULATES THE INITIAL T2*W CONTRIBUTION TO THE
C W5 HBAR INTERMEDIATES FOR RHF CASES ONLY.
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL TERM1,TERM2,bRedundant
      DOUBLE PRECISION ONE,ONEM,TWO,ZILCH,HALF,snrm2
      CHARACTER*3 GETTYP
      DIMENSION ICORE(MAXCOR),ILOCT(8),ILOCW(8,8),NDIMW(8)
      DIMENSION POPDUM(8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFEA(2),NFMI(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /FLAGS2/IFLAGS2(500)
C
      DATA ONE,ONEM,ZILCH,TWO,HALF /1.0D0,-1.0D0,0.0D0,2.0D0,0.5D0/

      bRedundant = IFLAGS2(155).EQ.0
C
C BEGIN BY REORDERING <Ce|Am> => <Ae|Cm> ON DISK
C
      LISTW=30
      CALL SYMTRLST(ICORE,MAXCOR,IUHF,VRT,VRT,VRT,POP,
     &               1,19,19,19,9,9,30,30,1,.TRUE.,.FALSE.,.FALSE.)
      CALL IZERO(POPDUM,8)
      POPDUM(1)=1
C
C CALCULATE SIZES OF W(Ae,m) ARRAYS FOR EACH SYMMETRY OF C   
C AND OFFSETS INTO REORDERED W(em,A) ARRAY
C
      DO 5 IRREPC=1,NIRREP
       NDIMW(IRREPC)=0
       INCREM2=0
       DO 6 IRREPM=1,NIRREP
        IRREPAE=DIRPRD(IRREPM,IRREPC)
        INCREM=POP(IRREPM,1)*IRPDPD(IRREPAE,19)
        ILOCW(IRREPM,IRREPC)=IINTFP*INCREM2
        INCREM2=VRT(IRREPM,1)*IRPDPD(IRREPAE,9)+INCREM2
        NDIMW(IRREPC)=NDIMW(IRREPC)+INCREM
6      CONTINUE
5     CONTINUE
C
C READ IN T VECTOR AS T(em,bi) AND SPIN ADAPT IT 
C
      NSIZET=ISYMSZ(ISYTYP(1,37),ISYTYP(2,37))
      I000=1
      I010=I000+IINTFP*NSIZET
      IOFFT=1
      DO 7 IRREP=1,NIRREP
       ILOCT(IRREP)=IOFFT
       DISSYT=IRPDPD(IRREP,ISYTYP(1,37))
       NUMDST=DISSYT
       if(bRedundant) then 
          CALL GETLST(ICORE(IOFFT),1,NUMDST,1,IRREP,37)
          CALL GETLST(ICORE(I010),1,NUMDST,1,IRREP,39)
       else
          CALL GETLST_NR(ICORE(IOFFT),ICORE(I010),MAXCOR-I010+1,
     &                   37,IRREP)
          CALL GETLST(ICORE(I010),1,NUMDST,1,IRREP,239)
       endif
       CALL SSCAL (NUMDST*DISSYT,TWO,ICORE(IOFFT),1)
       CALL SAXPY (NUMDST*DISSYT,ONEM,ICORE(I010),1,ICORE(IOFFT),1)
       IOFFT=IOFFT+IINTFP*NUMDST*DISSYT 
7     CONTINUE

C FIRST TERM:
C
C   Z(Ab,Ci) <= [2 * T(im,eb) - T(mi,eb)] * <Ce|Am>
C
C INTEGRALS STORED NOW AS W(Ae,Cm)
C T VECTOR HELD IN CORE AS T(em,bi)
C
C LOOP OVER IRREPS OF TARGET
C
      DO 10 IRREPCI=1,NIRREP
       DO 11 IRREPI=1,NIRREP
        IRREPC=DIRPRD(IRREPI,IRREPCI)
        NUMI=POP(IRREPI,1)
        NUMC=VRT(IRREPC,1)
        I020=I010+IINTFP*NDIMW(IRREPC)
        I030=I020+IINTFP*NDIMW(IRREPC)
        DO 12 C=1,NUMC
C
C READ IN ALL W(Ae;m) FOR THIS VALUE OF c
C
         CALL GET3(ICORE(I020),LISTW,1,C,IRREPC,VRT,VRT,
     &             VRT,POP,19,9,'124',.FALSE.,.FALSE.,ICORE(I020))

C
C REORDER TO W(em;A)
C
         CALL SSTGEN(ICORE(I020),ICORE(I010),NDIMW(IRREPC),VRT,
     &               VRT,POP,POPDUM,ICORE(I030),IRREPC,'2314') 
C
         DO 13 I=1,NUMI
          DO 14 IRREPB=1,NIRREP
           IRREPBI=DIRPRD(IRREPB,IRREPI)
           IRREPEM=IRREPBI
           IRREPA=DIRPRD(IRREPB,IRREPCI)
           NUMB=VRT(IRREPB,1)
           NUMA=VRT(IRREPA,1)
           NROW=NUMA
           NCOL=NUMB
           NSUM=IRPDPD(IRREPEM,9)
           IOFFT=ILOCT(IRREPEM)+IINTFP*
     &           (NSUM*(ISYMOFF(IRREPI,IRREPEM,9)-1)+NSUM*NUMB*(I-1))
           iotmp=(NSUM*(ISYMOFF(IRREPI,IRREPEM,9)-1)+NSUM*NUMB*(I-1))
           IOFFW=I010+ILOCW(IRREPA,IRREPC)
           IOFFZ=I020+IINTFP*(ISYMOFF(IRREPB,IRREPCI,19)-1)
C
C FORM Z(ab) = W(em,a) * T(em,bi) FOR FIXED i
C
           CALL XGEMM('T','N',NROW,NCOL,NSUM,HALF,ICORE(IOFFW),NSUM,
     &                ICORE(IOFFT),NSUM,ZILCH,ICORE(IOFFZ),NUMA)
14        CONTINUE
C
C CALCULATE ADDRESS OF CI RECORD
C          
          IREC=C+NUMC*(I-1)+ISYMOFF(IRREPI,IRREPCI,9)-1
          CALL PUTLST(ICORE(I020),IREC,1,1,IRREPCI,130+IOFFLIST)
13       CONTINUE
12      CONTINUE
11     CONTINUE
10    CONTINUE 
C
C END OF FIRST CONTRACTION
C
C REORDER INTEGRALS ON DISK AGAIN W(Ae,Cm)=>W(Ce,Am) AND UNSPIN-ADAPT
C THEM
C
      LISTW=30
      CALL SYMTRLST(ICORE,MAXCOR,IUHF,VRT,VRT,VRT,POP,
     &               1,19,19,19,9,9,30,30,1,.FALSE.,.FALSE.,.FALSE.)
c
c temporary code to unspinadapt integrals
c
      DO 1000 IRREPCI=1,NIRREP
       DISSYW=IRPDPD(IRREPCI,ISYTYP(1,30))
        NUMDSW=IRPDPD(IRREPCI,ISYTYP(2,30))
       DO 1001 IDIS=1,NUMDSW
        CALL GETLST(ICORE(I000),IDIS,1,1,IRREPCI,30)
        ITMP=I000+IINTFP*DISSYW
        CALL SYMTRA(IRREPCI,VRT,VRT,1,ICORE(I000),ICORE(ITMP))
        CALL SAXPY (DISSYW,TWO,ICORE(I000),1,ICORE(ITMP),1)
        CALL SSCAL (DISSYW,1.0D0/3.0D0,ICORE(ITMP),1)
        CALL PUTLST(ICORE(ITMP),IDIS,1,1,IRREPCI,30)
1001   CONTINUE
1000  CONTINUE
C
C NOW REORDER INTEGRALS W(Ae,Cm) => W(AC,em)
C      
      CALL SYMTRLST(ICORE,MAXCOR,IUHF,VRT,VRT,VRT,POP,
     &               1,19,19,19,9,9,30,30,2,.FALSE.,.FALSE.,.FALSE.)
C   
C SECOND CONTRACTION: T(Im,Ea)*<Ec|Bm>     
C
C READ IN T VECTOR AS T(Em,aI) 
C
      NSIZET=ISYMSZ(ISYTYP(1,39),ISYTYP(2,39))
      I000=1
      I010=I000+IINTFP*NSIZET
      IF(bRedundant) THEN
         CALL GETALL(ICORE(I000),NSIZET,1,39)
      ELSE
         CALL GETALL(ICORE(I000),NSIZET,1,239)
      ENDIF

C      
C
C   Z(Ab,Ci) <= T(Em,aI) * <Ec|Bm>
C
C INTEGRALS STORED NOW AS W(EB,cm)
C T VECTOR HELD IN CORE AS T(Em,aI)
C
C LOOP OVER IRREPS OF TARGET
C
      DO 110 IRREPCI=1,NIRREP
       DO 111 IRREPI=1,NIRREP
        IRREPC=DIRPRD(IRREPI,IRREPCI)
        NUMI=POP(IRREPI,1)
        NUMC=VRT(IRREPC,1)
        I020=I010+IINTFP*NDIMW(IRREPC)
        I030=I020+IINTFP*NDIMW(IRREPC)
        DO 112 C=1,NUMC
C
C READ IN ALL W(EB;m) FOR THIS VALUE OF c
C
         CALL GET3(ICORE(I020),LISTW,1,C,IRREPC,VRT,VRT,
     &             VRT,POP,19,9,'124',.FALSE.,.FALSE.,ICORE(I020))
C
C REORDER TO W(Em;B)
C
         CALL SSTGEN(ICORE(I020),ICORE(I010),NDIMW(IRREPC),VRT,
     &               VRT,POP,POPDUM,ICORE(I030),IRREPC,'1324')
C
         DO 113 I=1,NUMI
          LENAB=IRPDPD(IRREPCI,19)
          DO 114 IRREPB=1,NIRREP
           IRREPA=DIRPRD(IRREPB,IRREPCI)
           IRREPBI=DIRPRD(IRREPB,IRREPI)
           IRREPAI=DIRPRD(IRREPA,IRREPI)
           IRREPEM=IRREPAI
           NUMB=VRT(IRREPB,1)
           NUMA=VRT(IRREPA,1)
           NROW=NUMB
           NCOL=NUMA
           NSUM=IRPDPD(IRREPEM,9)
           IOFFT=ILOCT(IRREPEM)+IINTFP*
     &           (NSUM*(ISYMOFF(IRREPI,IRREPEM,9)-1)+NSUM*NUMA*(I-1))
           IOFFW=I010+ILOCW(IRREPB,IRREPC)
           IOFFZ=I020+IINTFP*(ISYMOFF(IRREPA,IRREPCI,19)-1)
C
C FORM Z(ba) = W(em,b) * T(em,ai) FOR FIXED i
C
           CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFW),NSUM,
     &                ICORE(IOFFT),NSUM,ZILCH,ICORE(IOFFZ),NUMB)
114       CONTINUE
C
C CALCULATE ADDRESS OF CI RECORD AND FORM 1/2 Z(ba)+Z(ab)
C          
          ITMP=I020+IINTFP*LENAB
          CALL SYMTRA(IRREPCI,VRT,VRT,1,ICORE(I020),ICORE(ITMP))
          CALL SAXPY (LENAB,HALF,ICORE(I020),1,ICORE(ITMP),1)
          IREC=C+NUMC*(I-1)+ISYMOFF(IRREPI,IRREPCI,9)-1
          CALL GETLST(ICORE(I020),IREC,1,1,IRREPCI,130+IOFFLIST)
          CALL SAXPY (LENAB,ONEM,ICORE(I020),1,ICORE(ITMP),1)
          CALL PUTLST(ICORE(ITMP),IREC,1,1,IRREPCI,130+IOFFLIST)
113      CONTINUE
112     CONTINUE
111    CONTINUE
c
c for debugging ...
c
c       numdis=irpdpd(irrepci,isytyp(2,30))
c       dissiz=irpdpd(irrepci,isytyp(1,30))
c       call getlst(icore(i010),1,numdis,1,irrepci,130+iofflst)
c       write(6,*)' checksum for irrep ',irrepci,' is ',
c     &            snrm2(numdis*dissiz,icore(i010),1)
110   CONTINUE 
C
C NOW REORDER INTEGRALS W(AC,em) => W(Ae,Cm)
C      
      CALL SYMTRLST(ICORE,MAXCOR,IUHF,VRT,VRT,VRT,POP,
     &               1,19,19,19,9,9,30,30,2,.FALSE.,.FALSE.,
     &               .FALSE.)
C
C NOW ADD LISTS 30 AND 130
C
      DO 2000 IRREP=1,NIRREP
       I000=1
       DISSIZ=IRPDPD(IRREP,ISYTYP(1,30))
       NUMDIS=IRPDPD(IRREP,ISYTYP(2,30))
       IF(DISSIZ.NE.0)THEN
        NINCOR=MAXCOR/(DISSIZ*IINTFP)
       ELSE
        NINCOR=2*NUMDIS
       ENDIF
       NBUNCH=MIN(NINCOR/2,NUMDIS)
       I010=I000+IINTFP*DISSIZ*NBUNCH
       ISTART=1
       NLEFT=NUMDIS
1      NREAD=MIN(NLEFT,NBUNCH)
       IF(TERM1)THEN
        CALL GETLST(ICORE(I000),ISTART,NREAD,1,IRREP,30)
       ELSE
        CALL ZERO  (ICORE(I000),NREAD*DISSIZ)
       ENDIF
       IF(TERM2)THEN
        CALL GETLST(ICORE(I010),ISTART,NREAD,1,IRREP,130+IOFFLIST)
       ELSE
        CALL ZERO  (ICORE(I010),NREAD*DISSIZ)
       ENDIF
       CALL SAXPY (NREAD*DISSIZ,ONEM,ICORE(I010),1,ICORE(I000),1)
       CALL PUTLST(ICORE(I000),ISTART,NREAD,1,IRREP,130+IOFFLIST)
       NLEFT=NLEFT-NREAD
       ISTART=ISTART+NREAD
       IF(NLEFT.NE.0)GOTO 1
2000  CONTINUE
C
      RETURN
      END
