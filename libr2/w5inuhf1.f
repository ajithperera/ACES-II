      SUBROUTINE W5INUHF1(ICORE,MAXCOR,IUHF,TERM1,TERM2,IOFFLIST)
C
C THIS ROUTINE CALCULATES THE INITIAL T2*W CONTRIBUTION TO THE
C W5 HBAR INTERMEDIATES FOR AAAA AND BBBB SPIN CASES.
C
C W(A<B,CI) =  T(IM,BE) * <AM||EC> + T(Im,Be) * <Am|Ce> [AAAA]
C
C W(a<b,ci) =  T(im,be) * <am||ec> + T(iM,bE) * <aM|cE> [BBBB]
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL TERM1,TERM2
      logical bRedundant
      CHARACTER*4 SPCASE
      DOUBLE PRECISION ONE,ONEM,TWO,ZILCH,HALF,snrm2
      CHARACTER*3 GETTYP
      DIMENSION ICORE(MAXCOR),ILOCT(8),ILOCW(8,8),NDIMW(8)
      DIMENSION POPDUM(8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFEA(2),NFMI(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /FLAGS2/IFLAGS2(500)
C
      DATA ONE,ONEM,ZILCH,TWO,HALF /1.0D0,-1.0D0,0.0D0,2.0D0,0.5D0/
      bRedundant = iflags2(155).eq.0

C
C LOOP OVER SPIN CASES
C
      DO 2 ISPIN=1,2
       CALL NEWTYP2(1,401,18+ISPIN,8+ISPIN,.TRUE.)
C
C FIRST DO AAAA->AAAA AND BBBB->BBBB CONTRACTIONS
C
C
C BEGIN BY REORDERING W(C<E,AM) => W(AE,CM) ON DISK
C
       LISTW=26+ISPIN
       LISTZ=126+ISPIN+IOFFLIST
       LSTSCR=401
       CALL SYMTRLST(ICORE,MAXCOR,IUHF,
     &               VRT(1,ISPIN),VRT(1,ISPIN),VRT(1,ISPIN),
     &               POP(1,ISPIN),1,ISPIN,18+ISPIN,18+ISPIN,
     &               8+ISPIN,8+ISPIN,LISTW,LSTSCR,1,
     &              .FALSE.,.TRUE.,.FALSE.)
       CALL IZERO(POPDUM,8)
       POPDUM(1)=1
C
C CALCULATE SIZES OF W(AE,M) ARRAYS FOR EACH SYMMETRY OF C   
C AND OFFSETS INTO REORDERED W(EM,A) ARRAY
C
       DO 5 IRREPC=1,NIRREP
        NDIMW(IRREPC)=0
        INCREM2=0
        DO 6 IRREPM=1,NIRREP
         IRREPAE=DIRPRD(IRREPM,IRREPC)
         INCREM=POP(IRREPM,ISPIN)*IRPDPD(IRREPAE,18+ISPIN)
         ILOCW(IRREPM,IRREPC)=IINTFP*INCREM2
         INCREM2=VRT(IRREPM,ISPIN)*IRPDPD(IRREPAE,8+ISPIN)+INCREM2
         NDIMW(IRREPC)=NDIMW(IRREPC)+INCREM
6       CONTINUE
5      CONTINUE
C
C READ IN T VECTOR AS T(EM,BI)
C
       LISTT =33+ISPIN
       NSIZET=ISYMSZ(ISYTYP(1,33+ISPIN),ISYTYP(2,33+ISPIN))
       I000=1
       I010=I000+IINTFP*NSIZET
       if(bRedundant) then
          CALL GETALL(ICORE(I000),NSIZET,1,LISTT)
       else
          CALL GETALL(ICORE(I000),NSIZET,1,LISTT+200)
       endif
       IOFFT=1
       DO 7 IRREP=1,NIRREP
        ILOCT(IRREP)=IOFFT
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
        NUMDST=DISSYT
        IOFFT=IOFFT+IINTFP*NUMDST*DISSYT 
7      CONTINUE
C      
C FIRST TERM:
C
C   Z(AB,CI) <= T(IM,EB)*<CE||AM>
C
C INTEGRALS STORED NOW AS W(AE,CM)
C T VECTOR HELD IN CORE AS T(EM,BI)
C
C LOOP OVER IRREPS OF TARGET
C
       DO 10 IRREPCI=1,NIRREP
        TARSIZ=IRPDPD(IRREPCI,ISPIN) 
        DO 11 IRREPI=1,NIRREP
         IRREPC=DIRPRD(IRREPI,IRREPCI)
         NUMI=POP(IRREPI,ISPIN)
         NUMC=VRT(IRREPC,ISPIN)
         I020=I010+IINTFP*NDIMW(IRREPC)
         I030=I020+IINTFP*NDIMW(IRREPC)
         DO 12 C=1,NUMC
C
C READ IN ALL W(AE;M) FOR THIS VALUE OF C
C
          CALL GET3(ICORE(I020),LSTSCR,1,C,IRREPC,
     &              VRT(1,ISPIN),VRT(1,ISPIN),VRT(1,ISPIN),
     &              POP(1,ISPIN),18+ISPIN,8+ISPIN,'124',
     &              .FALSE.,.FALSE.,ICORE(I020))
C
C REORDER TO W(EM;A)
C
          CALL SSTGEN(ICORE(I020),ICORE(I010),NDIMW(IRREPC),
     &                VRT(1,ISPIN),VRT(1,ISPIN),POP(1,ISPIN),
     &                POPDUM,ICORE(I030),IRREPC,'2314')
C
          DO 13 I=1,NUMI
           DO 14 IRREPB=1,NIRREP
            IRREPBI=DIRPRD(IRREPB,IRREPI)
            IRREPEM=IRREPBI
            IRREPA=DIRPRD(IRREPB,IRREPCI)
            NUMB=VRT(IRREPB,ISPIN)
            NUMA=VRT(IRREPA,ISPIN)
            NROW=NUMA
            NCOL=NUMB
            NSUM=IRPDPD(IRREPEM,8+ISPIN)
            IOFFT=ILOCT(IRREPEM)+IINTFP*
     &            (NSUM*(ISYMOFF(IRREPI,IRREPEM,8+ISPIN)-1)+
     &            NSUM*NUMB*(I-1))
            IOFFW=I010+ILOCW(IRREPA,IRREPC)
            IOFFZ=I020+IINTFP*(ISYMOFF(IRREPB,IRREPCI,18+ISPIN)-1)
C
C FORM Z(AB) = W(EM,A) * T(EM,BI) FOR FIXED i
C
            If (nsum .ne. 0) 
     &      CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFW),NSUM,
     &                 ICORE(IOFFT),NSUM,ZILCH,ICORE(IOFFZ),NUMA)
14         CONTINUE
C
C CALCULATE ADDRESS OF CI RECORD
C          
           IREC=C+NUMC*(I-1)+ISYMOFF(IRREPI,IRREPCI,8+ISPIN)-1
           CALL ASSYM2(IRREPCI,VRT(1,ISPIN),1,ICORE(I020))
           CALL PUTLST(ICORE(I020),IREC,1,1,IRREPCI,LISTZ)
13        CONTINUE
12       CONTINUE
11      CONTINUE
10     CONTINUE
C
C NOW DO THE ABAB AND BABA CONTRACTIONS
C
       IF(ISPIN.EQ.1)THEN
        LISTW=30
        ISPINE=2
        ISPINM=2 
        CALL NEWTYP2(1,401,13,11,.TRUE.)
C
C BEGIN BY REORDERING W(Ce,Am) => W(Ae,Cm)
C
        CALL SYMTRLST(ICORE,MAXCOR,IUHF,
     &                VRT(1,1),VRT(1,2),VRT(1,1),POP(1,2),
     &                1,13,13,13,11,11,LISTW,LSTSCR,1,
     &                .FALSE.,.FALSE.,.FALSE.)
        CALL IZERO(POPDUM,8)
        POPDUM(1)=1
C
C CALCULATE SIZES OF W(Ae,m) ARRAYS FOR EACH SYMMETRY OF C   
C AND OFFSETS INTO REORDERED W(em,A) ARRAY
C
        DO 105 IRREPC=1,NIRREP
         NDIMW(IRREPC)=0
         INCREM2=0
         DO 106 IRREPM=1,NIRREP
          IRREPAE=DIRPRD(IRREPM,IRREPC)
          INCREM=POP(IRREPM,ISPINM)*IRPDPD(IRREPAE,13)
          ILOCW(IRREPM,IRREPC)=IINTFP*INCREM2
          INCREM2=VRT(IRREPM,ISPIN)*IRPDPD(IRREPAE,10)+INCREM2
          NDIMW(IRREPC)=NDIMW(IRREPC)+INCREM
106      CONTINUE
105     CONTINUE
C
C READ IN T VECTOR AS T(em,BI)
C
        LISTT =36
        NSIZET=ISYMSZ(ISYTYP(1,36),ISYTYP(2,36))
        I000=1
        I010=I000+IINTFP*NSIZET
        if(bRedundant) then
           CALL GETALL(ICORE(I000),NSIZET,1,LISTT)
        else
           CALL GETALL_NR(ICORE(I000),ICORE(I010),MAXCOR-I010+1,LISTT)
c  endif bRedundant
        endif
        IOFFT=1
        DO 107 IRREP=1,NIRREP
         ILOCT(IRREP)=IOFFT
         DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
         NUMDST=IRPDPD(IRREP,ISYTYP(2,LISTT))
         IOFFT=IOFFT+IINTFP*NUMDST*DISSYT 
107     CONTINUE
C     
C SECOND TERM:
C
C   Z(AB,CI) <= T(Im,Be)*<Ce|Am>
C
C INTEGRALS STORED NOW AS W(Ae,Cm)
C T VECTOR HELD IN CORE AS T(em,BI)
C
C LOOP OVER IRREPS OF TARGET
C
        DO 110 IRREPCI=1,NIRREP
         TARSIZ=IRPDPD(IRREPCI,ISPIN) 
         DO 111 IRREPI=1,NIRREP
          IRREPC=DIRPRD(IRREPI,IRREPCI)
          NUMI=POP(IRREPI,ISPIN)
          NUMC=VRT(IRREPC,ISPIN)
          I020=I010+IINTFP*NDIMW(IRREPC)
          I030=I020+IINTFP*NDIMW(IRREPC)
          DO 112 C=1,NUMC
C
C READ IN ALL W(Ae;m) FOR THIS VALUE OF C
C
           CALL GET3(ICORE(I020),LSTSCR,1,C,IRREPC,
     &               VRT(1,1),VRT(1,2),VRT(1,1),
     &               POP(1,2),13,11,'124',
     &               .FALSE.,.FALSE.,ICORE(I020))
C
C REORDER TO W(em;A)
C
           CALL SSTGEN(ICORE(I020),ICORE(I010),NDIMW(IRREPC),
     &                 VRT(1,1),VRT(1,2),POP(1,2),
     &                 POPDUM,ICORE(I030),IRREPC,'2314')
C
           DO 113 I=1,NUMI
            DO 114 IRREPB=1,NIRREP
             IRREPBI=DIRPRD(IRREPB,IRREPI)
             IRREPEM=IRREPBI
             IRREPA=DIRPRD(IRREPB,IRREPCI)
             NUMB=VRT(IRREPB,ISPIN)
             NUMA=VRT(IRREPA,ISPIN)
             NROW=NUMA
             NCOL=NUMB
             NSUM=IRPDPD(IRREPEM,10)
             IOFFT=ILOCT(IRREPEM)+IINTFP*
     &             (NSUM*(ISYMOFF(IRREPI,IRREPEM,9)-1)+
     &             NSUM*NUMB*(I-1))
             IOFFW=I010+ILOCW(IRREPA,IRREPC)
             IOFFZ=I020+IINTFP*(ISYMOFF(IRREPB,IRREPCI,18+ISPIN)-1)
C
C FORM Z(AB) = W(em,A) * T(em,BI) FOR FIXED i
C
             If (nsum .ne. 0) 
     &       CALL XGEMM('T','N',NROW,NCOL,NSUM,ONEM,ICORE(IOFFW),NSUM,
     &                  ICORE(IOFFT),NSUM,ZILCH,ICORE(IOFFZ),NUMA)
114         CONTINUE
C
C CALCULATE ADDRESS OF CI RECORD
C          
            IREC=C+NUMC*(I-1)+ISYMOFF(IRREPI,IRREPCI,8+ISPIN)-1
            CALL ASSYM2(IRREPCI,VRT(1,ISPIN),1,ICORE(I020))
            CALL GETLST(ICORE(I030),IREC,1,1,IRREPCI,LISTZ)
            CALL SAXPY (IRPDPD(IRREPCI,ISPIN),ONE,ICORE(I030),1,
     &                  ICORE(I020),1)
            CALL PUTLST(ICORE(I020),IREC,1,1,IRREPCI,LISTZ)
113        CONTINUE
112       CONTINUE
111      CONTINUE
110     CONTINUE
C
       ELSEIF(ISPIN.EQ.2)THEN
        LISTW=29
        ISPINE=1
        ISPINM=1 
        CALL NEWTYP2(1,401,9,20,.TRUE.)
C
C BEGIN BY REORDERING W(Ec,Ma) => W(EM,ca)
C
        CALL SYMTRLST(ICORE,MAXCOR,IUHF,
     &                VRT(1,1),VRT(1,2),POP(1,1),VRT(1,2),
     &                1,13,9,9,18,20,LISTW,LSTSCR,2,
     &                .FALSE.,.FALSE.,.FALSE.)
        CALL IZERO(POPDUM,8)
        POPDUM(1)=1
C
C CALCULATE SIZES OF W(EM,a) ARRAYS FOR EACH SYMMETRY OF C   
C AND OFFSETS INTO REORDERED W(EM,a) ARRAY
C
        DO 205 IRREPC=1,NIRREP
         NDIMW(IRREPC)=0
         INCREM2=0
         DO 206 IRREPM=1,NIRREP
          IRREPAE=DIRPRD(IRREPM,IRREPC)
          INCREM=VRT(IRREPM,ISPIN)*IRPDPD(IRREPAE,9)
          ILOCW(IRREPM,IRREPC)=IINTFP*INCREM2
          INCREM2=VRT(IRREPM,ISPIN)*IRPDPD(IRREPAE,9)+INCREM2
          NDIMW(IRREPC)=NDIMW(IRREPC)+INCREM
206      CONTINUE
205     CONTINUE
C
C READ IN T VECTOR AS T(EM,bi)
C
        LISTT =37
        NSIZET=ISYMSZ(ISYTYP(1,37),ISYTYP(2,37))
        I000=1
        I010=I000+IINTFP*NSIZET
        if(.not.bRedundant) then
           CALL GETALL_NR(ICORE(I000),ICORE(I010),MAXCOR-I010+1,LISTT)
        else
           CALL GETALL(ICORE(I000),NSIZET,1,LISTT)
c endif bRedundant
        endif
        IOFFT=1
        DO 207 IRREP=1,NIRREP
         ILOCT(IRREP)=IOFFT
         DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
         NUMDST=IRPDPD(IRREP,ISYTYP(2,LISTT))
         IOFFT=IOFFT+IINTFP*NUMDST*DISSYT 
207     CONTINUE
C     
C SECOND TERM:
C
C   Z(AB,CI) <= T(Mi,Eb)*<Ec|Ma>
C
C INTEGRALS STORED NOW AS W(EM,ca)
C T VECTOR HELD IN CORE AS T(EM,bi)
C
C LOOP OVER IRREPS OF TARGET
C
        DO 210 IRREPCI=1,NIRREP
         TARSIZ=IRPDPD(IRREPCI,ISPIN) 
         DO 211 IRREPI=1,NIRREP
          IRREPC=DIRPRD(IRREPI,IRREPCI)
          NUMI=POP(IRREPI,ISPIN)
          NUMC=VRT(IRREPC,ISPIN)
          I020=I010+IINTFP*NDIMW(IRREPC)
          I030=I020+IINTFP*MAX(TARSIZ,NDIMW(IRREPC))
          DO 212 C=1,NUMC
C
C READ IN ALL W(EM;a) FOR THIS VALUE OF C
C
           CALL GET3(ICORE(I010),LSTSCR,1,C,IRREPC,
     &               VRT(1,1),POP(1,1),VRT(1,2),
     &               VRT(1,2),9,20,'124',
     &               .FALSE.,.FALSE.,ICORE(I010))
C
           DO 213 I=1,NUMI
            DO 214 IRREPB=1,NIRREP
             IRREPBI=DIRPRD(IRREPB,IRREPI)
             IRREPEM=IRREPBI
             IRREPA=DIRPRD(IRREPB,IRREPCI)
             NUMB=VRT(IRREPB,ISPIN)
             NUMA=VRT(IRREPA,ISPIN)
             NROW=NUMA
             NCOL=NUMB
             NSUM=IRPDPD(IRREPEM,9)
             IOFFT=ILOCT(IRREPEM)+IINTFP*
     &             (NSUM*(ISYMOFF(IRREPI,IRREPEM,10)-1)+
     &             NSUM*NUMB*(I-1))
             IOFFW=I010+ILOCW(IRREPA,IRREPC)
             IOFFZ=I020+IINTFP*(ISYMOFF(IRREPB,IRREPCI,18+ISPIN)-1)
C
C FORM Z(AB) = W(EM,a) * T(EM,bi) FOR FIXED i
C
             If (nsum .ne. 0)
     &       CALL XGEMM('T','N',NROW,NCOL,NSUM,ONEM,ICORE(IOFFW),NSUM,
     &                  ICORE(IOFFT),NSUM,ZILCH,ICORE(IOFFZ),NUMA)
214         CONTINUE
C
C CALCULATE ADDRESS OF CI RECORD
C          
            IREC=C+NUMC*(I-1)+ISYMOFF(IRREPI,IRREPCI,8+ISPIN)-1
            CALL ASSYM2(IRREPCI,VRT(1,ISPIN),1,ICORE(I020))
            CALL GETLST(ICORE(I030),IREC,1,1,IRREPCI,LISTZ)
            CALL SAXPY (IRPDPD(IRREPCI,ISPIN),ONE,ICORE(I030),1,
     &                  ICORE(I020),1)
            CALL PUTLST(ICORE(I020),IREC,1,1,IRREPCI,LISTZ)
213        CONTINUE
212       CONTINUE
211      CONTINUE
210     CONTINUE
C
       ENDIF
C
C NOW ADD LISTS 26+ISPIN AND 126+ISPIN
C
       DO 2000 IRREP=1,NIRREP
        LISTW=26+ISPIN
        LISTZ=126+ISPIN+IOFFLIST
        I000=1
        DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
        NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
        IF(DISSIZ.NE.0)THEN
         NINCOR=MAXCOR/(DISSIZ*IINTFP)
        ELSE
         NINCOR=2*NUMDIS
        ENDIF
        NBUNCH=MIN(NINCOR/2,NUMDIS)
        I010=I000+IINTFP*DISSIZ*NBUNCH
        ISTART=1
        NLEFT=NUMDIS
1       NREAD=MIN(NLEFT,NBUNCH)
        IF(TERM1)THEN
         CALL GETLST(ICORE(I000),ISTART,NREAD,1,IRREP,LISTW)
        ELSE
         CALL ZERO  (ICORE(I000),NREAD*DISSIZ)
        ENDIF
        IF(TERM2)THEN
         CALL GETLST(ICORE(I010),ISTART,NREAD,1,IRREP,LISTZ)
        ELSE
         CALL ZERO  (ICORE(I010),NREAD*DISSIZ)
        ENDIF
        CALL SAXPY (NREAD*DISSIZ,ONEM,ICORE(I010),1,ICORE(I000),1)
        CALL PUTLST(ICORE(I000),ISTART,NREAD,1,IRREP,LISTZ)
        NLEFT=NLEFT-NREAD
        ISTART=ISTART+NREAD
        IF(NLEFT.NE.0)GOTO 1
2000   CONTINUE
2     CONTINUE 
      RETURN
      END
