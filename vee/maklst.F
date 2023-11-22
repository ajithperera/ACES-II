
C CREATES LISTS USED IN EOM CALCULATIONS

#include "disttype.h"
#include "flags.h"

      SUBROUTINE MAKLST(ICORE,MAXCOR,IUHF)
      IMPLICIT INTEGER (A-Z)
      LOGICAL CIS,RPA,EOMCC,CISD,FULDIAG,INCORE,READGUES,ESPROP
      LOGICAL MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      LOGICAL NODAVID
      LOGICAL DAVID,BLOCK_DAVID,LANCZOS,G_TDENS 
c
CJDW
      LOGICAL TRIPNI,TRIPNI1,TRIPIT,T3STOR
      DIMENSION ICORE(MAXCOR),MAXOO(2),MAXVV(2),MAXVO(2)
      DIMENSION TYPEL(10),TYPER(10),NROOT(8)
      COMMON /METH/ CIS,RPA,EOMCC,CISD,FULDIAG,INCORE,READGUES
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /EXTRAP/ MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON /SYM/    POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /PROPGRAD/ ESPROP,IOPTROOT,IOPTSYM
      COMMON /PARTEOM/ NODAVID
      COMMON /REFTYPE/MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
CJDW 2/8/2008.
      COMMON /FLAGS/   IFLAGS(100)
      COMMON /FLAGS2/  IFLAGS2(500)
CJDW
C     Common blocks for triples i/o parameters.
C
      COMMON /TRIPLES/TRIPNI,TRIPNI1,TRIPIT,T3STOR
      COMMON /AUXIO / DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     &                NRECS,LUAUX
      COMMON /T3IOOF/ IJKPOS(8,8,8,2),IJKLEN(36,8,4),IJKOFF(36,8,4),
     &                NCOMB(4)
      COMMON /ROOTSEARCH/DAVID,BLOCK_DAVID,LANCZOS
C
      DATA TYPEL/13, 9,19,13,18,20,19, 9,19,20/
      DATA TYPER/11,20, 9,11,13,10,10,20, 9,10/
C
c      INEWFIL=1
      call aces_io_remove(54,'DERGAM')
      INEWFIL=0
      IZILCH = 0
C
C CALCULATE MAXIMUM SIZES OF VV, OO AND VO VECTORS
C
      MAXVV(1)=0
      MAXOO(1)=0
      MAXVO(1)=0
      DO IRREP = 1, NIRREP
         MAXVV(1)=MAX(MAXVV(1),IRPDPD(IRREP,19))
         MAXOO(1)=MAX(MAXOO(1),IRPDPD(IRREP,21))
         MAXVO(1)=MAX(MAXVO(1),IRPDPD(IRREP,9))
      END DO
      if (iuhf.eq.1) then
      MAXVV(2)=0
      MAXOO(2)=0
      MAXVO(2)=0
      DO IRREP = 1, NIRREP
         MAXVV(2)=MAX(MAXVV(2),IRPDPD(IRREP,20))
         MAXOO(2)=MAX(MAXOO(2),IRPDPD(IRREP,22))
         MAXVO(2)=MAX(MAXVO(2),IRPDPD(IRREP,10))
      END DO
      end if
C
C CREATE SINGLES VECTOR LISTS AND DIPOLE MOMENT LISTS
C
      DO 11 ISPIN=1,1+IUHF
C
CMN       IF(.NOT.CIS)THEN
        CALL UPDMOI(1,MAXVO(ISPIN),ISPIN,490,INEWFIL,0)
        INEWFIL=0
        CALL UPDMOI(1,MAXVO(ISPIN),ISPIN,493,INEWFIL,0)
        IF (LANCZOS) CALL UPDMOI(1,MAXVO(ISPIN),ISPIN+2,493,
     &                           INEWFIL,0)
        CALL UPDMOI(1,MAXVO(ISPIN),ISPIN+2,490,INEWFIL,0)
        CALL UPDMOI(1,MAXVO(ISPIN),9,447+ISPIN,INEWFIL,0)
        call aces_list_memset(ispin,  490,0)
        call aces_list_memset(ispin,  493,0)
        if (lanczos) call aces_list_memset(ispin+2,  493,0)
        call aces_list_memset(ispin+2,490,0)
        call aces_list_memset(9,      447+ISPIN,0)
CJDW 2/9/2008 Not really new. Need logic.
C
C        write(6,*) ' @MAKLST-I, setting up list 410 '
        CALL UPDMOI(1,MAXVO(ISPIN),ISPIN  ,410,INEWFIL,0)
        CALL UPDMOI(1,MAXVO(ISPIN),ISPIN+2,410,INEWFIL,0)
        CALL UPDMOI(1,MAXVO(ISPIN),ISPIN+4,410,INEWFIL,0)
C        write(6,*) ' @MAKLST-I, done with list 410 '
C
        IF(CISD.OR.CIS)THEN
c YAU : old
c         CALL INIPCK(1,ISYTYP(1,23),ISYTYP(2,23),40,IZILCH,IZILCH,1)
c         CALL INIPCK(1,ISYTYP(1,24),ISYTYP(2,24),41,IZILCH,IZILCH,1)
c         CALL INIPCK(1,ISYTYP(1,18),ISYTYP(2,18),42,IZILCH,IZILCH,1)
c YAU : new
          CALL INIPCK(1,h_DistType_AI,h_DistType_AI,40,IZILCH,0,1)
          CALL INIPCK(1,h_DistType_ai,h_DistType_ai,41,IZILCH,0,1)
          CALL INIPCK(1,h_DistType_AI,h_DistType_ai,42,IZILCH,0,1)
c YAU : end
        ENDIF
C
       DO 17 IRREP=1,NIRREP
        LENOO=IRPDPD(IRREP,20+ISPIN) 
        LENVV=IRPDPD(IRREP,18+ISPIN) 
        LENVO=IRPDPD(IRREP,8+ISPIN) 
        IF(EOMCC)THEN
         NUMRECS=6
        ELSE
         NUMRECS=3
        ENDIF
        CALL UPDMOI(NUMRECS,LENOO,IRREP,475+ISPIN,INEWFIL,0)
        INEWFIL=0
        CALL UPDMOI(NUMRECS,LENVV,IRREP,477+ISPIN,INEWFIL,0)
        CALL UPDMOI(NUMRECS,LENVO,IRREP,479+ISPIN,INEWFIL,0)
17     CONTINUE
11    CONTINUE
C
      IENTER=0
      IOFF=0
      IF(CIS)THEN
c       IENTER=1
c       IOFF=-1
       call aces_io_remove(51,'GAMLAM')
      ENDIF
      CALL UPDMOI(3,NFMI(1),1,160,IENTER,IOFF)
      IENTER=0
      IOFF=0
      CALL UPDMOI(3,NFEA(1),3,160,IENTER,IOFF)
      CALL UPDMOI(3,NT(1),5,160,IENTER,IOFF)
      IF(IUHF.NE.0)THEN
       CALL UPDMOI(3,NFMI(2),2,160,IENTER,IOFF)
       CALL UPDMOI(3,NFEA(2),4,160,IENTER,IOFF)
       CALL UPDMOI(3,NT(2),6,160,IENTER,IOFF)
      ENDIF
C
C CREATE AREA FOR Q(AB) AND Q(IJ) THREE-BODY INTERMEDIATES
C
      DO 110 ISPIN=1,1+IUHF
       MAXAB=0
       MAXIJ=0
       MAXIA=0
       DO 120 IRREP=1,NIRREP
        MAXAB=MAX(MAXAB,IRPDPD(IRREP,18+ISPIN))
        MAXIJ=MAX(MAXIJ,IRPDPD(IRREP,20+ISPIN))
        MAXIA=MAX(MAXIA,IRPDPD(IRREP,8+ISPIN))
120    CONTINUE
       CALL UPDMOI(1,MAXVO(ISPIN),ISPIN,490,INEWFIL,0)
       INEWFIL=0
       CALL UPDMOI(1,MAXAB,ISPIN,492,INEWFIL,0)
       CALL UPDMOI(1,MAXIJ,ISPIN,491,INEWFIL,0)
       CALL UPDMOI(1,MAXVO(ISPIN),ISPIN+2,490,INEWFIL,0)
       CALL UPDMOI(1,MAXAB,ISPIN+2,492,INEWFIL,0)
       CALL UPDMOI(1,MAXIJ,ISPIN+2,491,INEWFIL,0)
110   CONTINUE
C
C CALCULATE MAXIMUM SIZE OF DAVIDSON LISTS
C
C Check to see whether the genralized transition densities C (<L_k|p+q|R_l>) 
C are requested. In this case both R and L vecors for each state need to be stored 
C and the maximume allowed is 4 states per irrep. The calculation terminates in joda 
C if the user asked more than four (this amount maximum of 64 R and L
C vectors on the disk. Hope no one dare to try! Ajith Perera, Gunsan, South Korea, 10/2023
C 
      G_TDENS = (IFLAGS2(h_IFLAGS2_genrlizd_tdens) .Gt. 0) 

      If (G_TDENS) Then 
      Write(6,*)
      Write(6,"(2x,a,a)")"The geralized transition densities", 
     +            "(<L_k|p+q|R_l>) are requested. The R and L vectors"
      Write(6,"(2x,2a)") "for each state is stored (maximum of 4 per",
     +                   " irrep)."
      Endif 
    
      MAXLEN=0
      DO 25 IRREPX=1,NIRREP
       LEN=0
       DO 13 ISPIN=1,1+IUHF
        LEN=LEN+IRPDPD(IRREPX,8+ISPIN)
13     CONTINUE
       LEN=LEN+IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
       IF(IUHF.NE.0)THEN
        LEN=LEN+IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
        LEN=LEN+IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
       ENDIF
       MAXLEN=MAX(MAXLEN,LEN)
25    CONTINUE
      IF(ESPROP)THEN
       NUMLST=2
      ELSE
       NUMLST=1
      ENDIF
      IF (NODAVID) THEN
        MAXEXP0 = 0
        NREDUCE0 = 0
      ELSE
        MAXEXP0 = MAXEXP
        NREDUCE0 = NREDUCE
      ENDIF
      IF (G_TDENS) THEN 
           CALL GETREC(-1,'JOBARC','EESYMINF',NIRREP,NROOT)
           NTOT = 0
           DO IRREP= 1, NIRREP
              NTOT = NTOT + NROOT(IRREP)
           ENDDO
           IF (NTOT .GT. 32) THEN
              Write(6,"(2a)") "The maximum number of root requested",
     +                        " exceed the maximum (32) allowed for"
              Write(6,"(2a)") "generalized transition density",
     +                        "calculations."
              CALL ERREX
           ENDIF 
           WRITE(6,*) "NREDUCE0,NTOT",NREDUCE0,NTOT,NUMLST
      ENDIF 
      DO 27 I=1,NUMLST
       CALL UPDMOI(MAXEXP0,MAXLEN,I,470,0,0)
       CALL UPDMOI(MAXEXP0,MAXLEN,I,471,0,0)
       IF (G_TDENS) THEN
          CALL UPDMOI(6+NREDUCE0,MAXLEN,I,472,0,0)
          CALL UPDMOI(6+NREDUCE0+2*NTOT,MAXLEN,3,472,0,0)
       ELSE 
          CALL UPDMOI(6+NREDUCE0,MAXLEN,I,472,0,0)
       ENDIF 
27    CONTINUE

C
C NOW MAKE DENOMINATOR AND T2 LISTS   
C       
      DO 30 ISPIN=3,3-2*IUHF,-1
C
C DENOMINATOR AND T2 LISTS
C
      TTYPEL=ISYTYP(1,43+ISPIN)
      TTYPER=ISYTYP(2,43+ISPIN)
      IF ((RCCD .OR. DRCCD) .AND. ISPIN .NE.3) THEN
         TTYPEL= ISPIN
         TTYPER= 2 + ISPIN 
      ENDIF 
       CALL INIPCK2(1,TTYPEL,TTYPER,460+ISPIN,IZILCH,IZILCH,1)
       IF (LANCZOS) CALL INIPCK2(1,TTYPEL,TTYPER,466+ISPIN,IZILCH,
     &                           IZILCH,1)
       CALL INIPCK2(1,TTYPEL,TTYPER,447+ISPIN,IZILCH,IZILCH,1)
       CALL INIPCK2(1,TTYPEL,TTYPER,443+ISPIN,IZILCH,IZILCH,1)
C
CJDW 2/9/2008 [not new really]
C Make a set of lists for R2, L2, and increments for noniterative
C triples calculations. NOTE: need decision structure: we don't need
C these for iterative triples or CCSD.
       CALL INIPCK2(1,TTYPEL,TTYPER,400+ISPIN,IZILCH,IZILCH,1)
       CALL INIPCK2(1,TTYPEL,TTYPER,403+ISPIN,IZILCH,IZILCH,1)
       CALL INIPCK2(1,TTYPEL,TTYPER,406+ISPIN,IZILCH,IZILCH,1)
C
C FOR AO-BASED ALGORITHMS, NEED AO T2 LISTS
C
       IF(IFLAGS(93).EQ.2)THEN
        CALL INIPCK2(1,TTYPEL,TTYPER,280+ISPIN,IZILCH,IZILCH,1)
        TTYPEL=15
        CALL INIPCK2(1,TTYPEL,TTYPER,213+ISPIN,IZILCH,IZILCH,1)
c        CALL INIPCK2(1,TTYPEL,TTYPER,413+ISPIN,IZILCH,IZILCH,1)
        CALL INIPCK2(1,TTYPEL,TTYPER,463+ISPIN,IZILCH,IZILCH,1)
       ENDIF
C
30    CONTINUE
C
C MAKE RESORTED R2 AND L2 LISTS
C
      CALL INIPCK2(1,h_DistType_AI,h_DistType_ai,426,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_Ai,h_DistType_Ai,428,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_AI,h_DistType_AI,434,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_AI,h_DistType_ai,437,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_Ai,h_DistType_aI,439,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_AI,h_DistType_AI,440,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_ai,h_DistType_ai,441,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_AI,h_DistType_ai,442,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_Ai,h_DistType_aI,443,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_AI,h_DistType_AI,454,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_AI,h_DistType_ai,457,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_Ai,h_DistType_aI,459,IZILCH,0,1)
      do iGrp = 1, nirrep
         call aces_list_memset(iGrp,426,0)
         call aces_list_memset(iGrp,428,0)
         call aces_list_memset(iGrp,434,0)
         call aces_list_memset(iGrp,437,0)
         call aces_list_memset(iGrp,439,0)
         call aces_list_memset(iGrp,440,0)
         call aces_list_memset(iGrp,441,0)
         call aces_list_memset(iGrp,442,0)
         call aces_list_memset(iGrp,443,0)
         call aces_list_memset(iGrp,454,0)
         call aces_list_memset(iGrp,457,0)
         call aces_list_memset(iGrp,459,0)
      end do
      IF (IUHF.NE.0) THEN
      CALL INIPCK2(1,h_DistType_AI,h_DistType_AI,424,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_ai,h_DistType_ai,425,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_ai,h_DistType_AI,427,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_aI,h_DistType_aI,429,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_ai,h_DistType_ai,435,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_ai,h_DistType_AI,436,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_aI,h_DistType_Ai,438,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_ai,h_DistType_ai,455,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_ai,h_DistType_AI,456,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_aI,h_DistType_Ai,458,IZILCH,0,1)
      do iGrp = 1, nirrep
         call aces_list_memset(iGrp,424,0)
         call aces_list_memset(iGrp,425,0)
         call aces_list_memset(iGrp,427,0)
         call aces_list_memset(iGrp,429,0)
         call aces_list_memset(iGrp,435,0)
         call aces_list_memset(iGrp,436,0)
         call aces_list_memset(iGrp,438,0)
         call aces_list_memset(iGrp,455,0)
         call aces_list_memset(iGrp,456,0)
         call aces_list_memset(iGrp,458,0)
      end do
      END IF
C
      IF (IFLAGS(91).GT.1) THEN
      CALL INIPCK2(1,h_DistType_AltB,h_DistType_IltJ,114,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_altb,h_DistType_iltj,115,IZILCH,0,1)
      CALL INIPCK2(1,h_DistType_Ab,  h_DistType_Ij,  116,IZILCH,0,1)
      END IF
C
C-----------------------------------------------------------------------
C     Create lists for "differentiated" abci, mcjk quantities. These
C     handle the three-body terms in the TS and TD blocks. They are
C     needed in the CCSDT-3, CC3, and higher iterative methods, as
C     noniterative versions thereof.
C-----------------------------------------------------------------------
C
C      write(6,*) ' @MAKLST-I, IFLAGS2 24, 124',iflags2(24),iflags2(124)

      IF(IFLAGS(2).EQ.16.OR.IFLAGS(2).EQ.18.OR. IFLAGS(2).EQ.33.OR.
     &                                          IFLAGS(2).EQ.34.OR.
     &                     (IFLAGS(87).EQ.11.AND.IFLAGS2(124).GE.5))THEN
      IMODE = 0
C
      IF(IUHF.NE.0)THEN
        CALL INIPCK2(1, 3,16,377,imode,0,1)
        CALL INIPCK2(1, 4,17,378,imode,0,1)
        CALL INIPCK2(1,14,11,379,imode,0,1)
C
        CALL INIPCK2(1, 1, 9,381,IMODE,0,1)
        CALL INIPCK2(1, 2,10,382,IMODE,0,1)
        CALL INIPCK2(1,13,18,383,IMODE,0,1)
      ENDIF
C
        CALL INIPCK2(1,14,18,380,IMODE,0,1)
        CALL INIPCK2(1,13,11,384,IMODE,0,1)
C
      WRITE(6,*) ' @MAKLST-I, Creating mcjk lists 377-380 '
      WRITE(6,*) ' @MAKLST-I, Creating abci lists 381-384 '
C
      ENDIF
C-----------------------------------------------------------------------
C     Special scratch list
C-----------------------------------------------------------------------
      IF(
     &   (IFLAGS(2).EQ.22 .AND. IFLAGS2(124).GE.11 .AND. IUHF.NE.0) .OR.
     &   (IFLAGS(2).EQ.16 .AND. IUHF.NE.0)                         )THEN
C
      CALL INIPCK(1,IDOL,IDOR,411,IMODE,0,1)
C
      WRITE(6,*) ' @MAKLST-I, Creating scratch list 411 '
C
      ENDIF
C-----------------------------------------------------------------------
C
      SIZMAX=0
      DO 15 I=1,10
       SIZE=ISYMSZ(TYPEL(I),TYPER(I))
       IF(SIZE.GT.SIZMAX)THEN
        IDOL=TYPEL(I)
        IDOR=TYPER(I)
        SIZMAX=SIZE
       ENDIF
   15 CONTINUE

C
C      IF(NODAVID)THEN
C      DO ILIST=18,20
C        TTYPEL=ISYTYP(1,ILIST)
C        TTYPER=ISYTYP(2,ILIST)
C        CALL INIPCK2(1,TTYPEL,TTYPER,400+ILIST,IZILCH,IZILCH,1)
C      ENDDO
C      ENDIF
C
      RETURN
      END
