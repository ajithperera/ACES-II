      PROGRAM FS
C
C PROGRAM FOR FOCK-SPACE COUPLED-CLUSTER CALCULATIONS.  WRITTEN
C AND DEVELOPED BY J.F. STANTON AND M. RITTBY.
C
      IMPLICIT INTEGER (A-Z)
      double precision eref,efinal
#include "maxbasfn.par"
      dimension ido(22),nactive(2),iactive(MAXBASFN,2)
      double precision x,snrm2,sdot
      COMMON // ICORE(1)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /ISTART/ I0,ICRSIZ
      COMMON /LISTS/ MOIO(10,500),MOIOWD(10,500),MOIOSZ(10,500),
     &               MOIODS(10,500),MOIOFL(10,500)
      COMMON /FSSYMPOP/ FSDPDAN(8,22),FSDPDNA(8,22),FSDPDAA(8,22),
     &                  FSDPDIN(8,22),FSDPDNI(8,22),FSDPDII(8,22),
     &                  FSDPDAI(8,22),FSDPDIA(8,22)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      common /refstate/ eref,efinal

      MXBAS = MAXBASFN 

      CALL CRAPSI(ICORE,IUHF,0)
      MAXCOR=ICRSIZ-IAND(ICRSIZ,1)
      call getrec(20,'JOBARC','TOTENERG',iintfp,eref)
c
      call izero(iactive,2*mxbas)

      do 10 i=1,22
       ido(i)=1
10    continue

      call getinput(nactive,iactive,mxbas,IUHF)

      nbas=nocco(1)+nvrto(1)
C
      CALL DUMPACT(ICORE(I0),NACTIVE,IACTIVE,NBAS,IUHF)
      
      CALL FSIODRV(ICORE,MAXCOR,IUHF,IDO,IACTIVE,NACTIVE,NBAS,I0)
C
C SOLVE 0,1 SECTOR
C
      Write(6,*) "The MAXCOR", MAXCOR
      CALL FS01 (ICORE(I0),MAXCOR,IUHF)

c      CALL FS02 (ICORE(I0),MAXCOR,IUHF)
c
      call putrec(20,'JOBARC','TOTENERG',iintfp,efinal)
      call aces_fin 
      stop
      end
