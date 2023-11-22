      
c**********************************************************
      SUBROUTINE RFIELD2(HCORE,P,DX,DY,DZ,DXT,DYT,DZT,BUF,IBUF,
     &                   PFULL,SCR,REPULS,DSIZE,NBAS,ILNBUF,
     &                   IUHF,LDIM2,ITER)
C-----------------------------------------------------------------------
C rfield2 calculates the reaction field contribution to the
C one-electron hamiltonian and to the energy.
c
c Piotr Rozyczko, fall '95
c
c Modification : JDW 3/11/96.
c   As noted by SB, the call to GETREC to get the nuclear
C   repulsion was made with incorrect length. The length
C   should be IINTFP, rather than 1. /MACHSP/ had to be
C   added, along with variable declarations.
c Modification : JDW 4/18/96.
c   I goofed up. Elements of /MACHSP/ were already declared (although
c   /MACHSP/ was not used and IBITWO was given instead of IBITWD) and
c   I unnecessarily added to declarations. Cray gave error messages
c   about multiple specifications. Now we try to get it right.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER DSIZE,ILENGTH,IERR,ISPIN,NBAS,N,L,K,IPRINT,
     $   INDEX,NUT,IINTLN,IFLTLN,IINTFP,IALONE,IBITWD,IUHF,
     $   IOFF,I,IBUF,IZERO,ILNBUF,IFLAGS,ITER,
     $   NPROTON,MXIRR2,J,ICNT,
     $   NNM102,IX,IEXTI,IEXTJ,INDI,INDJ,LDIM2
      INTEGER NIRREP,NBFIRR,IRPSZ1,IRPSZ2,IRPDS1,IRPDS2,
     &        IRPOFF,IREPS,DIRPRD,IWOFF1,IWOFF2,INEWVC,IDXVEC,
     &        ITRILN,ITRIOF,ISQRLN,ISQROF,ITHRU,ISTART,
     &        IRREP,IDIM,JDIM,IEND
      DOUBLE PRECISION PDX,PDY,PDZ,DIPNUX,DIPNUY,DIPNUZ,
     $   PX,PY,PZ,E1,E2,REPULS,FFACT,GFACT,EBORN,HCORE,CAVIT,
     $   BUF,P,DX,DY,DZ,REPOLD,PEX,PEY,PEZ,CHANGE,SCR,TRCPRD,
     $   PFULL,DXT,DYT,DZT
      CHARACTER*8 LABEL
      CHARACTER*32 JUNK
      CHARACTER*80 FNAME
      DIMENSION HCORE(DSIZE),P((IUHF+1)*DSIZE),SCR(LDIM2),
     &          PFULL(NBAS,NBAS)
      DIMENSION DX(NBAS,NBAS),DY(NBAS,NBAS),DZ(NBAS,NBAS),DXT(DSIZE),
     &          DYT(DSIZE),DZT(DSIZE)
      DIMENSION BUF(ILNBUF),IBUF(ILNBUF)
C
      DIMENSION IRPSZ1(36),IRPSZ2(28),DIRPRD(8,8)
      DIMENSION IRPDS1(36),IRPDS2(56)
      DIMENSION NBFIRR(8),IRPOFF(9),IREPS(9),INEWVC(1000)
      DIMENSION IWOFF1(37),IWOFF2(29)
      DIMENSION IDXVEC(1200),ITRILN(9),ITRIOF(8)
      DIMENSION ISQRLN(9),ISQROF(8)
C
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /SCRF /FFACT,GFACT,EBORN,NPROTON,E1,E2
C
      COMMON /SYMM2/NIRREP,NBFIRR,IRPSZ1,IRPSZ2,IRPDS1,IRPDS2,
     &              IRPOFF,IREPS,DIRPRD,IWOFF1,IWOFF2,INEWVC,IDXVEC,
     &              ITRILN,ITRIOF,ISQRLN,ISQROF,MXIRR2
      COMMON /MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      NNM102(IX)=(IX*(IX-1))/2
      IEXTI(IX)=1+(-1+INT(DSQRT(8.D0*IX+0.999D0)))/2
      IEXTJ(IX)=IX-NNM102(IEXTI(IX))
C
      IPRINT=IFLAGS(1)
      CALL GFNAME('VPOUT   ',FNAME,ILENGTH)
      OPEN(UNIT=30, FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED', 
     &     STATUS='OLD')
      REWIND(30)
C
c
c since we change the repulsion energy at each iter, get 
c the pristine repuls

      CALL ZERO(PFULL,NBAS*NBAS)
      CALL ZERO(SCR,MXIRR2)
    
      CALL GETREC(20,'JOBARC','NUCREP',IINTFP,REPULS)
C
c get nuclear dipoles (x,y,z) and dipole integrals Dx(m,n),
c Dy(m,n),Dz(m,n)
c
   3  LABEL='     X  '
C
      CALL SEEKLB(LABEL,IERR,0)
      BACKSPACE(30)
      READ(30)JUNK,DIPNUX
      DIPNUX=-DIPNUX
 197  READ(30)BUF,IBUF,NUT
      IF(NUT.EQ.-1)GOTO 4
      DO 17 INDEX=1,NUT
         INDI=IEXTI(IBUF(INDEX))
         INDJ=IEXTJ(IBUF(INDEX))
         DX(INDI,INDJ)=BUF(INDEX)
         DX(INDJ,INDI)=BUF(INDEX)
 17   CONTINUE
      GOTO 197
C
 4    LABEL='     Y  '
C
      CALL SEEKLB(LABEL,IERR,0)
      BACKSPACE(30)
      READ(30)JUNK,DIPNUY
      DIPNUY=-DIPNUY
 191  READ(30)BUF,IBUF,NUT
      IF(NUT.EQ.-1)GOTO 5
      DO 18 INDEX=1,NUT
         INDI=IEXTI(IBUF(INDEX))
         INDJ=IEXTJ(IBUF(INDEX))
         DY(INDI,INDJ)=BUF(INDEX)
         DY(INDJ,INDI)=BUF(INDEX)
 18   CONTINUE
      GOTO 191
c
 5    LABEL='     Z  '
C
      CALL SEEKLB(LABEL,IERR,0)
      BACKSPACE(30)
      READ(30)JUNK,DIPNUZ
      DIPNUZ=-DIPNUZ
 192  READ(30)BUF,IBUF,NUT
      IF(NUT.EQ.-1)GOTO 20
      DO 19 INDEX=1,NUT
         INDI=IEXTI(IBUF(INDEX))
         INDJ=IEXTJ(IBUF(INDEX))
         DZ(INDI,INDJ)=BUF(INDEX)
         DZ(INDJ,INDI)=BUF(INDEX)
 19   CONTINUE
      GOTO 192
C
 20   CONTINUE
      CLOSE(30)
C
      PDX=0.0D0
      PDY=0.0D0
      PDZ=0.0D0
c
c density matrix must be expanded to calculate Tr(P*D)
c    
      DO 100 ISPIN=1,(IUHF+1)
         CALL MKFULL(PFULL,P((ISPIN-1)*DSIZE+1),SCR,DSIZE,MXIRR2,NBAS,1)
CSSS         call expnd2(p((ISPIN-1)*DSIZE+1),pfull,nbas)
CSSS         IF (ITER .GT. 9) THEN
CSSS         Print*, "RFDEN MAT"
CSSS         CALL OUTPUT(PFULL, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
CSSS         Print*, "RFDX MAT"
CSSS         CALL OUTPUT(DX, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
CSSS         Print*, "RFDY MAT"
CSSD         CALL OUTPUT(DY, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
CSSS         Print*, "RFDZ MAT"
CSSS         CALL OUTPUT(DZ, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
CSSS         ENDIF
CSSS         call checksum("DEN-RF  ", pfull, nbas*nbas)
CSSS         call checksum("DX -RF  ", DX, nbas*nbas)
C
         PDX=PDX-TRCPRD(PFULL,DX,NBAS,NBAS)*FFACT
         PDY=PDY-TRCPRD(PFULL,DY,NBAS,NBAS)*FFACT
         PDZ=PDZ-TRCPRD(PFULL,DZ,NBAS,NBAS)*FFACT
 100  CONTINUE

C
      PX=DIPNUX+PDX
      PY=DIPNUY+PDY
      PZ=DIPNUZ+PDZ
Cdebug
      if (iprint.gt.10) then
         WRITE(6,*)'PDZ= Tr',PDX,PDY,PDZ
         WRITE(6,*)'PZ= Tr+DIPNUZ',PX,PY,PZ
      endif
C
c PEx,y,z- actual response to the reaction field 
c       
      PEX=PX*GFACT
      PEY=PY*GFACT
      PEZ=PZ*GFACT
cdebug
      if(iprint.gt.10) WRITE(6,*)'PEZ=(Tr+DIPNUZ)*GFACT',PEZ
C
c E1: nuclear contribution 
c E2: electronic contribution
c
      E1=(PEX*DIPNUX+PEY*DIPNUY+PEZ*DIPNUZ)
      E2=(PEX*PX+PEY*PY+PEZ*PZ)*0.5
C
      REPULS=REPULS-E1+E2+EBORN
Cdebug
      if (iprint.gt.10) then
         WRITE(6,*)'E1=',E1
         WRITE(6,*)'E2=',E2
         WRITE(6,*)'REPULS=',REPULS
CSSS         CALL OUTPUT(DZ, 1, NBAS, 1, NBAS, NBAS, NBAS, 1) 
      endif
c
c Now update the H matrix
c
      ITHRU  = 0
      ISTART = 1 
      DO IRREP=1, NIRREP
         IEND = ISTART + NBFIRR(IRREP) - 1
         DO IDIM = ISTART,  IEND
            DO JDIM = ISTART, IDIM
               ITHRU = ITHRU + 1
               DXT(ITHRU) = DX(IDIM, JDIM)
               DYT(ITHRU) = DY(IDIM, JDIM)
               DZT(ITHRU) = DZ(IDIM, JDIM)
            END DO
         END DO
         ISTART  = ISTART + NBFIRR(IRREP)
      END DO
  
      DO 30 I=1,DSIZE
         CHANGE=(PEX*DXT(I)+PEY*DYT(I)+PEZ*DZT(I))*1.0E-6
         HCORE(I)=HCORE(I)+CHANGE*FFACT
 30   CONTINUE
C
CSSS      Print*, "The core Hamiltonian"
CSSS      Write(*, "4(F10.5)") (HCORE(I), I=1, DSIZE)
      CALL ZERO(SCR,MXIRR2)
      CALL ZERO(DX,NBAS*NBAS)
      CALL ZERO(DY,NBAS*NBAS)
      CALL ZERO(DZ,NBAS*NBAS)
c
      Return
      END
