
      SUBROUTINE OOFCOR(ICORE,MAXMEM,IRPPMEM,ILNBUF,NAOBK,MXLSSZ1,
     &                  MXLSSZ2,TOL,REPULS,NBAS,IUHF,YESNO2,
     &                  DAMP,DMPFLG,DMPTOL,DEAVG,SEWARD,IORDGRP,NSIZE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER(LUINT=10)
      PARAMETER(LUHF2=15)
C
      INTEGER SCFPRNT
      LOGICAL NOCONV,YESNO2,CHGOCC
      LOGICAL GSSOPT,GSSALW,GSSALT,GSSLOK,GSSRED,GSSWRT,GSSUFR
      LOGICAL DMPFLG,SEWARD,POSTSCF
      LOGICAL ROHFMO
      LOGICAL CHKPHASE, FIRSTSCF, YESNO3
      DIMENSION IOLDOCC(16)
      DIMENSION ICORE(1)
      DIMENSION NOCC(16)
      DIMENSION IADERR(200),IADFCK(200,2)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /KILLS/ KILLFLG(3)
      COMMON /POPUL/ NOCC
#include "symm2.com"
      COMMON /RPPINF/ ISTART,IORDER,IADERR,IADFCK
      COMMON /SCFPRINT/ SCFPRNT
      COMMON /GSCOMA/ GSSOPT,GSSALW,GSSALT,GSSLOK,GSSRED,GSSWRT,GSSUFR
      COMMON /GSCOMC/ ROHFMO
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
C
      DATA IONE /1/
      DATA I16 /16/
C
C  Allocate space for the AO integral sorting routine
C
C  I000 = Starting address for integral buffer (BUF)
C  I010 = Starting address for integral index buffer (IBUF)
C
      I000=1
C
      CALL DRVSRT(ICORE(I000),MAXMEM,ILNBUF,NBAS,IUHF,YESNO2)
C
C  Now for the simple stuff.  Pick up the core Hamiltonian and overlap
C  matrices.
C
      I000=1
      I010=I000+ITRILN(NIRREP+1)*IINTFP
      I020=I010+ITRILN(NIRREP+1)*IINTFP
      I030=I020+ILNBUF*IINTFP
      I040=I030+ILNBUF
      I050=I040+NBAS
      IF(I050-I000.GT.MAXMEM) THEN
        CALL NOMEM('One electron integrals','{GET1E} <-- OOFCOR <-- '//
     &             'VSCF',I050-I000,MAXMEM)
      ENDIF
C
      CALL GET1E(ICORE(I000),ICORE(I010),ICORE(I020),ICORE(I030),
     &           ICORE(I040),ITRILN(NIRREP+1),ILNBUF,REPULS,NBAS,
     &           NIRREP,NBFIRR)
C
C  Now we diagonlize the overlap matrix (stored starting at I010) and
C  obtain the S*exp(-1/2) othogonalization matrix.  This matrix is
C  returned in place of the overlap matrix.
C
      I030=I020+MXIRR2*IINTFP
      I040=I030+MXIRR2*IINTFP
      I050=I040+MXIRR2*IINTFP
      I060=I050+MXIRR2*IINTFP
      I070=I060+NBAS*NBAS*IINTFP
      I080=I070+NBAS*IINTFP
      I090=I080+NBAS
      IF(I090-I000.GT.MAXMEM) THEN
        CALL NOMEM('Form S**-1/2','{DOVRLP} <-- OOFCOR <-- VSCF',
     &             I090-I000,MAXMEM)
      ENDIF
C
      CALL DOVRLP(ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040),
     &            ICORE(I050),ICORE(I060),ICORE(I070),ICORE(I080),
     &            ITRILN(NIRREP+1),MXIRR2,NBAS)
C
C  Initialize the density matrix to zero.
C
      I030=I020+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP
      CALL ZERO(ICORE(I020),(IUHF+1)*ITRILN(NIRREP+1))

      CALL GETREC(20,'JOBARC','NAOBASFN',1,NBASX)
      CALL GETREC(20,'JOBARC','NATOMS',1,NATOMS)
      I040 = I030 + IINTFP*NBAS*NBAS
      I050 = I040 + IINTFP*MXIRR2
      I060 = I050 + IINTFP*(IUHF+1)*NBAS
      I070 = I060 + IINTFP*MXIRR2
      I080 = I070 + IINTFP*MXIRR2
      I090 = I080 + IINTFP*MXIRR2
      I100 = I090 + IINTFP*NBAS*NBAS
      I110 = I100 + IINTFP*NBAS*NBAS*(IUHF+1)
      I120 = I110 + IINTFP*NBASX*NBAS*(IUHF+1)
      I130 = I120 + IINTFP*NBASX*NBAS
      I140 = I130 + IINTFP*NBASX
      I150 = I140 + IINTFP*NBASX
      I160 = I150 + IINTFP*NBASX
      I170 = I160 + IINTFP*NBAS*(1+IUHF)
      I180 = I170 + IINTFP*NATOMS*3
      I190 = I180 + IINTFP*NATOMS*3
      I200 = I190 + IINTFP*NATOMS
      I210 = I200 + 3*NSIZE*NSIZE        * IINTFP
     &            + MAX(NBASX*10,NATOMS) * IINTFP
     &            + MAX(NBASX+MOD(NBASX,2),IINTFP*3*NATOMS)
     &            + MAX(NBASX+MOD(NBASX,2),NATOMS+MOD(NATOMS,2))
     &            + MAX(IINTFP*NBASX,NATOMS+MOD(NATOMS,2))
     &            + MAX(1080,NATOMS+MOD(NATOMS,2),IINTFP*NBASX*5)
     &            + MAX(2000,NATOMS+MOD(NATOMS,2),NBASX+MOD(NBASX,2),
     &                  NBASX*NBASX*IINTFP)
     &            + MAX(NATOMS*5+MOD(NATOMS,2),IINTFP*NBASX)
     &            + NATOMS*5+MOD(NATOMS,2)
     &            + IORDGRP              * IINTFP
      IF (I210-I000.GT.MAXMEM) THEN
         CALL NOMEM('Initial guess','{INITGES} <-- OOFCOR',
     &              I210-I000,MAXMEM)
      ELSE
         CALL INITGES(ICORE(I000),ICORE(I010),
     &                ICORE(I020),ICORE(I030),
     &                ICORE(I040),ICORE(I050),
     &                ICORE(I060),ICORE(I070),
     &                ICORE(I080),ICORE(I090),
     &                ICORE(I100),ICORE(I110),
     &                ICORE(I120),ICORE(I130),
     &                ICORE(I140),ICORE(I150),
     &                ICORE(I160),ICORE(I170),
     &                ICORE(I180),ICORE(I190),
     &                ITRILN(NIRREP+1),MXIRR2,
     &                NBAS,NBASX,NATOMS,IUHF,IOLDOCC,
     &                ICORE(I200),(MAXMEM+1-I200)/IINTFP)
      END IF
C
      NOCONV=.TRUE.
      MAXIT=IFLAGS(16)
      ITER=0
C
      if(maxit.ne.0) then
        WRITE(LUOUT,5011)
        WRITE(LUOUT,5010)
 5010   FORMAT(T3,'Iteration',9X,'Total Energy',12X,'Largest Density',
     &            ' Difference')
        WRITE(LUOUT,5011)
 5011   FORMAT(T3,68('-'))
      endif
C
C  Now we do the hard part:  construct the Fock matrix.
C
C    I000 = Core Hamiltonian
C    I010 = Orthogonalization matrix
C    I020 = Old density matrix (is updated at end of SCF cycle)
C    I030 = Fock matrix
C    I040 = Scratch space for PK lists
C
C  We also iterate over this section.  This is an SCF cycle.
C
 2000 ITER=ITER+1
      I040=I030+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP
      I050=I040+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP
      I060=I050+MXIRR2*IINTFP
      i070=i060+nbas*nbas*iintfp
      IMEM=MAXMEM-(I070-I000)
      IF(I070-I000.GT.MAXMEM) THEN
        CALL NOMEM('Fock matrix construction','{MKFOCK2} <-- OOFCOR'//
     &            ' <-- VSCF',I070-I000,MAXMEM)
      ENDIF
C
      CALL MKFOCK2(ICORE(I000),ICORE(I030),ICORE(I020),ICORE(I040),
     &             ICORE(I050),icore(i060),ICORE(I070),
     &             ITRILN(NIRREP+1),MXIRR2,IMEM,nbas,IUHF)
C
C  Now compute the total energy.
C
      I050=I040+ITRILN(NIRREP+1)*IINTFP
      I060=I050+MXIRR2*IINTFP
      I070=I060+MXIRR2*IINTFP
      IF(I070-I000.GT.MAXMEM) THEN
        CALL NOMEM('Calculate energy','{MKENER} <-- OOFCOR <-- VSCF',
     &             I070-I000,MAXMEM)
      ENDIF
C
      CALL MKENER(ICORE(I000),ICORE(I020),ICORE(I030),ICORE(I040),
     &            ICORE(I050),ICORE(I060),ITRILN(NIRREP+1),MXIRR2,
     &            REPULS,DMAX,ITER-1,IUHF,ETOT,0,TOL)
C
C  Do Davidson or static damping.
C
C  Bug Fix, Ajith Perera 08/2001. In previous incarnations the 
C  DMPSCF was called with wrong number of arguments (DEAVG 
C  variable was missing). It get initialized in vscf.F and passed 
C  in oofcor.F. 

      IF(DMPFLG)THEN
        CALL DMPSCF(ICORE(I030),ICORE(I040),ITRILN(NIRREP+1),ETOT,
     &              ITER,IUHF,DAMP,DMPFLG,DMPTOL,DEAVG)
      ENDIF
C
C  If this calculation uses DIIS, then go ahead and do it.
C
      IF(IFLAGS(10).NE.0.AND.ITER.GT.1) THEN
        NDIMA=MAX(MXIRR2,(IORDER+1)*(IORDER+2))
c       NDIMB=MAX(MXIRR2,IORDER+1)
C To screen for singularities in DIIS matrix
        NDIMB=MAX(MXIRR2,(IORDER+1)*(IORDER+1))
        I050=I040+NBAS*NBAS*IINTFP
        I060=I050+(IUHF+1)*ISQRLN(NIRREP+1)*IINTFP
        I070=I060+MXIRR2*IINTFP
        I080=I070+MXIRR2*IINTFP
        I090=I080+MXIRR2*IINTFP
        I100=I090+MAX(MXIRR2,(IORDER+1)*(IORDER+2))*IINTFP
c       I110=I100+MAX(MXIRR2,IORDER+1)*IINTFP
c To match NDIMB, above
        I110=I100+MAX(MXIRR2,(IRODER+1)*(IORDER+1))*IINTFP
        I120=I110+MXIRR2*IINTFP
        IF(I120-I000.GT.MAXMEM) THEN
          CALL NOMEM('Convergence acceleration','{DORPP} <-- OOFCOR'//
     &               ' <-- VSCF',I120-I000,MAXMEM)
        ENDIF
C
        CALL DORPP(ICORE(I000+MAXMEM),ICORE(I030),ICORE(I020),
     &             ICORE(I040),ICORE(I050),ICORE(I060),ICORE(I070),
     &             ICORE(I080),ICORE(I090),ICORE(I100),ICORE(I110),
     &             ITRILN(NIRREP+1),ISQRLN(NIRREP+1),IRPPMEM,
     &             MXIRR2,NDIMA,NDIMB,NBAS,ITER-1,IUHF)
      ENDIF
C
c     IF(ITER.GT.1.AND.IFLAGS(11).EQ.2) THEN
      IF(ROHFMO.AND.IFLAGS(11).EQ.2) THEN
C
C MAKE ROHF-FOCK MATRIX
C
       I031=I030+ITRILN(NIRREP+1)*IINTFP
       I011=I010+ITRILN(NIRREP+1)*IINTFP
       I051=I040+IINTFP*ITRILN(NIRREP+1)
       I052=I051+IINTFP*ITRILN(NIRREP+1)
       I053=I052+IINTFP*ISQRLN(NIRREP+1)
       I054=I053+IINTFP*ISQRLN(NIRREP+1)
       IEND=I054+IINTFP*ISQRLN(NIRREP+1)
       CALL MKROHFF(ICORE(I030),ICORE(I031),
     &              ICORE(I040),ICORE(I051),
     &              ICORE(I052),ICORE(I053),
     &              ICORE(I054),NBFIRR,NOCC,ITRILN(NIRREP+1),
     &              ISQRLN(NIRREP+1))
      ENDIF
C
C  Now we have fun.  Let's do the rest of the SCF iteration now.
C  Pass all this wonderful stuff into the next routine.
C
C     I000 = Core Hamiltonian
C     I010 = Orthogonalization matrix
C     I020 = Density matrix
C     I030 = Fock matrix
C     I040 = Vector for orbital energies (both alpha and beta)
C     I050 = Orbital eigenvectors
C     I060 = Scratch density matrix
C     I070 = Scratch matrix
C     I080 = Scratch matrix
C     I090 = Scratch matrix
C     I100 = Scratch matrix
C
      I050=I040+(IUHF+1)*NBAS*IINTFP
      I060=I050+(IUHF+1)*ISQRLN(NIRREP+1)*IINTFP
      I070=I060+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP
      I080=I070+MAX(2*NBAS,NBAS*NBAS,MXIRR2)*IINTFP
      I090=I080+MAX(2*NBAS,NBAS*NBAS,MXIRR2)*IINTFP
      I100=I090+MAX(2*NBAS,NBAS*NBAS,MXIRR2)*IINTFP
      I110=I100+MAX(2*NBAS,MXIRR2)*IINTFP
C
      IF(GSSLOK)THEN
       I120 = I110 + NBAS * NBAS * IINTFP
       I130 = I120 + NBAS * NBAS * IINTFP
      ELSE
       I110 = I090
       I120 = I100
       I130 = I100 + MAX(2*NBAS,MXIRR2)*IINTFP
      ENDIF
C
      IF(I130-I000.GT.MAXMEM) THEN
        CALL NOMEM('SCF iteration','{SCFIT} <-- OOFCOR <-- VSCF',
     &             I130-I000,MAXMEM)
      ENDIF
C
cjp really do not perform any single iteration when requested
      if(iflags(16).ne.0)then
cjp
      CALL SCFIT(ICORE(I000),ICORE(I010),ICORE(I020),ICORE(I030),
     &           ICORE(I040),ICORE(I050),ICORE(I060),ICORE(I070),
     &           ICORE(I080),ICORE(I090),ICORE(I100),ICORE(I110),
     &           ICORE(I120),REPULS,
     &           ITRILN(NIRREP+1),MXIRR2,ISQRLN(NIRREP+1),NOCONV,
     &           TOL,ITER,DMAX,IUHF,NBAS)
cjp
      else
      call reload(icore(i050),icore(i040),icore(i070),icore(i080),
     &            icore(i090),isqrln(nirrep+1),nbas,iuhf)
       noconv=.false.
       endif
cjp

C
      IF(NOCONV.AND.ITER.LT.MAXIT) GOTO 2000
C
cjp
cjp modified diagnostics output
cjp
      IF (NOCONV.or.iflags(16).eq.0) THEN
cjp scf process skipped when maxit==0
         if (iflags(16).eq.0) then
            write(*,*) '\n@VSCF: SCF has been omitted.\n'
         else
            if (iflags(16).eq.1) then
               write(*,*) '\n@VSCF: SCF 1 iteration performed.\n'
            else
               write(*,*) '\n@VSCF: SCF has not converged!!!\n'
               KILLFLG(1)=1
            end if
         end if
cjp
      ELSE
         WRITE(*,*) '\n@OOFCOR: SCF has converged.\n'
      END IF
C
C  Construct the final fock matrix and the final energy.
C
      I070=I060+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP
      I080=I070+MXIRR2*IINTFP
      i090=i080+nbas*nbas*iintfp
      IMEM=MAXMEM-(I090-I000)
      IF(I090-I000.GT.MAXMEM) THEN
        CALL NOMEM('Final Fock matrix','{MKFOCK2} <-- OOFCOR '//
     &             '<-- VSCF',I100-I000,MAXMEM)
      ENDIF
C
      CALL MKFOCK2(ICORE(I000),ICORE(I030),ICORE(I020),ICORE(I060),
     &             ICORE(I070),icore(i080),ICORE(I090),
     &             ITRILN(NIRREP+1),MXIRR2,IMEM,nbas,IUHF)
C
      I070=I060+ITRILN(NIRREP+1)*IINTFP
      I080=I070+MXIRR2*IINTFP
      I090=I080+MXIRR2*IINTFP
      IF(I090-I000.GT.MAXMEM) THEN
        CALL NOMEM('Final energy','{MKENER} <-- OOFCOR <-- VSCF',
     &             I090-I000,MAXMEM)
      ENDIF
C
C GET ORBITAL ENERGIES FOR ROHF
C
      IF(IFLAGS(11).EQ.2) THEN
       I031=I030+ITRILN(NIRREP+1)*IINTFP
       I051=I070
       I052=I051+ISQRLN(NIRREP+1)*IINTFP
       I053=I052+ISQRLN(NIRREP+1)*IINTFP
       I054=I053+ISQRLN(NIRREP+1)*IINTFP
       CALL GETEVAL(ICORE(I030),ICORE(I031),ICORE(I040),
     &              ICORE(I051),ICORE(I052),ICORE(I053),
     &              ITRILN(NIRREP+1),ISQRLN(NIRREP+1),NBFIRR)
      ENDIF

C
      CALL MKENER(ICORE(I000),ICORE(I020),ICORE(I030),ICORE(I060),
     &            ICORE(I070),ICORE(I080),ITRILN(NIRREP+1),MXIRR2,
     &            REPULS,DMAX,ITER,IUHF,ETOT,0,TOL)
C
C  If the SCF occupancies have changed, print 'em out.
C
      CHGOCC=.FALSE.
      DO 600 I=1,16
        IF(NOCC(I).NE.IOLDOCC(I)) CHGOCC=.TRUE.
  600 CONTINUE
      IF(CHGOCC) THEN
        WRITE(LUOUT,5100)
        WRITE(LUOUT,5101)(NOCC(I),I=1,NIRREP)
        WRITE(LUOUT,5102)(NOCC(8+I),I=1,NIRREP)
        WRITE(LUOUT,5103)
 5100   FORMAT(/,T3,'@OOFCOR-I, Final occupancies:',/)
 5101   FORMAT(T8,'   Alpha population by irrep: ',8(I3,2X))
 5102   FORMAT(T8,'    Beta population by irrep: ',8(I3,2X))
 5103   FORMAT(/)
      ENDIF
C
C  For open-shell singlet calculations checks phase of orbitals
C
CJDW 6/24/96. Modification : redefinition of CHKPHASE. See VSCF for
C                            explanation.
C
c& new lines
      IDELOCC = 0
      DO I=1, NIRREP
        IDELOCC = IDELOCC + IABS(NOCC(I) - NOCC(8+I))
      ENDDO
      CALL GETREC(-1,'JOBARC','OSCALC  ',1,IOS)
cjdw  CHKPHASE = (IOS .NE. 0) .OR. ((IUHF .GE. 1) .AND.
cjdw &   (IDELOCC .EQ. 0))
      CHKPHASE =  IOS .NE. 0
      IF(CHKPHASE) THEN
         I070=I060+NBAS*NBAS*IINTFP
         I080=I070+MXIRR2*IINTFP
         I090=I080+MXIRR2*IINTFP
         MXMEM=MAXMEM-(I090-I000)
         IF(I090-I000.GT.MAXMEM) THEN
          CALL NOMEM('Rephase orbitals','{OSPHASE} <-- OOFCOR'//
     &         ' <-- VSCF',I090-I000,MAXMEM)
          ENDIF
C
          CALL OSPHASE(ICORE(I050),ICORE(I060),ICORE(I070),
     &       ICORE(I080),ISQRLN(NIRREP+1),MXIRR2,ITRILN(NIRREP+1),NBAS,
     &       ICORE(I090),MXMEM,IUHF,ICORE(I020),DUMMY,ICORE(I000),
     &       ICORE(I030),IPKSIZ,.FALSE.)
       ENDIF
c& end
C
C  Now print out the eigenvectors, if the user wants them.
C
      IF (SCFPRNT.NE.0)
     &   CALL OUTEVC(ICORE(I040),ICORE(I050),
     &               NBAS,ISQRLN(NIRREP+1),IUHF,.FALSE.)
C
C  Now calculate the value of S**2, if UHF.
C
      IF(IUHF.EQ.1) THEN
        I070=I060+NBAS*NBAS*IINTFP
        I080=I070+MXIRR2*IINTFP
        I090=I080+MXIRR2*IINTFP
        IF(I090-I000.GT.MAXMEM) THEN
          CALL NOMEM('Calculate S**2 and multiplicity','{S2CALC} <--'//
     &               ' OOFCOR <-- VSCF',I090-I000,MAXMEM)
        ENDIF
C
        CALL S2CALC(ICORE(I050),ICORE(I060),ICORE(I070),
     &              ICORE(I080),ISQRLN(NIRREP+1),MXIRR2,NBAS)
      ENDIF
C
C  Dump out the SCF occupations for use in subsequent steps.
C
      CALL PUTREC(20,'JOBARC','OCCUPYA0',NIRREP,NOCC(1))
      IF(IUHF.NE.0) THEN
        CALL PUTREC(20,'JOBARC','OCCUPYB0',NIRREP,NOCC(9))
      ENDIF
CMN
CMN  transform occ and virt orbitals to v(N-1) potential if
CMN  requested
CMN
      IF (IFLAGS2(121) .GE. 1) THEN
        IO70 = I060
      I080=I070+MAX(2*NIRREP,NBAS+MOD(NBAS,2))
      I090=I080+5*ISQRLN(NIRREP+1)*IINTFP
      I100=I090+ITRILN(NIRREP+1)*IINTFP
      IF(i100-I000.GT.MAXMEM) THEN
        CALL NOMEM('NEW VIRTUALS','{DRNEWVRT} <-- OOFCOR',i100-I000,
     &             MAXMEM)
      ENDIF
C
      CALL DRNEWVRT(icore(i000),icore(i010),ICORE(I020),ICORE(I030),
     &            ICORE(I040),
     &            ICORE(I050),ICORE(I070),ICORE(I080),ICORE(I090),
     &            ITRILN(NIRREP+1),
     &            ISQRLN(NIRREP+1),MXIRR2,NBAS,IUHF)
C
C
C  For open-shell singlet calculations checks phase of orbitals
C
c& new lines
      CALL GETREC(-1,'JOBARC','OSCALC  ',1,IOS)
      IF(IOS.NE.0) THEN
         I070=I060+NBAS*NBAS*IINTFP
         I080=I070+MXIRR2*IINTFP
         I090=I080+MXIRR2*IINTFP
         MXMEM=MAXMEM-(I090-I000)
         IF(I090-I000.GT.MAXMEM) THEN
          CALL NOMEM('Rephase orbitals','{OSPHASE} <-- OOFCOR'//
     &         ' <-- VSCF',I090-I000,MAXMEM)
          ENDIF
C
          CALL OSPHASE(ICORE(I050),ICORE(I060),ICORE(I070),
     &       ICORE(I080),ISQRLN(NIRREP+1),MXIRR2,ITRILN(NIRREP+1),NBAS,
     &       ICORE(I090),MXMEM,IUHF,ICORE(I020),DUMMY,ICORE(I000),
     &       ICORE(I030),IPKSIZ,.FALSE.)
       ENDIF
c& end
C
C  Now print out the eigenvectors, if the user wants them.
C
      IF (SCFPRNT.NE.0)
     &   CALL OUTEVC(ICORE(I040),ICORE(I050),
     &               NBAS,ISQRLN(NIRREP+1),IUHF,.FALSE.)

      ENDIF
C
CMN END
C
C  If this is a QRHF calculation, set it up.
C
      IF(IFLAGS(77).NE.0) THEN
        CALL QRHFIT2(ICORE(I000),ICORE(I020),ICORE(I030),ICORE(I040),
     &               ICORE(I050),
     &               ICORE(I060),MAXMEM-(I060-I000),ITRILN(NIRREP+1),
     &               ISQRLN(NIRREP+1),NBAS,MXLSSZ1,REPULS,ETOT,IUHF,IOS)
      ENDIF
C
CMN Add UNO-REF possibility
C
C     I000 = Core Hamiltonian
C     I010 = Orthogonalization matrix
C     I020 = Density matrix
C     I030 = Fock matrix
C     I040 = Vector for orbital energies (both alpha and beta)
C     I050 = Orbital eigenvectors
C     I060 = Scratch density matrix
c
      LDIM2 = ISQRLN(NIRREP+1)
      IDENSNEW=I060
      IDENSSCR=IDENSNEW+2*LDIM2*IINTFP
      ISCRA=IDENSSCR+2*LDIM2*IINTFP
      ISCRB=ISCRA+2*LDIM2*IINTFP
      IEND=ISCRB+2*LDIM2*IINTFP

      IF(IFLAGS2(148).NE.0) THEN
        CALL UNO_REF2(ICORE(I000),ICORE(I020),ICORE(I030),
     &              ICORE(I040),ICORE(I050),ICORE(IEND),
     &              MAXMEM-(IEND-I000),ITRILN(NIRREP+1),
     &              ISQRLN(NIRREP+1),NBAS,MXLSSZ1,REPULS,ETOT,IUHF,
     $        ICORE(I010), ICORE(IDENSNEW), ICORE(IDENSSCR),
     $        ICORE(ISCRA), ICORE(ISCRB))
      ENDIF
C
C  Now dump the stuff to JOBARC.
C
      I070=I060+NBAS*NBAS*IINTFP
      I080=I070+MXIRR2*IINTFP
      I090=I080+NBAS+MOD(NBAS,2)
      I100=I090+4*ISQRLN(NIRREP+1)*IINTFP
      I110=I100+NBAS*NBAS*IINTFP
      IF(I110-I000.GT.MAXMEM) THEN
        CALL NOMEM('Dump to JOBARC','{DMPJOB} <-- OOFCOR <-- VSCF',
     &             I110-I000,MAXMEM)
      ENDIF
C
      CALL DMPJOB(icore(i000),icore(i010),ICORE(I020),ICORE(I030),
     &            ICORE(I040),ICORE(I050),ICORE(I060),ICORE(I070),
     &            ICORE(I080),ICORE(I090),icore(i100),ITRILN(NIRREP+1),
     &            ISQRLN(NIRREP+1),MXIRR2,NBAS,ETOT,
     &            IUHF,SEWARD)
C
C  Now form the eigenvectors in the original AO basis (not the SO basis).
C
       CALL GETREC(20,'JOBARC','NAOBASFN',1,NBASX)
       CALL GETREC(20,'JOBARC','NATOMS',1,NATOMS)
       I010 = I000 + NBAS *NBAS*IINTFP
       I020 = I010 + NBASX*NBAS*IINTFP
       I030 = I020 + NBASX*NBAS*IINTFP
       I040 = I030 + NBAS *  2 *IINTFP
       I050 = I040 + NATOMS    *IINTFP
       I060 = I050 + NBASX     *IINTFP
       I070 = I060 + NBASX     *IINTFP
       I080 = I070 + NATOMS    *IINTFP
       I090 = I080 + NATOMS    *IINTFP
       I100 = I090 + NATOMS    *IINTFP
       I110 = I100 + NATOMS    *IINTFP
       I140 = I110 + NATOMS    *IINTFP
       I150 = I140 + NBASX     *IINTFP
       I160 = I150 + NBAS * 2  *IINTFP
       I170 = I160 + 3*NATOMS  *IINTFP
       IF(I170-I000.GT.MAXMEM) THEN
         CALL NOMEM('Eigenvectors (SO --> AO)','{EVCAO2} <-- OOFCOR',
     &              I140-I000,MAXMEM)
       ENDIF
       CALL EVCAO2(ICORE(I000),ICORE(I010),ICORE(I020),ICORE(I030),
     &             ICORE(I040),ICORE(I050),ICORE(I060),ICORE(I070),
     &             ICORE(I080),ICORE(I090),ICORE(I100),ICORE(I110),
     &                                     ICORE(I140),NOCC,
     &             ICORE(I150),ICORE(I160),NBASX,      NBAS,
     &             NATOMS,     IUHF)

C  Now we do some fun stuff.  We determine the irreps of each of
C  the eigenvectors.
C
C
CJDW 3/28/97. Check whether this is the first entry to xvscf.
C             FIRSTSCF is used to control whether we print certain
C             things or not.
C
      CALL GETREC(1,'JOBARC','FIRSTRUN',1,I)
      FIRSTSCF = (i.ne.0)

      CALL GETREC(20,'JOBARC','JODAFLAG',I16,ICORE(I000))
      IF ((ICORE(I000+8).NE.2) .AND. FIRSTSCF) THEN
         CALL IRRORB_SHELL(ICORE(I000),MAXMEM,NBAS,IUHF)
         CALL PRJDENS(ICORE(I000),MAXMEM,NBAS,IUHF)
c dump the reordered eigenvalues picking them up from JOBARC
         I010 = I000 + (IUHF+1)*NBAS*IINTFP
         I020 = I010 + NBAS+MOD(NBAS,2)
         I030 = I020 + NBAS+MOD(NBAS,2)
         IF (I030-I000.GT.MAXMEM) THEN
            CALL NOMEM('Dump eigenvalues',
     &                 '{EVLOUT} <-- OOFCOR <-- VSCF',I030-I000,MAXMEM)
         ELSE
            CALL EVLOUT(ICORE(I000),ICORE(I010),ICORE(I020),
     &                  IREPS,NBAS,NIRREP,IUHF)
         END IF
      END IF

C The SA-AO occupied-occupied block of the Fock matrix is need for the 
C gradients. Built it here. This is need for SCF/Correlated gradients
C with Alaska code. 07/2000, Ajith Perera
C
      if (seward) then
         postscf = (iflags(2).ne.0)
         I010 = I000 + NBAS*NBAS*IINTFP
         I020 = I010 + NBAS*NBAS*IINTFP
         I030 = I020 + NBAS*NBAS*IINTFP
         I040 = I030 + 2*NBAS*NBAS*IINTFP
         I050 = I040 + NBAS*NBAS*IINTFP
         IF ((I050-I000).GE.MAXMEM) THEN
            CALL INSMEM("vscf", (I050-I000), MAXMEM)
         ELSE
            CALL BLT_EFCTVFOCK(ICORE(I000), ICORE(I010), ICORE(I020),
     &                         ICORE(I030), ICORE(I040), NBAS, NOCC,
     &                         NIRREP, NBFIRR, IUHF, POSTSCF)
         END IF
      end if

C
      RETURN
      END
