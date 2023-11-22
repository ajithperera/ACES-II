      SUBROUTINE EHTGSS(DENS,FOCK,SOVRLP,SCR1,SCR2,SCR3,EIGVAL,
     &                  NBAS,NBASX,LDIM1,LDIM2,IUHF)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Arguments.
C-----------------------------------------------------------------------
      DOUBLE PRECISION DENS,FOCK,SOVRLP,SCR1,SCR2,SCR3,EIGVAL
      INTEGER NBAS,NBASX,LDIM1,LDIM2,IUHF
C-----------------------------------------------------------------------
C     Parameter.
C-----------------------------------------------------------------------
      INTEGER MXNATOM
C-----------------------------------------------------------------------
C     Local variables.
C-----------------------------------------------------------------------
      CHARACTER*5 ITYPE
      DOUBLE PRECISION HALF
      DOUBLE PRECISION ECOR,EVAL,BETA0, DOCC
      INTEGER NATOMS,IZ,NCOR,NVAL,IATOM,JATOM,IOFF,IVAL,JVAL,IBAS,JBAS
      INTEGER ICNT,I,J,N,IROFF
      integer ehtnbfirr(8),ehtireps(9)
      INTEGER NBASXTRUE
      LOGICAL ISTHER,MBS
C-----------------------------------------------------------------------
C     Common block variables.
C-----------------------------------------------------------------------
      INTEGER IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      INTEGER NOCC, NSUM
C-----------------------------------------------------------------------
C     Functions.
C-----------------------------------------------------------------------
      INTEGER INDX2
C-----------------------------------------------------------------------
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
#include "mxatms.par"
#include "symm2.com"
#include "maxbasfn.par"
      COMMON /POPUL/  NOCC(16)
C-----------------------------------------------------------------------
      DIMENSION DENS((IUHF+1)*LDIM1),FOCK(NBASX*NBASX),
     &          SOVRLP(NBASX*NBASX), DOCC(MAXBASFN*2),
     &          SCR1(NBASX*NBASX),SCR2(NBASX*NBASX),SCR3(LDIM2),
     &          EIGVAL(NBAS*(IUHF+1))
      DIMENSION IZ(MXATMS),NCOR(MXATMS),NVAL(MXATMS),IOFF(MXATMS)
      DIMENSION ECOR(5,18),EVAL(4,18),BETA0(18),NSUM(16)
C-----------------------------------------------------------------------
      DATA HALF / 0.5D+00 /
C-----------------------------------------------------------------------
      INDX2(I,J,N)=I+(J-1)*N
C-----------------------------------------------------------------------
C
C     Important:
C
C     Despite the notation, this routine expects NBAS and NBASX to refer
C     to the EHT basis set, i.e. a minimal basis set. The true NBAS and
C     NBASX (as well as related quantities such as NBFIRR and various
C     offsets) refer to the ab initio basis set being used.
C
C     This routine has various limitations when the ab initio basis set
C     is not STO-3G, particularly CMP2ZMAT is not available, so only C1
C     calculations can be run.
C
C     If the ab initio basis set is STO-3G, this routine can run with
C     symmetry and all sorts of benefits arise. One strategy is the
C     two-step procedure: run STO-3G, using EHT as the initial guess.
C     The MOs will be written to either MOFILE or EHTMOFILE, and they
C     can be read in a subsequent larger basis set calculation.
C
C-----------------------------------------------------------------------
C     Determine if this is a minimal basis set run. [This is not a
C     rock-solid test, but it should get us by for now.]
C-----------------------------------------------------------------------
     
      CALL GETREC(20,'JOBARC','NAOBASFN',1,NBASXTRUE)
      IF(NBASXTRUE.EQ.NBASX)THEN
       write(6,*) ' @EHTGSS-I, Apparently this is an MBS calculation '
       MBS = .TRUE.
      ELSE
       write(6,*) ' @EHTGSS-I, This is not an MBS calculation '
       MBS = .FALSE.
      ENDIF
C
      IF(.NOT.MBS)THEN
       IF(NIRREP.NE.1)THEN
        write(6,*) ' @EHTGSS-F, Program can only handle C1 currently '
        call errex
       ENDIF
      ENDIF
C     
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C     Set atomic parameters.
C-----------------------------------------------------------------------
      CALL EHTPARAM(ECOR,EVAL,BETA0)
C
C-----------------------------------------------------------------------
C     Get the overlap matrix and put it in the ZMAT-ordered AO basis set.
C
C     We expect the caller to supply the overlap matrix in SOVRLP (in
C     the ZMAT-ordered AO basis set).
C
C     The following block reads the overlap matrix from JOBARC. This is
C     only meaningful in an STO-3G calculation. The idea is that the
C     caller-supplied overlap matrix is independent of basis set: calc_s
C     calculates the overlap matrix for EHT STOs, and this is what
C     we normally expect to use. However, the following block may be
C     useful for debug/checking purposes.
C-----------------------------------------------------------------------
c      CALL GETREC(20,'JOBARC','AOOVRLAP',NBAS*NBAS*IINTFP,SOVRLP)
c      write(6,*) ' @EHTGSS-I, Overlap matrix in SO basis '
c      call output(sovrlp,1,nbas,1,nbas,nbas,nbas,1)
c      CALL GETREC(20,'JOBARC','ZMAT2CMP',NBAS*NBASX*IINTFP,SCR1)
c      write(6,*) ' @EHTGSS-I, Transformation matrix '
c      call output(scr1,1,nbas,1,nbasx,nbas,nbasx,1)
c      CALL XGEMM('N','N',NBAS,NBASX,NBAS,1.0D+00,SOVRLP,NBAS,
c     &           SCR1,NBAS,0.0D+00,SCR2,NBAS)
c      write(6,*) ' @EHTGSS-I, SCR2 after first XGEMM '
c      call output(scr2,1,nbas,1,nbas,nbas,nbas,1)
c      CALL XGEMM('T','N',NBASX,NBASX,NBAS,1.0D+00,SCR1,NBAS,
c     &           SCR2,NBAS,0.0D+00,SOVRLP,NBASX)
c      write(6,*) ' @EHTGSS-I, Overlap matrix in AO basis '
c      call output(sovrlp,1,nbasx,1,nbasx,nbasx,nbasx,1)
C
C-----------------------------------------------------------------------
C     Set number of core and valence orbitals on each atom.
C-----------------------------------------------------------------------
C
      CALL GETREC(20,'JOBARC','NATOMS',1,NATOMS)
      CALL GETREC(20,'JOBARC','ATOMCHRG',NATOMS,IZ)
      DO  50 IATOM=1,NATOMS
C
      IF(IZ(IATOM).GE. 1 .AND. IZ(IATOM).LE. 2) ITYPE = 'FIRST'
      IF(IZ(IATOM).GE. 3 .AND. IZ(IATOM).LE.10) ITYPE = 'SECND'
      IF(IZ(IATOM).GE.11 .AND. IZ(IATOM).LE.18) ITYPE = 'THIRD'
C
      IF(ITYPE.EQ.'FIRST')THEN
       NCOR(IATOM) = 0
       NVAL(IATOM) = 1
      ELSEIF(ITYPE.EQ.'SECND')THEN
       NCOR(IATOM) = 1
       NVAL(IATOM) = 4
      ELSEIF(ITYPE.EQ.'THIRD')THEN
       NCOR(IATOM) = 5
       NVAL(IATOM) = 4
      ENDIF
   50 CONTINUE
C
C-----------------------------------------------------------------------
C     Atomic offsets in basis set list.
C-----------------------------------------------------------------------
      IOFF(1) = 0
      IF(NATOMS.GT.1)THEN
       DO 60 IATOM=2,NATOMS
       IOFF(IATOM) = IOFF(IATOM-1) + NCOR(IATOM-1) + NVAL(IATOM-1)
   60  CONTINUE
      ENDIF
C
C-----------------------------------------------------------------------
C     Build the EHT Fock matrix in the AO basis.
C-----------------------------------------------------------------------
      CALL ZERO(FOCK,NBASX*NBASX)
      DO 100 JATOM=1,NATOMS
      DO  90 IATOM=1,NATOMS
C
      DO  80 JVAL =1,NCOR(JATOM)+NVAL(JATOM)
      DO  70 IVAL =1,NCOR(IATOM)+NVAL(IATOM)
C
      IBAS = IOFF(IATOM) + IVAL
      JBAS = IOFF(JATOM) + JVAL
C
C     Core-core :
C
      IF(IVAL.LE.NCOR(IATOM) .AND. JVAL.LE.NCOR(JATOM))THEN
C
       IF(IBAS .EQ. JBAS)THEN
        FOCK( (JBAS-1)*NBASX + IBAS ) = ECOR(IVAL,IZ(IATOM))
       ELSE
        FOCK( (JBAS-1)*NBASX + IBAS ) = SOVRLP( (JBAS-1)*NBASX + IBAS )
     &       * ( ECOR(JVAL,IZ(JATOM)) + ECOR(IVAL,IZ(IATOM)) ) * HALF
       ENDIF
C
      ENDIF
C
C     Valence-valence :
C
      IF(IVAL.GT.NCOR(IATOM) .AND. JVAL.GT.NCOR(JATOM))THEN
C
       IF(IBAS .EQ. JBAS)THEN
        FOCK( (JBAS-1)*NBASX + IBAS ) = EVAL(IVAL-NCOR(IATOM),IZ(IATOM))
       ELSE
        FOCK( (JBAS-1)*NBASX + IBAS ) = SOVRLP( (JBAS-1)*NBASX + IBAS )
     &       * ( BETA0(IZ(IATOM)) + BETA0(IZ(JATOM)) ) * HALF
       ENDIF
C
      ENDIF
C
   70 CONTINUE
   80 CONTINUE
C
   90 CONTINUE
  100 CONTINUE
C
      write(6,*) ' @EHTGSS-I, Extended Huckel Fock matrix. AO basis '
      call output(FOCK,1,NBASX,1,NBASX,NBASX,NBASX,1)

C-----------------------------------------------------------------------
C     Transform the matrix to the SO basis set.
C
C     Serious symmetry restrictions here!!!!!!!!!!!
C-----------------------------------------------------------------------
C
      EHTNBFIRR(1) = NBAS
      EHTIREPS(1)  = 1
      EHTIREPS(2)  = 1 + NBAS




      IF(MBS)THEN
       CALL GETREC(20,'JOBARC','CMP2ZMAT',NBASX*NBAS*IINTFP,SCR1)
      ELSE
       CALL EHTTRAN(SCR1,SCR2,IZ,NBASX,NATOMS)
      ENDIF
      write(6,*) ' @EHTGSS-I, CMP2ZMAT Transformation matrix '
      call output(scr1,1,nbas,1,nbasx,nbas,nbasx,1)
      CALL XGEMM('N','N',NBASX,NBAS,NBASX,1.0D+00,FOCK,NBASX,
     &           SCR1,NBASX,0.0D+00,SCR2,NBASX)
c     write(6,*) ' @EHTGSS-I, SCR2 after first XGEMM '
c     call output(scr2,1,nbas,1,nbas,nbas,nbas,1)
      CALL XGEMM('T','N',NBAS,NBAS,NBASX,1.0D+00,SCR1,NBASX,
     &           SCR2,NBASX,0.0D+00,FOCK,NBAS)
C
      write(6,*) ' @EHTGSS-I, Extended Huckel Fock matrix. SO basis '
      call output(FOCK,1,NBAS,1,NBAS,NBAS,NBAS,1)
C
C-----------------------------------------------------------------------
C     Note that we are throwing away the EHT overlap matrix here. This
C     might be worth changing later (or the matrix could be saved).
C-----------------------------------------------------------------------
C
      ICNT=0
      CALL ZERO(SOVRLP,NBAS*NBAS)
      IF(MBS)THEN
       DO 120 I=1,NIRREP
        IF(NBFIRR(I).EQ.0)GOTO 120
        CALL GETBLK(FOCK,SCR2,NBFIRR(I),NBAS,IREPS(I))
        CALL ZERO(SCR1,NBAS*NBAS)
        CALL EIG(SCR2,SCR1,NBFIRR(I),NBFIRR(I),0)
C
C     Save MOs in SOVRLP.
C
        CALL PUTBLK(SOVRLP,SCR1,NBFIRR(I),NBAS,IREPS(I))
C
        DO 110 J=1,NBFIRR(I)
         ICNT=ICNT+1
         EIGVAL(ICNT)=SCR2(INDX2(J,J,NBFIRR(I)))
  110   CONTINUE
  120  CONTINUE
      ELSE
       DO 200 I=1,NIRREP
        IF(EHTNBFIRR(I).EQ.0)GOTO 200
        CALL GETBLK(FOCK,SCR2,EHTNBFIRR(I),NBAS,EHTIREPS(I))
        CALL ZERO(SCR1,NBAS*NBAS)
        CALL EIG(SCR2,SCR1,EHTNBFIRR(I),EHTNBFIRR(I),0)
C
C     Save MOs in SOVRLP.
C
        CALL PUTBLK(SOVRLP,SCR1,EHTNBFIRR(I),NBAS,EHTIREPS(I))
C
        DO 190 J=1,EHTNBFIRR(I)
         ICNT=ICNT+1
         EIGVAL(ICNT)=SCR2(INDX2(J,J,EHTNBFIRR(I)))
  190   CONTINUE
  200  CONTINUE
      ENDIF
      IF(IUHF.EQ.1) THEN
       CALL ICOPY(NBAS*IINTFP,EIGVAL(1),1,EIGVAL(NBAS+1),1)
      ENDIF
C
C-----------------------------------------------------------------------
C     Determine the EHT occupation if none has been specified.
C     A little messy right now. Need to clean up here and in caller.
C-----------------------------------------------------------------------
C
      ISTHER = .FALSE.
      DO 210 I=1,(IUHF+1)*8
      IF(NOCC(I).NE.0) ISTHER = .TRUE.
  210 CONTINUE
C
      IF(.NOT.ISTHER)THEN
       IF(MBS)THEN
        CALL OCCUPY(NIRREP,NBFIRR,NBAS,EIGVAL,SCR1,IUHF)
       ELSE
        CALL OCCUPY(NIRREP,EHTNBFIRR,NBAS,EIGVAL,SCR1,IUHF)
       ENDIF
       WRITE(6,5000)
      ENDIF
C
      WRITE(6,5001)(NOCC(I),I=1,NIRREP)
      WRITE(6,5002)(NOCC(8+I),I=1,NIRREP)
      WRITE(6,5010)
 5000 FORMAT(T3,'@EHTGSS-I, Occupancies from EHT Hamiltonian:',/)
 5001 FORMAT(T8,'   Alpha population by irrep: ',8(I3,2X))
 5002 FORMAT(T8,'    Beta population by irrep: ',8(I3,2X))
 5010 FORMAT(/)
C
C-----------------------------------------------------------------------
C     Orthonormalize the MOs.
C-----------------------------------------------------------------------

      CALL ICOPY(NBAS*NBAS*IINTFP,SOVRLP,1,SCR1,1)
C
      write(6,*) ' @EHTGSS-I, MO coefficients from EHT '
      call output(SCR1,1,NBAS,1,NBAS,NBAS,NBAS,1)
C
C     Dangerous right now --- SORTHO presumably expects different dim.
C     Even for STO-3G, how meaningful is this? In favor: we are taking
C     EHT MO coefficients and creating a set that is orthonormal with
C     respect to the STO-3G overlap matrix. Then we use this set to
C     create an STO-3G initial density matrix.
C      
      IF(MBS)THEN
       CALL SORTHO(SCR1,SOVRLP,SCR2,SCR3,FOCK,LDIM2,NBAS,.TRUE.,NOCC)
       write(6,*) ' @EHTGSS-I, MOs from EHT after sortho. SO basis '
       call output(SCR1,1,NBAS,1,NBAS,NBAS,NBAS,1)
      ENDIF

C
      IF(MBS)THEN
       CALL GETREC(20,'JOBARC','CMP2ZMAT',NBASX*NBAS*IINTFP,SOVRLP)
      ELSE
       CALL EHTTRAN(SOVRLP,SCR2,IZ,NBASX,NATOMS)
      ENDIF
      CALL XGEMM('N','N',NBASX,NBAS,NBAS,1.0D+00,SOVRLP,NBASX,
     &           SCR1,NBAS,0.0D+00,SCR2,NBASX)
      write(6,*) ' @EHTGSS-I, MOs from EHT (sortho). AO basis '
      call output(SCR2,1,NBASX,1,NBASX,NBASX,NBASX,1)
C
C-----------------------------------------------------------------------
C     Write the EHT MOs and other information to EHTMOFILE. The format
C     is exactly the same as for AOBASMOS (written by evcao2), but the
C     dimensions are for a minimal basis set.
C-----------------------------------------------------------------------
      write(6,*) ' @EHTGSS-I, Calling EHTDMPMOS '
      CALL EHTDMPMOS(SCR2,EIGVAL,NOCC,EHTNBFIRR,NBASX,NBAS,NIRREP,
     &               NATOMS,IUHF)
      write(6,*) ' @EHTGSS-I, Back from EHTDMPMOS '
C
C-----------------------------------------------------------------------
C     Make the density matrix. [Note: this is the density matrix for
C     the EHT Huckel orbitals. Without some projection scheme it is
C     NOT an initial guess for anything but a minimal basis set]
C-----------------------------------------------------------------------
      IF(MBS)THEN 
#ifdef _DPOCC
        CALL MKDP_OCCNOS(DOCC, NOCC, NSUM, NBAS, IUHF, .TRUE.)
#endif
        IROFF = 0
        DO 810 I=1,NIRREP
        CALL GETBLK(SCR1,SCR2,NBFIRR(I),NBAS,IREPS(I))
C
#ifdef _DPOCC
        IOFF = IROFF + 1
        CALL MKDDEN_4IRREP(SCR2,SCR3,SCR1,NSUM(I),DOCC(IOFF),
     &                     NBFIRR(I),1,IUHF)
        IROFF = IROFF + NBFIRR(I)
#else
        CALL MKDEN(SCR2,SCR3,NOCC(I),NBFIRR(I),IUHF)
#endif
        CALL SQUEZ2(SCR3,DENS(ITRIOF(I)),NBFIRR(I))
  810  CONTINUE
C
       IF(IUHF.EQ.1)THEN
        IROFF = 0
        DO 820 I=1,NIRREP
         IOFF = IROFF  + NBAS + 1
         CALL GETBLK(SCR1,SCR2,NBFIRR(I),NBAS,IREPS(I))
#ifdef _DPOCC
         CALL MKDDEN_4IRREP(SCR2,SCR3,SCR1, NSUM(8+I), 
     &                      DOCC(IOFF),NBFIRR(I),2,IUHF)
         IROFF = IROFF + NBFIRR(I)
#else
         CALL MKDEN(SCR2,SCR3,NOCC(8+I),NBFIRR(I),IUHF)
#endif
         CALL SQUEZ2(SCR3,DENS(LDIM1+ITRIOF(I)),NBFIRR(I))
  820   CONTINUE
       ENDIF
      ENDIF
C
      RETURN
      END
