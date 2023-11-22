      SUBROUTINE RELOAD (EVEC,EVAL,SCR1,SCRN,SCR2,LDIM3,NBAS,IUHF)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION SCRN(NBAS*NBAS), SCR2(NBAS*NBAS)
      DIMENSION SCR1(NBAS*NBAS),EVEC(LDIM3*(1+IUHF)),EVAL(NBAS*(1+IUHF))
      DIMENSION JFLAG(16)
      LOGICAL REORDER
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
#include "symm2.com"
      COMMON /POPUL/ NOCC(16)
C
      one   = 1.0d0 
      zilch = 0.0d0
c
      CALL GETREC(-1,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCR1)
      CALL GETREC(-1, 'JOBARC', 'FOCKA   ', NBAS*NBAS*IINTFP,SCRN)
      call xgemm('N', 'N', nbas, nbas, nbas, one, scrn, nbas,
     $     scr1, nbas, zilch, scr2, nbas)
      call xgemm('T', 'N', nbas, nbas, nbas, one, scr1, nbas,
     $     scr2, nbas, zilch, scrn, nbas)
c
c pick up diagonal elements
c
      call SCOPY(nbas, scrn, nbas+1, scr2, 1)
      call putrec(20, 'JOBARC', 'SCFEVLA0', NBAS*IINTFP, scr2)
#ifdef _DEBUG_LVLM
      write(*,*) "The eigenvalues"
      call checksum("MO-FOCK", scrn, Nbas*Nbas)
      write(6,"(6(1x,F12.6))") (scr2(i),i=1,nbas)
#endif 
      if (iuhf .eq. 0) then
      call putrec(20, 'JOBARC', 'SCFEVLB0', NBAS*IINTFP, scr2)
      else
      CALL GETREC(-1,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,SCR1)
      CALL GETREC(-1, 'JOBARC', 'FOCKB   ', NBAS*NBAS*IINTFP,SCRN)
      call xgemm('N', 'N', nbas, nbas, nbas, one, scrn, nbas,
     $     scr1, nbas, zilch, scr2, nbas)
      call xgemm('T', 'N', nbas, nbas, nbas, one, scr1, nbas,
     $     scr2, nbas, zilch, scrn, nbas)
            call SCOPY(nbas, scrn, nbas+1, scr2, 1)
      call putrec(20, 'JOBARC', 'SCFEVLB0', NBAS*IINTFP, scr2)
      endif
#ifdef _DEBUG_LVLM
      write(*,*) "The eigenvalues"
      write(6,"(6(1x,F12.6))") (scr2(i),i=1,nbas)
#endif 
cjp this is somewhat dirty trick, but the eigenvalues in this record
cjp are needed just for sorting purposes when dropping in vtran, so
cjp create 'dummy' energies of MOs for the DROP section in transformation,
cjp which sorts according to MO energies
cjp alternative would be to really compute fock in MO basis and take the diagonal
cjp but it is not necessary
cjp
cAP      do i=1,nbas
cAP      EVAL(i)=i
cAP      enddo
cAP      call putrec(-1,'JOBARC','SCFEVLA0',nbas*iintfp,eval)
cjp
C
C Not so fast Jiri. This is not only dirty but also a problem
C for people who want to feed eigenvectors and eigenvalues from a
C separate calculation into a Fock matrix and want to use
C those eigenvalues as they are in vtran. Ajith  Perera, 08/2004.
C Do the reordering of eiegnvectors if necessary.
C
C We need to have a conditional control over whether we need to
C do the reordering or not. There are two orderings. The SCF
C order: The eigenvalues and the corresponding eigenvectors are in
C ascending order within each irrep. Also, we can order the eigen
C values and the vectors in ascending order regardless of the irrep
C (all the occupied and then virtuals). The latter order is useful
C for DROPMO. When we read eigenvalues and vectors from elsewhere,
C we need to know what order they are in. If they are not in proper
C SCF order we need to reorder them with the call to REOREVL. The
C success of REOREVL depends on the JOBARC recored REORDRA(B) which
C is the mapping to SCF order from another. Whoever use the REORDER
C option please keep a log here as I have done for Brueckner.
C Ajith Perera, 07/2005.
C
       REORDER = .FALSE.
       REORDER = (IFLAGS(22).EQ.0)
C
C Brueckner: The rotated orbital generated in rott1 is in SCF
C            order. No need to do any kind of reordering.
C
      CALL GETREC(0,'JOBARC','SCFEVLA0',LENGTH,EVAL)
      IF (LENGTH.GT.0) THEN
         CALL GETREC(20,'JOBARC','SCFEVLA0',NBAS*IINTFP,EVAL)
      ELSE
         DO I = 1, NBAS
            EVAL(I) = DBLE(I)
         ENDDO
         CALL PUTREC(20,'JOBARC','SCFEVLA0',NBAS*IINTFP,EVAL)
      ENDIF

      CALL GETREC(-1,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCR1)
      IF (REORDER) CALL REOREVL(EVAL,SCRN,NBAS,1)
#ifdef _DEBUG_LVLM
      write(*,*) "The reordered eigenvalues"
      write(6,"(6(1x,F12.6))") (eval(i),i=1,nbas)
#endif
      DO 100 I=1,NIRREP
        IOFFEVC=ISQROF(I)
        CALL GETBLK(SCR1,EVEC(IOFFEVC),NBFIRR(I),NBAS,IREPS(I))
  100 CONTINUE
      IF(IUHF.NE.0)THEN
       CALL GETREC(0,'JOBARC','SCFEVLB0',LENGTH,EVAL)
       IF (LENGTH.GT.0) THEN
          CALL GETREC(-1,'JOBARC','SCFEVLB0',NBAS*IINTFP,EVAL(NBAS+1))
       ELSE
          DO I = 1, NBAS
             EVAL(I+NBAS) = DBLE(I)
          ENDDO
          CALL PUTREC(-1,'JOBARC','SCFEVLB0',NBAS*IINTFP,EVAL(NBAS+1))
       ENDIF
cjp the same for beta part
cjp
cAP      do i=1,nbas
cAP      EVAL(nbas+i)=i
cAP      enddo
cAP      call putrec(-1,'JOBARC','SCFEVLB0',nbas*iintfp,eval(nbas+1))
cjp
       IF (REORDER) CALL REOREVL(EVAL(NBAS+1),SCRN,NBAS,1)
       CALL GETREC(-1,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,SCR1)
       DO 200 I=1,NIRREP
         IOFFEVC=LDIM3+ISQROF(I)
         CALL GETBLK(SCR1,EVEC(IOFFEVC),NBFIRR(I),NBAS,IREPS(I))
  200  CONTINUE
      ENDIF
      RETURN
      END
