      SUBROUTINE QRHFIT2(ONEH,DENS,FOCK,EVAL,EVEC,ICORE,MAXMEM,LDIM1,
     &                  LDIM2,NBAS,MXLSSZ1,REPULS,ETOT,IUHF,IOS)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "maxbasfn.par"
C
      DIMENSION ONEH(LDIM1),DENS((IUHF+1)*LDIM1),FOCK((IUHF+1)*LDIM1)
      DIMENSION EVEC((IUHF+1)*LDIM2),EVAL((IUHF+1)*NBAS),
     &          DOCC(MAXBASFN*2)
      DIMENSION ICORE(1)
      DIMENSION IDUMMY(MAXBASFN),IDUMMY2(MAXBASFN),ILOCATE(MAXBASFN)
      DIMENSION IDUMMY3(MAXBASFN)
      DIMENSION IADD(8,2),IREM(8,2), NSUM(8, 2)
      LOGICAL CLOSED
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /POPUL/ NOCC(8,2)
#include "symm2.com"
C
      DATA ONE /1.0/
      DATA TWO /2.0/
C
      INDX2(I,J,N)=I+(J-1)*N
C
      IONE=1
      CALL IZERO(IREM,16)
      CALL IZERO(IADD,16)
C
C SET UP ARRAY WHICH MATCHES EIGENVALUE POSITION WITH ORIGINAL
C POSITION IN ARRAY.
C
      DO 1000 I=1,NBAS
       ILOCATE(I)=I
1000  CONTINUE
C
C  Make the adjustment in the occupation vector based on the QRHF
C  flags.
C
      CALL GETREC(20,'JOBARC','QRHFTOT ',IONE,NMODIFY)
#ifdef _ASSERT
cYAU - This is probably unnecessary since one cannot drop more orbitals
c      than actually exist.
      if (nmodify.gt.maxbasfn) then
         print *, '@QRHFIT2: Assertion failed.'
         print *, '          nmodify = ',nmodify
         print *, '          allowed = ',maxbasfn
         call errex
      end if
#endif
      CALL GETREC(20,'JOBARC','QRHFIRR ',NMODIFY,IDUMMY)
      CALL GETREC(-1,'JOBARC','QRHFLOC ',NMODIFY,IDUMMY2)
      CALL GETREC(-1,'JOBARC','QRHFSPN ',NMODIFY,IDUMMY3)
      IF(IFLAGS(77).NE.0) THEN
       DO 401 I=1,NMODIFY 
        IRRP=IDUMMY(I)
        ILOC=IDUMMY2(I)
        ISPN=IDUMMY3(I)
        IF(IRRP.LT.0)THEN 
         IF(ISPN.NE.1) ISPN=2
         IRRP=-IRRP
         IPOSORG=IRPOFF(IRRP)+NOCC(IRRP,ISPN)+1-MAX(ILOC,1)+
     &                IREM(IRRP,ISPN)-IADD(IRRP,ISPN)
         IPOSABS=ILOCATE(IPOSORG)
         IREM(IRRP,ISPN)=IREM(IRRP,ISPN)+1
         CALL QRHFPOP(IRRP,ILOC,NBAS,EVAL,EVEC,IPOSABS,LDIM2,
     &                ILOCATE,ICORE,ISPN,IOS)
        ELSE
         IF(ISPN.NE.2) ISPN=1
         IPOSORG=IRPOFF(IRRP)+NOCC(IRRP,ISPN)+MAX(ILOC,1)-
     &                IADD(IRRP,ISPN)+IREM(IRRP,ISPN)
         IPOSABS=ILOCATE(IPOSORG)
         IADD(IRRP,ISPN)=IADD(IRRP,ISPN)+1
         CALL QRHFADD(IRRP,ILOC,NBAS,EVAL,EVEC,IPOSABS,LDIM2,
     &                ILOCATE,ICORE,ISPN,IOS)
        ENDIF
401    CONTINUE
C
      ENDIF
C
      CALL PUTREC(20,'JOBARC','QRHFEVAL',NBAS*IINTFP,EVAL)
C
C  Now construct a new density matrix.
C
#ifdef _DPOCC_QRHF
        CALL MKDP_OCCNOS(DOCC, NOCC, NSUM, NBAS, IUHF, .TRUE.)
#endif
      I000=1
      I010=I000+MXIRR2*IINTFP
      I015=I010+MXIRR2*IINTFP
      IF(I010-I000.GT.MAXMEM) THEN
        CALL NOMEM('Construct density','{MKDEN} <-- QRHFIT2 <-- '//
     &             'OOFCOR <-- VSCF',I010-I000,MAXMEM)
      ENDIF
C
      DO 100 ISPIN=1,(IUHF+1)
        IROFF = (ISPIN - 1)*NBAS
        DO 101 I=1,NIRREP
          IF(NBFIRR(I).EQ.0) GOTO 101
#ifdef _DPOCC_QRHF
          IOFF = IROFF + 1
          CALL MKDDEN_4IRREP(EVEC((ISPIN-1)*LDIM2+ISQROF(I)),
     &                       ICORE(I000),ICORE(I010), NSUM(I,ISPIN),
     &                       DOCC(IOFF),NBFIRR(I),ISPIN,IUHF)
          IROFF = IROFF + NBFIRR(I)
#else
          CALL MKDEN(EVEC((ISPIN-1)*LDIM2+ISQROF(I)),ICORE(I000),
     &               NOCC(I,ISPIN),NBFIRR(I),IUHF)
#endif
          CALL SQUEZ2(ICORE(I000),DENS((ISPIN-1)*LDIM1+ITRIOF(I)),
     &                NBFIRR(I))
  101   CONTINUE
  100 CONTINUE
C
C  Now construct the new Fock matrix.
C
      I010=I000+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP
      I020=I010+MXIRR2*IINTFP
      I030=I020+NBAS*NBAS*IINTFP
      IMEM=MAXMEM-(I030-I000)
      IF(I030-I000.GT.MAXMEM) THEN
        CALL NOMEM('Make Fock matrix','{MKFOCK2} <-- QRHFIT2 <-- '//
     &             'OOFCOR <-- VSCF',I030-I000,MAXMEM)
      ENDIF
C
      CALL MKFOCK2(ONEH,FOCK,DENS,ICORE(I000),ICORE(I010),
     &            ICORE(I020),ICORE(I030),ITRILN(NIRREP+1),
     &            MXIRR2,IMEM,NBAS,IUHF)
C
C  Now determine the QRHF energy.  The value for TOL is irrelevant, since
C  it is used to determine whether to use a different format statement for
C  the energy.  This section in MKENER is skipped for a QRHF case.
C
      TOL=0.0
      DMAX=0.0
      ITER=0
      I010=I000+ITRILN(NIRREP+1)*IINTFP
      I020=I010+MXIRR2*IINTFP
      I030=I020+MXIRR2*IINTFP
      IF(I030-I000.GT.MAXMEM) THEN
        CALL NOMEM('QRHF energy','{MKENER} <-- QRHFIT2 <-- OOFCOR'//
     &             ' <--- VSCF',I030-I000,MAXMEM)
      ENDIF
C
      CALL MKENER(ONEH,DENS,FOCK,ICORE(I000),ICORE(I010),ICORE(I020),
     &            ITRILN(NIRREP+1),MXIRR2,REPULS,DMAX,ITER,IUHF,ETOT,
     &            1,TOL)
C
      WRITE(LUOUT,5000)ETOT
 5000 FORMAT(T3,'The QRHF reference energy is ',F20.10,/)
C
C  At this point, the values needed for the QRHF calculation have
C  been placed in the appropriate arrays and will be dumped out
C  with everything else in DMPJOB.
C
C  We also need to calculate the new S**2 value, based on the new
C  occupancies.  This reference function is an eigenfunction of
C  spin, so we only determine the uncontaminated value.
C
      NALPHA=0
      NBETA=0
      DO 200 I=1,NIRREP
        NALPHA=NALPHA+NOCC(I,1)
        NBETA=NBETA+NOCC(I,2)
  200 CONTINUE
      SAVG=ABS(NALPHA-NBETA)/TWO
      S2=SAVG*(SAVG+ONE)
      AMULT=SQRT(1.0+4.0*S2)
      WRITE(LUOUT,5100)AMULT,S2
 5100 FORMAT(/,T3,'     The QRHF average multiplicity is ',F12.7,/,
     &         T3,'The QRHF expectation value of S**2 is ',F12.7,/)
      CALL PUTREC(20,'JOBARC','S2SCF   ',IINTFP,S2)
C
c mn print some further information
C
      write(6,*)
      write(6,*) ' Final QRHF population'
      WRITE(LUOUT,5302)(NOCC(I, 1),I=1,NIRREP)
      WRITE(LUOUT,5303)(NOCC(I, 2),I=1,NIRREP)
 5302 FORMAT(T8,'   Alpha population by irrep: ',8(I3,2X))
 5303 FORMAT(T8,'    Beta population by irrep: ',8(I3,2X))
      write(6,*)
C
C CHANGE TO CLOSED SHELL CALCULATION IF THE STATE IS REOCCUPIED TO
C CLOSED SHELL, AND SCF-REF WAS RHF OR ROHF  AND MAKERHF=ON
C I Think  better logic is to check whether it is actually a
C closed shell and MAKERHF=ON, 03/2011, Ajith Perera.

CSSS      IF ((IFLAGS(11) .EQ. 0 .OR. IFLAGS(11) .EQ. 2)
CSSC     $     .AND. IFLAGS2(140) .eq. 1) THEN
         Nalpha = 0
         Nbeta  = 0
         CLOSED = .FALSE.
         DO IRREP = 1, NIRREP
            Nalpha = Nalpha + NOCC(IRREP,1)
            Nbeta  = Nbeta  + NOCC(IRREP,2)
         ENDDO
         If (Nalpha .EQ. Nbeta) CLOSED=.TRUE.

         IF (CLOSED .AND. IFLAGS2(140) .eq. 1) THEN
            WRITE(6,*) ' The QRHF reference state is closed shell'
            WRITE(6,*)
            IUHF = 0
            iflags(11) = iuhf
            nflags=100
            call putrec(20, 'JOBARC','IFLAGS  ',nflags,iflags)
         endif
CSSS      endif

      RETURN
      END
