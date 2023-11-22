      SUBROUTINE EVLOUT(EVAL,IRERDR,IVMLSYM,
     &                  IREPS,NBAS,NIRREP,IUHF,NBFIRR)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "maxbasfn.par"
      INTEGER NBAS
      CHARACTER*8 IRRSYM1(MAXBASFN),IRRSYM2(MAXBASFN)
      CHARACTER*4 IRPCHAR
      CHARACTER*5 SPIN(2)
      LOGICAL FANDUP, SYMINFO_EXIST
      DIMENSION EVAL((IUHF+1)*NBAS),IRERDR(NBAS)
      DIMENSION IVMLSYM(NBAS),IREPS(9),NBFIRR(8)
      DIMENSION IRPCHAR(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /POPUL/ NOCC(8,2)
C
      DATA FACTOR /27.2113957D0/
      DATA SPIN /'ALPHA',' BETA'/
cSSSC
cSSSC  If this is a run which is done in spherical orbitals, the orbital
cSSSC  symmetry analysis will be screwed up.  Let the user know this
cSSSC  about the stuff which follows in this routine.
cSSSC
cSSS      IF(IFLAGS(62).EQ.1) THEN
cSSS        WRITE(LUOUT,5100)
cSSS 5100   FORMAT(T3,'@EVLOUT-W,',/,/,
cSSS     &         T10,23('*'),'   WARNING   ',24('*'),/,
cSSS     &         T10,'*',58X,'*',/,
cSSS     &         T10,'*   This calculation was done using spherical ',
cSSS     &            'functions.   *',/,
cSSS     &         T10,'*   Symmetry analysis for these types of ',
cSSS     &            'functions has     *',/,
cSSS     &         T10,'*   not been implemented.  The following ',
cSSS     &            'analysis of       *',/,
cSSS     &         T10,'*   symmetry may be incorrect.           ',
cSSS     &            '                  *',/,
cSSS     &         T10,'*',58X,'*',/,
cSSS     &         T10,60('*'),/)
cSSS      ENDIF
C
C  Now go get the eigenvalues.
C
      CALL GETREC(20,'JOBARC','SCFEVLA0',NBAS*IINTFP,EVAL(1))
      IF(IUHF.EQ.1) THEN
        CALL GETREC(20,'JOBARC','SCFEVLB0',NBAS*IINTFP,EVAL(NBAS+1))
      ENDIF

C
      DO 100 ISPIN=1,(IUHF+1)
        if(ispin.eq.1) then
          CALL GETREC(20,'JOBARC','EVCSYMAF',NBAS*IINTFP,IRRSYM1)
          CALL GETREC(20,'JOBARC','EVCSYMAC',NBAS*IINTFP,IRRSYM2)
        else
          CALL GETREC(20,'JOBARC','EVCSYMBF',NBAS*IINTFP,IRRSYM1)
          CALL GETREC(20,'JOBARC','EVCSYMBC',NBAS*IINTFP,IRRSYM2)
        endif
        DO 101 I=1,NBAS
          IRERDR(I)=I
  101   CONTINUE
C
C  Check to see if we have MOs which are constitued solely from f
C  or higher functions.  Write out a message stating this fact so
C  the user doesn't get worried.
C
        IF(ISPIN.EQ.1) THEN
          FANDUP=.FALSE.
          DO 102 I=1,NBAS
            IF(IRRSYM1(I).EQ.'----') FANDUP=.TRUE.
  102     CONTINUE
          IF(FANDUP) THEN
            WRITE(LUOUT,5500)
 5500       FORMAT(T3,'@EVLOUT-I, There are MOs which are solely made ',
     &                'up of f and/or higher',/,
     &             T14,'angular momentum AOs.  Analysis for these ',
     &                 'cases is not yet',/,
     &             T14,'implemented.  These orbitals will have "----"',
     &                 ' as the symmetry.',/)
          ENDIF
        ENDIF
C
        DO 110 J=1,NBAS
          EMIN=EVAL((ISPIN-1)*NBAS+J)
          IMIN=J
          DO 115 K=J,NBAS
            IF(EVAL((ISPIN-1)*NBAS+K).LT.EMIN) THEN
              EMIN=EVAL((ISPIN-1)*NBAS+K)
              IMIN=K
            ENDIF
  115     CONTINUE
          SJUNK=EVAL((ISPIN-1)*NBAS+J)
          EVAL((ISPIN-1)*NBAS+J)=EVAL((ISPIN-1)*NBAS+IMIN)
          EVAL((ISPIN-1)*NBAS+IMIN)=SJUNK
          IJUNK=IRERDR(J)
          IRERDR(J)=IRERDR(IMIN)
          IRERDR(IMIN)=IJUNK
  110   CONTINUE

        CALL PUTREC(20,"JOBARC","IRPORDMO",NBAS,IRERDR)
C
        DO 200 I=1,NBAS
          DO 201 J=1,NIRREP
            IF(IRERDR(I).GE.IREPS(J).AND.IRERDR(I).LT.IREPS(J+1))
     &        IVMLSYM(I)=J
  201     CONTINUE
  200   CONTINUE
C
        CALL PUTREC(20,'JOBARC','IVMLSYM',NBAS*IINTFP,IVMLSYM)
C
        ITOTOCC=0
        DO 210 I=1,NIRREP
          ITOTOCC=ITOTOCC+NOCC(I,ISPIN)
  210   CONTINUE
C
        IF(IFLAGS(77).NE.0)THEN
          WRITE(LUOUT,5000)SPIN(ISPIN)
 5000     FORMAT(T3,'ORBITALS EIGENVALUES OF THE PARENT SATATE (',A5,')',
     &              '  (1H = 27.2113957 eV)')
          WRITE(6,"(A,A)")" CAUTION: THESE ARE NOT THE DIAGNOL",
     &                    " FOCK OF THE REFERENCE STATE"
          WRITE(6,*)
        ELSE
          WRITE(LUOUT,5002)SPIN(ISPIN)
 5002     FORMAT(T3,'ORBITAL EIGENVALUES (',A5,')',
     &              '  (1H = 27.2113957 eV)',/)
        ENDIF
        WRITE(LUOUT,5010)
 5010   FORMAT(T8,'MO #',8X,'E(hartree)',15X,'E(eV)',11X,'FULLSYM',4X,
     &            'COMPSYM')
        WRITE(LUOUT,5015)
 5015   FORMAT(T8,'----',3X,'--------------------',3X,
     &            '--------------------',3X,'-------',3X,'---------')
        DO 125 J=1,NBAS
          WRITE(LUOUT,5020)J,IRERDR(J),EVAL((ISPIN-1)*NBAS+J),
     &                     EVAL((ISPIN-1)*NBAS+J)*FACTOR,
     &                     IRRSYM1(J)(1:4),IRRSYM2(J)(1:4),
     &                     IVMLSYM(J)
 5020     FORMAT(T3,I3,3X,I3,3X,F20.10,3X,F20.10,4X,A4,6X,A4,
     &           1X,'(',I1,')')
          IF(J.EQ.ITOTOCC) WRITE(LUOUT,5021)
 5021     FORMAT(T3,77('+'))
  125   CONTINUE
c
cmn
c Determine mapping between comp symmetry labels and internal number of irrep.
c
        do j=1,8
           irpchar(j)='XXXX'
        enddo
        do j= 1, nbas
           k = ivmlsym(j)
           if (irpchar(k) .eq. 'XXXX') then
              do i = 1, 4
                 irpchar(k)(i:i) = IRRSYM2(J)(i:i)
              enddo
           endif
        enddo
c
c write this info to file SYMINFO
c
        inquire(file='SYMINFO', exist=syminfo_exist)
        if (.not. syminfo_exist) then
        open(unit=10, FILE='SYMINFO',STATUS='NEW')
        do i = 1, 8
           write(10,599) i,  (irpchar(i)(k:k),k=1,4)
        enddo
 599    format(i4,2x,4a1)
 589    format(i4,a4)
        close(unit=10, status='KEEP')
        endif
cmn end
C  
        Call Map_orb(Nbas, Nirrep, Nocc(1,Ispin), Nbfirr, Ivmlsym, 
     &               Irerdr, Ispin) 

        WRITE(LUOUT,5025)
 5025   FORMAT(/)
  100 CONTINUE        
C      MO #        E(hartree)               E(eV)           FULLSYM    COMPSYM
C      ----   --------------------   --------------------   -------   ---------
C XXX   XXX   000000000.0000000000   000000000.0000000000     XXX       XXX (X)
C
      RETURN
      END
