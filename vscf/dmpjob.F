#include <flags.h>
      SUBROUTINE DMPJOB(ONEH,  SHALF, DENS,  FOCK,
     &                  EVAL,  EVEC,  SCR1,  SCR2,
     &                  ISCR,  SCR3,  SCRAB, 
     &                  DENSAO, TMO2AO, TMO2AOSCR, 
     &                  LDIM1,
     &                  LDIM2, LDIM3, NBAS, NBASX, 
     &                  ETOT, IUHF,  SEWARD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      LOGICAL GSSOPT,GSSALW,GSSALT,GSSLOK,GSSRED,GSSWRT,GSSUFR,
     &        SEWARD,CLOSED,AODENS_PRESENT,UNOS
      CHARACTER*8 CDENS(2),CEVAL(2),CEVEC(2),CFOCK(2),CPOPO(2),CPOPV(2)
      CHARACTER*8 COAO(2),C1EHMO(2)
      DIMENSION ONEH(LDIM1),SHALF(LDIM1),SCRAB(NBAS*NBAS)
      DIMENSION DENS((IUHF+1)*LDIM1),FOCK((IUHF+1)*LDIM1)
      DIMENSION EVEC((IUHF+1)*LDIM2),EVAL((IUHF+1)*NBAS)
      DIMENSION SCR1(NBAS*NBAS),SCR2(LDIM3),ISCR(NBAS),SCR3(4*LDIM2)
      DIMENSION DENSAO(NBASX*NBASX), TMO2AO(NBASX*NBAS), 
     &          TMO2AOSCR(NBASX*NBASX)
      DIMENSION ITOT1(2),ITOT2(2)
      DIMENSION NOCC(16)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /POPUL/ NOCC
#include "symm2.com"
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /GSCOMA/ GSSOPT,GSSALW,GSSALT,GSSLOK,GSSRED,GSSWRT,GSSUFR
C
      DATA ITWO /2/
      DATA ZILCH /0.0/
      DATA ONE /1.0/
      DATA CDENS /'SCFDENSA','SCFDENSB'/
      DATA CEVAL /'SCFEVLA0','SCFEVLB0'/
      DATA CEVEC /'SCFEVCA0','SCFEVCB0'/
      DATA CFOCK /'FOCKA   ','FOCKB   '/
      DATA CPOPO /'SYMPOPOA','SYMPOPOB'/
      DATA CPOPV /'SYMPOPVA','SYMPOPVB'/
      DATA COAO  /'EVECOAOA','EVECOAOB'/
      DATA C1EHMO /'ONEHMOA ','ONEHMOB '/
C
C     Write converged MOs in SO basis to formatted file (NEWMOS) which 
C     (after renaming to OLDMOS) can be read in subsequent jobs by GETMOS.
C     Now we always write this file.
C
C     IF(GSSWRT)THEN
#ifdef _DEBUG_LVLM
      Write(6,*) "The eigenvectors written to NEMOS file"
      write(6,"(4(1x,f15.7))") (EVEC(i), i=1, (IUHF+1)*LDIM2)
#endif
      CALL PUTMOS(EVEC,LDIM2,IUHF,NIRREP,NBFIRR)
C     ENDIF
C
c mn put orbital eigenvalues from active_smican to jobarc if they
c were defined.
c
      UNOS = .FALSE.
      UNOS = (iflags2(h_IFLAGS2_uno_ref) .NE. 0) 

      if (iflags(16) .eq. 0) then
         nsiz = 0
         call getrec(-1, 'JOBARC', 'NSMICAN ', 1, nsiz)
         if (nsiz .ne. 0) then
            write(6,*) ' Putting active_smican evals to JOBARC ', nsiz
            call getrec(20,'JOBARC','EVALSEMI', nsiz*IINTFP, eval)
            if (.false.) then
               write(6,*) ' Eval read from evalsemi record '
               call output(eval, 1, nsiz, 1, 1, nsiz, 1, 1)
            endif
            do ispin = 1, 1+IUHF
               CALL PUTREC(20,'JOBARC',CEVAL(ISPIN),nsiz*IINTFP,
     &              EVAL)
            enddo
         endif
      endif
C
      DO 10 I=1,2
        ITOT1(I)=0
        ITOT2(I)=0
   10 CONTINUE
      IOFF1=1
      IOFF2=1
      IOFF3=1
      INQUIRE(FILE='AODENS',EXIST=AODENS_PRESENT)
      DO 100 ISPIN=1,(IUHF+1)
        CALL IZERO(ISCR,MAX(2*NIRREP,NBAS))
        DO 110 J=1,NIRREP
          ISCR(J)=NOCC((ISPIN-1)*8+J)
          ISCR(NIRREP+J)=NBFIRR(J)-ISCR(J)
          ITOT1(ISPIN)=ITOT1(ISPIN)+ISCR(J)
          ITOT2(ISPIN)=ITOT2(ISPIN)+ISCR(NIRREP+J)
  110   CONTINUE
        CALL PUTREC(20,'JOBARC',CPOPO(ISPIN),NIRREP,ISCR(1))
        CALL PUTREC(20,'JOBARC',CPOPV(ISPIN),NIRREP,ISCR(NIRREP+1))
C
C IF SEMICANONICAL ORBITALS ARE TO BE USED, THEN MAKE THAT
C  TRANSFORMATION NOW
C
CSSS        IF(IFLAGS(39).EQ.1 .and. (.not. iflags2(148) .eq. 1
CSSS     $    .and. iflags2(150) .eq. 1) )THEN
C
        IF(IFLAGS(39).EQ.1) THEN
         IF(ISPIN.EQ.1) THEN
          WRITE(LUOUT,5000)
 5000     FORMAT(T3,'@DMPJOB-I, Rotating orbitals to semicanonical ',
     &              'basis.',/,
     &           T14,'*** Eigenvalues correspond to new basis. ***',/)
         ENDIF
         I000=1
         I010=I000+LDIM2
         I020=I010+LDIM2
         I030=I020+LDIM2
         I040=I030+LDIM2
         CALL SMICAN(FOCK(IOFF1),SCR3(I000),SCR3(I010),SCR3(I020),
     &               SCR3(I030),SCRAB,EVEC(IOFF2),EVAL(IOFF3),
     &               ISCR,ISCR(NIRREP+1),SCR1,SCR2,NBAS,LDIM1,
     &               LDIM3,NIRREP,ISPIN)
C
        ENDIF
        CALL PUTREC(20,'JOBARC',CEVAL(ISPIN),NBAS*IINTFP,
     &              EVAL((ISPIN-1)*NBAS+1))
#ifdef _DEBUG_LVLM
        Write(6,*)
        Write(6,"(a)") "@-DMPJOB; The eigenvalues"
        Write(6,"(6F12.6)") (Eval(i),i=1, (1+iuhf)*Nbas)
#endif 
        CALL ZERO(SCR1,NBAS*NBAS)
        CALL MKFULL(SCR1,DENS((ISPIN-1)*LDIM1+1),SCR2,LDIM1,LDIM3,
     &              NBAS,1)
        CALL PUTREC(20,'JOBARC',CDENS(ISPIN),NBAS*NBAS*IINTFP,SCR1)
c
#ifdef _DEBUG_LVLM
        if (ispin .eq.1) then
        Write(6,*) "The density matrix in sym alpha AO basis"
        call output(scr1, 1, nbas, 1, nbas, nbas, nbas, 1)
        else
        Write(6,*) "The density matrix in sym beta AO basis"
        call output(scr1(nbas*nbas+1), 1, nbas, 1, nbas, nbas, 
     &  nbas, 1)
        endif 
#endif
C
C Built the non-symmetry adapted AO density matrix.
C
         IF ((Iflags(h_IFLAGS_vib) .GT. 1) .or.
     $        (iFlags2(h_IFLAGS2_prep_vibron).EQ.1) .or.
     $        (iFlags2(h_IFLAGS2_scan_normal).EQ.1)) then

           IF (.NOT. AODENS_PRESENT) THEN
               IF (ISPIN .EQ. 1) THEN
                   OPEN(UNIT=61,FILE='AODENS',STATUS='NEW',
     &              ACCESS='SEQUENTIAL',FORM='FORMATTED')
               ELSE
                   OPEN(UNIT=61,FILE='AODENS',STATUS='OLD',
     &              ACCESS='SEQUENTIAL',FORM='FORMATTED')
                    ITEST = 0
                    DO WHILE (ITEST .EQ. 0) 
                       READ(61, "(A)", IOSTAT=IERR) IJUNK
                       IF (IERR .LT. 0) ITEST = 1
                       IF (IERR .GT. 0) THEN
                          WRITE (*,*) ' Problem reading AODENS file '
                          CALL  ERREX
                       ENDIF
                    ENDDO
                    WRITE(61,*) 
               ENDIF
              CALL GETREC(20,'JOBARC','CMP2ZMAT',NBASX*NBAS*IINTFP,
     &                    TMO2AO)
              CALL XGEMM('N','N',NBASX,NBAS,NBAS,1.0D+00,
     &                    TMO2AO,NBASX,SCR1,NBAS,0.0D+00,TMO2AOSCR,
     &                    NBASX)
              CALL XGEMM('N','T',NBASX,NBASX,NBAS,1.0D+00,
     &                    TMO2AO,NBASX,TMO2AOSCR,NBASX,0.0D+00,
     &                    DENSAO,NBASX)

#ifdef _DEBUG_LVLM
        if (ispin .eq.1) then
        Write(6,*) "The density matrix in no-sym AO basis"
        call output(densao, 1, nbasx, 1, nbasx, nbasx, nbasx, 1)
        call checksum("@-DMPJOB-DENSA", densao, NBAS*NBAS)
         
        else
        Write(6,*) "The density matrix in no-sym beta AO basis"
        call output(densao, 1, nbasx, 1, nbasx, nbasx, nbasx, 1)
        call checksum("@-DMPJOB-DENSB", densao, NBAS*NBAS)
        endif 
#endif
              CALL PUTAODENS(DENSAO, 61, NBASX)
              CLOSE(UNIT=61)

#ifdef _DEBUG_LVLM2
        CALL GETREC(20,'JOBARC','ZMAT2CMP',NBASX*NBAS*IINTFP,TMO2AO)
        CALL XGEMM('N','N',NBAS,NBASX,NBASX,1.0D+00,
     &             TMO2AO,NBAS,DENSAO,NBASX,0.0D+00,TMO2AOSCR,NBAS)
        CALL XGEMM('N','T',NBAS,NBAS,NBASX,1.0D+00,
     &             TMO2AO,NBAS,TMO2AOSCR,NBAS,0.0D+00,DENSAO,NBAS)
        Write(6,*) "The density matrix in MO basis (rebuilt)"
        call output(densao, 1, nbas, 1, nbas, nbas, nbas, 1)
#endif
           ENDIF
        ENDIF
C
        CALL ZERO(SCR1,NBAS*NBAS)
        CALL MKFULL(SCR1,FOCK((ISPIN-1)*LDIM1+1),SCR2,LDIM1,LDIM3,
     &              NBAS,1)
        CALL PUTREC(20,'JOBARC',CFOCK(ISPIN),NBAS*NBAS*IINTFP,SCR1)
        
#ifdef _DEBUG_LVLM
      write(6,*) "The Fock matrix (FOCKA and B) from dmpjob"
      call output(scr1, 1,nbas,1,nbas,nbas, nbas,1)
      call checksum("@DMPJOB;"//CFOCK(ISPIN), scr1, nbas*nbas)
#endif 
C If UNOs are the eigenvalues array does not have the correct data.
C So the printing looks funny. Lets copy the diagonal of the Fock
C matrix to EVAL and put it into the JOBARC here.
C
        IF (UNOS) THEN
            DO I = 1, NBAS
               K = I + (I-1)*NBAS 
               EVAL((ISPIN-1)*NBAS+I) = SCR1(K)
            ENDDO

            CALL PUTREC(20,'JOBARC',CEVAL(ISPIN),NBAS*IINTFP,
     &                  EVAL((ISPIN-1)*NBAS+1))
        ENDIF 
C
        CALL ZERO(SCR1,NBAS*NBAS)
        CALL MKFULL(SCR1,EVEC((ISPIN-1)*LDIM2+1),SCR2,LDIM2,LDIM3,
     &              NBAS,0)
        CALL PUTREC(20,'JOBARC',CEVEC(ISPIN),NBAS*NBAS*IINTFP,SCR1)

#ifdef _DEBUG_LVLM
      write(6,*) "The eigenvectors (SCFEVCA0) from dmpjob"
      call output(scr1, 1,nbas,1,nbas,nbas, nbas,1)
      call checksum("@-DMPJOB", scr1, NBAS*NBAS)
#endif 
C
CJDW 8/26/97. Write SCF SO vectors also to SCFEVCAS and SCFEVCBS. These
C             records should NEVER be reordered. They can then be used as
C             initial guess MOs without reordering.
C
        IF(ISPIN.EQ.1)THEN
         CALL PUTREC(20,'JOBARC','SCFEVCAS',NBAS*NBAS*IINTFP,SCR1)
C
C For RHF MOLCAS runs, also write the redundent beta vectors to JOBARC
C file. Note that this will make Roland's life simpler since he
C doesn't need to check whether it is UHF/ROHF or RHF in his code;
C
         IF (SEWARD) CALL PUTREC(20, 'JOBARC', 'SCFEVCBS', 
     &                           NBAS*NBAS*IINTFP,SCR1)
        ELSE
         CALL PUTREC(20,'JOBARC','SCFEVCBS',NBAS*NBAS*IINTFP,SCR1)
        ENDIF
C
        CALL ZERO(SCR1,NBAS*NBAS)
        CALL MKFULL(SCR1,SHALF,SCR2,LDIM1,LDIM3,NBAS,1)
        CALL PUTREC(20,'JOBARC','SMHALF  ',NBAS*NBAS*IINTFP,SCR1)
C
        CALL ZERO(SCR1,NBAS*NBAS)
        CALL ZERO(SCR2,LDIM3)
        CALL ZERO(SCR3,4*LDIM2)
        DO 170 ICNT=1,NIRREP
          IF(NBFIRR(ICNT).EQ.0) GOTO 170
          CALL EXPND2(SHALF(ITRIOF(ICNT)),SCR2,NBFIRR(ICNT))
          CALL MINV(SCR2,NBFIRR(ICNT),NBFIRR(ICNT),SCR3,DET,ZILCH,0,1)
          CALL XGEMM('T','N',NBFIRR(ICNT),NBFIRR(ICNT),NBFIRR(ICNT),
     &               ONE,SCR2,NBFIRR(ICNT),EVEC((ISPIN-1)*LDIM2+
     &               ISQROF(ICNT)),NBFIRR(ICNT),ZILCH,
     &               SCR3(ISQROF(ICNT)),NBFIRR(ICNT))
  170   CONTINUE
        CALL MKFULL(SCR1,SCR3,SCR2,LDIM1,LDIM3,NBAS,0)
        CALL PUTREC(20,'JOBARC',COAO(ISPIN),NBAS*NBAS*IINTFP,SCR1)
C
C
C mn: write one-electron integrals in AO basis in full triangular format
c Used in mrcc module
c
        call mkfull(scrab, oneh, scr1, ldim1, ldim2, nbas, 1)
        ithru = 0
        do i=1, nbas
           do j= 1, i
              ithru = ithru + 1
              scr1(ithru) = scrab(i+(j-1)*nbas)
           enddo
        enddo
        CALL PUTREC(20,'JOBARC','ONEHAO  ',ithru*IINTFP,scr1)
c
c mn: Write overlap integrals in the same full triangular format
c
        CALL GETREC(20,'JOBARC','AOOVRLAP',nbas*nbas*IINTFP,scrab)
        ithru = 0
        do i=1, nbas
           do j= 1, i
              ithru = ithru + 1
              scr1(ithru) = scrab(i+(j-1)*nbas)
           enddo
        enddo
        CALL PUTREC(20,'JOBARC','SAO_TRI ',ithru*IINTFP,scr1)
cmn end
        CALL ZERO(SCR1,NBAS*NBAS)
        CALL ZERO(SCR2,LDIM3)
        CALL ZERO(SCR3,4*LDIM2)
        CALL ZERO(SCRAB,NBAS*NBAS)
        DO 180 ICNT=1,NIRREP
          IF(NBFIRR(ICNT).EQ.0) GOTO 180
          CALL EXPND2(ONEH(ITRIOF(ICNT)),SCR2,NBFIRR(ICNT))
          CALL XGEMM('T','N',NBFIRR(ICNT),NBFIRR(ICNT),NBFIRR(ICNT),
     &               ONE,EVEC((ISPIN-1)*LDIM2+ISQROF(ICNT)),
     &               NBFIRR(ICNT),SCR2,NBFIRR(ICNT),ZILCH,SCR3,
     &               NBFIRR(ICNT))
          CALL XGEMM('N','N',NBFIRR(ICNT),NBFIRR(ICNT),NBFIRR(ICNT),
     &               ONE,SCR3,NBFIRR(ICNT),EVEC((ISPIN-1)*LDIM2+
     &               ISQROF(ICNT)),NBFIRR(ICNT),ZILCH,
     &               SCR3(LDIM2+ISQROF(ICNT)),NBFIRR(ICNT))
  180   CONTINUE
        CALL MKFULL(SCR1,SCR3(LDIM2+1),SCR2,LDIM1,LDIM3,NBAS,0)
        CALL SQUEZ2(SCR1,SCRAB,NBAS)
        IABTRI=(NBAS*(NBAS+1))/2
        CALL PUTREC(20,'JOBARC',C1EHMO(ISPIN),IABTRI*IINTFP,SCRAB)
C
        IOFF1=IOFF1+LDIM1
        IOFF2=IOFF2+LDIM2
        IOFF3=IOFF3+NBAS 
  100 CONTINUE
C
      CALL PUTREC(20,'JOBARC','NOCCORB ',ITWO,ITOT1)
      CALL PUTREC(20,'JOBARC','NVRTORB ',ITWO,ITOT2)
C
      CALL PUTREC(20,'JOBARC','SCFENEG ',IINTFP,ETOT)
      CALL PUTREC(20,'JOBARC','TOTENERG',IINTFP,ETOT)
      CALL PUTREC(20,'JOBARC','NUMBASIR',NIRREP,NBFIRR)
C
C  Dump out the values for OCCUPYA and OCCUPYB; if this is a QRHF
C  calculation, these will be different.
C
      CALL PUTREC(20,'JOBARC','OCCUPYA ',NIRREP,NOCC(1))
      IF(IUHF.EQ.1) THEN
        CALL PUTREC(20,'JOBARC','OCCUPYB ',NIRREP,NOCC(9))
      ENDIF
C
      ICNT=0
      CALL IZERO(ISCR,NBAS)
      DO 150 I=1,NIRREP
        DO 151 J=1,NBFIRR(I)
          ICNT=ICNT+1
          ISCR(ICNT)=I
  151   CONTINUE
  150 CONTINUE
      CALL PUTREC(20,'JOBARC','IRREPALP',NBAS,ISCR)
      IF(IUHF.EQ.1) CALL PUTREC(20,'JOBARC','IRREPBET',NBAS,ISCR)
C
C  Now calculate spin density, if this is an open-shell calculation.
C
      IF(IUHF.EQ.1) THEN
        CALL SCOPY(LDIM1,DENS,1,SCRAB,1)
        CALL SAXPY(LDIM1,-ONE,DENS(LDIM1+1),1,SCRAB,1)
        CALL ZERO(SCR1,NBAS*NBAS)
        CALL MKFULL(SCR1,SCRAB,SCR2,LDIM1,LDIM3,NBAS,1)
        CALL PUTREC(20,'JOBARC','SCFSPDEN',NBAS*NBAS*IINTFP,SCR1)
      ENDIF
C
C  Now compute occupation number vectors.
C
      FACT=DFLOAT(2-IUHF)
      IOFF1=0
      DO 350 ISPIN=1,1+IUHF
       IOFF=0
       CALL ZERO(SCR1,NBAS)
       DO 351 IRREP=1,NIRREP
        IOCC=NOCC(IOFF1+IRREP)
        NBF=NBFIRR(IRREP)
        DO 352 I=1,IOCC
         SCR1(IOFF+I)=FACT
352     CONTINUE
        IOFF=IOFF+NBF
351    CONTINUE
       IOFF1=8
       IF(ISPIN.EQ.1)CALL PUTREC(20,'JOBARC','OMEGA_A ',
     &                           IINTFP*NBAS,SCR1) 
       IF(ISPIN.EQ.2)CALL PUTREC(20,'JOBARC','OMEGA_B ',
     &                           IINTFP*NBAS,SCR1) 
350   CONTINUE
C
CJDW 10/24/95. Write energy and coordinates in POLYRATE jobs.
C
      IF(IFLAGS2(113).EQ.1)THEN
        IONE = 1
        CALL GETREC(20,'JOBARC','NATOMS  ',IONE,NATOMS)
        CALL POLYPRT0(NATOMS,IINTFP,IFLAGS(1),
     &                .TRUE.,.TRUE.,.FALSE.,.FALSE.)
      ENDIF
C
C CHANGE TO CLOSED SHELL CALCULATION IF THE STATE IS REOCCUPIED TO
C CLOSED SHELL when AND MAKERHF=ON. This applies only for Marcel's
C MRCC code runs. Ajith Perera, 03/2001
C
         CLOSED = .FALSE.
         Nalpha = 0 
         Nbeta  = 0
         DO IRREP = 1, NIRREP
            Nalpha = Nalpha + NOCC(IRREP) 
            Nbeta  = Nbeta  + NOCC(IRREP+8)
         ENDDO
         If (Nalpha .EQ. Nbeta) CLOSED=.TRUE.

         IF (CLOSED .AND. IFLAGS2(140) .eq. 1) THEN
            WRITE(6,*) ' The reference state is closed shell'
            WRITE(6,*)
            IUHF = 0
            iflags(11) = iuhf
            nflags=100
            call putrec(20, 'JOBARC','IFLAGS  ',nflags,iflags)
         endif
C
      RETURN
      END
