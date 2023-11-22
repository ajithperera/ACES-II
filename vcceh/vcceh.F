C
C THIS PROGRAM SOLVES FOR SECOND ORDER PROPERTIES USING THE
C COUPLED-CLUSTER EFFECTIVE HAMILTONIAN APPROXIMATION.
C
C PROGRAMMED BY J.F. STANTON, GAINESVILLE, 1993
C
      PROGRAM VCCEH
      IMPLICIT INTEGER (A-Z)

      PARAMETER(MAXPRT=600)

      LOGICAL CIS,EOMCC,FULDIAG,INCORE,READGUES,VPROP
      LOGICAL CCSD, MBPT, PARTEOM, TRANABCD, TRANABCI
      LOGICAL SS, SD, DS, DD, NODAVID
      LOGICAL REALFREQ, DOMORE, FREQ_EXIST
      LOGICAL EXTRAPOL
      LOGICAL iSOTOPES_PRESENT

      CHARACTER*31 STRING1, STRING2, STRING3
      CHARACTER*80 FNAME
      DOUBLE PRECISION ONEP, ONEM, TDAVID, TMULT, TQUAD, FREQ
      DOUBLE PRECISION XFREQ(200)
      DOUBLE PRECISION TENSOR(MAXPRT,MAXPRT)
      INTEGER IMAG_FREQ(200)
C
      COMMON / / ICORE(1)
      COMMON /ISTART/ I0,ICRSIZ
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO /  NOCCO(2),NVRTO(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /MOPOPS/ MOPOP(8)
      COMMON /AOOFST/ INDOCC(8,2),INDVRT(8,2)
      COMMON /RMAT/ R(10000)
      COMMON /METH/CIS,EOMCC,FULDIAG,INCORE,READGUES   
      COMMON /EIGPROB/ISIDE
      COMMON /FILES/LUOUT, MOINTS
      COMMON /PERT1/NOPERT
      COMMON /PERT/NTPERT, NPERT(8), IPERT(8)
      COMMON/EOMINFO/CCSD, MBPT, PARTEOM, NODAVID
      COMMON/DRVHBAR/SS, SD, DS, DD
      COMMON/TIMSUB/TDAVID, TMULT, TQUAD
      COMMON /NEWSPINAD/ USESPINF
      COMMON /SPINF/ T1F, T2FD, T2FX, T1FI, T2FDI, T2FXI      
      COMMON/LISTDENS/LDENS
C   
C Common blocks used in the quadratic term
C 
      COMMON /QADINTI/ INTIMI, INTIAE, INTIME
      COMMON /QADINTG/ INGMNAA, INGMNBB, INGMNAB, INGABAA,
     &                 INGABBB, INGABAB, INGMCAA, INGMCBB,
     &                 INGMCABAB, INGMCBABA, INGMCBAAB,
     &                 INGMCABBA
      COMMON /APERTT1/ IAPRT1AA             
      COMMON /BPERTT1/ IBPRT1AA
      COMMON /APERTT2/ IAPRT2AA1, IAPRT2BB1, IAPRT2AB1, IAPRT2AB2, 
     &                 IAPRT2AB3, IAPRT2AB4, IAPRT2AA2, IAPRT2BB2,
     &                 IAPRT2AB5 
C
      IONE = 1
      ONEP = 1.00D+00
      ONEM = -1.00D+00
C      
      TDAVID = 0.0D0
      TMULT = 0.0D0
      TQUAD = 0.0D0

      CALL CRAPSI(ICORE,IUHF,0)
      CALL SETMET(IUHF,ICORE(I0))
      imultipole = iflags2(162)
c
      MAXCOR = ICRSIZ
#ifdef _DEBUG_LVLM
      CALL CHECKHBAR(ICORE(I0),MAXCOR,IUHF)
#endif 

C Check whether the NMR spin-spin coupling constant calculation has been
C alrady done with thje standard isotopes and the isotopic mass file exist.
C If these are satisfied, then carry out the recalculating spin-spin coupling
C constants.

      CALL GETREC(0,"JOBARC","JJTENSOR",IPRESENT,TENSOR)
      CALL GFNAME("iSOTOPES",FNAME,LENTH)
      INQUIRE(FILE=FNAME(1:LENTH),EXIST=iSOTOPES_PRESENT)
      IUNITIS=1
      IF (iSOTOPES_PRESENT) THEN
         OPEN(UNIT=IUNITIS,FILE=FNAME(1:LENTH),FORM="FORMATTED")
      ENDIF
      Print*, Ipresent, iSOTOPES_PRESENT
      IF (IPRESENT .GT. 0 .AND. iSOTOPES_PRESENT) THEN
         CALL GETREC(20,'JOBARC','NREALATM',IONE,NATOMS)
         IEOMPROP = IFLAGS2(106)
         WRITE(6,"(2a)") " NMR spin-spin coupling constants of the",
     &                   " isotopic substituted molecule is"
         WRITE(6,"(a)")  " obtained."
         CALL DO_ISOTOPES_JJ(IUNITIS,TENSOR,IPRESENT,NATOMS,IEOMPROP,
     &                       IUNITIS)
         Go TO 911
      ENDIF
C
C Originally vcceh module was set up to run with the density code.
C However, continuing changes in density code encourage us to make the 
C vcceh code independent of the density code. Ajith 09/96.
C
c create a listnumber for original 160 list
c
      LDENS = 260

C During deveoping codes to do amplitude extrapolation this call 
C is moved (see around line 193). This should not affect anything!
CSSS      CALL GENRSPNSDEN(ICORE(I0), MAXCOR, IUHF)
C
      CALL ORBANAL(ICORE(I0),MAXCOR,IUHF)
C     
      CALL INCOR(I0,IUHF)
      MAXCOR = ICRSIZ
C
      IF (IUHF.EQ.0) CALL SPINAD56(ICORE(I0),MAXCOR,IUHF)
C
C For RHF calculations Hbar(Ab,Cd) integrals are not created in
C the lambda code. Because of the large number, these integrals
C are not stored in the disk unless it is required. The quadratic
C term requires these integrals and the following subroutine 
C call and it's dependencies create these integrals and store them
C in the list 233. Ajith 09/21/94
C
      TRANABCD = (ISYTYP(1, 233) .NE. 5
     &   .AND. IFLAGS(93) .NE. 2 .AND. IFLAGS2(122) .EQ. 1
     &   .AND. IFLAGS2(106) .GE. 2) 
C
      TRANABCI = ( CCSD .AND. IFLAGS2(123) .EQ. 1
     &   .AND. (ISYTYP(1,233) .NE. 5) .AND. IFLAGS(93) .NE. 2
     &   .AND. IFLAGS2(106) .GE. 2) 
C
      CALL GETREC(0,"JOBARC","LAMBDA  ",ILAMB,JUNK)
      CALL GETREC(0,"JOBARC","HBAR    ",IHBAR,JUNK)
      CALL GETREC(0,"JOBARC","ETA_ITER",IEXTR,JUNK)
     
      IF (CCSD .AND. IFLAGS2(123) .EQ. 1) THEN
        write(6,*) 'ABCD contribution to HBARABCI is not included !'
      ENDIF
C     
      IF (ILAMB .LT. 0 .AND. IHBAR .EQ. 1 .AND. IEXTR .EQ. 1) THEN
          IF (IFLAGS(93) .NE. 2) THEN
          Write(6,"(3a)") " This is mix PEOM/EOM property calculation",
     &                    " HBARABCI and HBARABCD are"
          Write(6,"(a)")  " augmented with the missing contributions."
          Write(6,*) 
          TRANABCD = .TRUE.
          TRANABCI = .TRUE.
          CALL MODHBAR(ICORE(I0), MAXCOR, IUHF, TRANABCD, TRANABCI)
          ENDIF 
      ENDIF 
C
      IF (IFLAGS2(106) .GE. 2 .AND. IFLAGS2(122) .NE. 2) THEN
        WRITE(6,*) ' QUADRATIC PROPERTIES REQUIRE FULL ABCD INTEGRALS '
CSSS        CALL ERREX
      ENDIF
C
      IF (IFLAGS(1) .GE. 50) THEN
         WRITE(6,*) ' Integral [Ab|Cd] list No. 233 '
         WRITE(6,*)
         CALL PRTLIST(ICORE(I0), MAXCOR, 3, 233, 1, LUOUT, .FALSE.,
     &                DIRPRD, IRPDPD, ISYTYP, NIRREP)
         WRITE(6,*)
      ENDIF
C
      CALL MODF(ICORE(I0),MAXCOR,IUHF,1) 
C
C  IF (IFLAGS(18) .EQ. 13 LOOP OVER DIFFERENT NMR SPIN COUPLING CONSTANTS
C
      IF (IFLAGS(18) .EQ. 13) THEN
        NCASE = 3
      ELSE
        NCASE = 1
      ENDIF
C
      NFREQ = 0
      do i = 1, 200
         imag_freq(i) = 0
      enddo
c
      IF (IFLAGS(18) .EQ. 11) THEN
C
        INQUIRE(FILE='frequency',EXIST=FREQ_EXIST)
        IF (FREQ_EXIST) THEN
            OPEN(UNIT=10, FILE='frequency', FORM='FORMATTED',
     &           STATUS='OLD')
            IJUNK = 0
            READ(10, *, END=800) FREQ
            READ(10, *, END=999) IJUNK
  999       REALFREQ = IJUNK .EQ. 0 .or. abs(freq) .lt. 1.0d-6
            NFREQ = NFREQ + 1
            XFREQ(NFREQ) = FREQ
            if (.not. realfreq) imag_freq(nfreq) = 1
        ELSE
            FREQ     = 0.0D0
            REALFREQ = .TRUE.
        ENDIF
      ENDIF
c
      if ((freq_exist) .and. iflags(18) .eq. 11) then
         ieomprop = iflags2(106)
         if (ieomprop .eq. 0 .or. ieomprop .eq. 3) then
            open(80, file='alpha_cilike', status='unknown')
            rewind(80)
            write(80,801)
         endif
         if (ieomprop .eq. 1 .or. ieomprop .eq. 3) then
            open(81, file='alpha_linear', status='unknown')
            rewind(81)
            write(81,801)
         endif
         if (ieomprop .eq. 2 .or. ieomprop .eq. 3) then
            open(82, file='alpha_quadratic', status='unknown')
            rewind(82)
            write(82,801)
         endif
      endif
c
 801     format(10x, 'Freq',4x, 8x, 'XX',14x, 'YY',14x, 'ZZ', 14x,
     $     'XY', 14x,'XZ',14x,'YZ')
C
  100 CONTINUE

      CALL GETREC(0,"JOBARC","ETA_ITER",ILENGTH,IJUNK)
      I_XP = 0
      IF (ILENGTH .GT. 0) THEN
          EXTRAPOL = .TRUE.
          CALL GETREC(20,"JOBARC","ETA_ITER",IONE,I_XP)
          IF (I_XP .GT. IONE)  THEN
             IF (IFLAGS(93) .EQ. 2) CALL ACES_IO_REMOVE(52,'MOABCD')
             CALL ACES_IO_REMOVE(53,'DERINT')
             CALL ACES_IO_REMOVE(54,'DERGAM')
          ENDIF
      ENDIF 
      CALL GENRSPNSDEN(ICORE(I0), MAXCOR, IUHF)
C
      DO 10 ICASE = 1 , NCASE
        IF (NCASE .EQ. 3) IFLAGS(18) = 7 + ICASE
C
      CALL MAKLST(ICORE(I0),MAXCOR,IUHF)
C
C If this is NMR spin-spin coupling constant calculation
C then we need symmetry adapted NMR perturbation integrals and
C also we need to force the program to go into the SECPRP2.
C
      IF (IFLAGS(18) .EQ. 8 .OR. IFLAGS(18) .EQ. 9 .OR. 
     &    IFLAGS(18) .EQ. 10) CALL SYMADAPT(ICORE(I0), MAXCOR, IUHF)
C
      CALL GETREC(-1, 'JOBARC', 'NTOTPERT', 1, JUNK)
C
      VPROP = JUNK .EQ. 0
C
      IF (VPROP) THEN 
         CALL SECPRP(ICORE(I0), MAXCOR, IUHF, FREQ, REALFREQ,
     $        freq_exist)
      ELSE   
         CALL SECPRP2(ICORE(I0), MAXCOR, IUHF)
      ENDIF
C     
C Calculate diamagnetic spin-orbit contribution to the NMR spin-spin 
C coupling constant based on response-density. Ajith 09/21/94
C

      IF (IFLAGS(18) .EQ. 8) CALL CMPDSO(ICORE(I0), MAXCOR, IUHF)
C
      IF (ICASE .NE. NCASE) THEN
        IF (IFLAGS(93) .EQ. 2) CALL ACES_IO_REMOVE(52,'GAMLAM')
        CALL ACES_IO_REMOVE(53,'DERINT')
        CALL ACES_IO_REMOVE(54,'DERGAM')
      ENDIF
C
      IF (IFLAGS(1) .GE. 1) THEN
C
         WRITE(6,*) 
         STRING1 = 'Total time spent in davidson  :'
         STRING2 = 'Total time spent in hbarxc    :'
         STRING3 = 'Total time spent in quadratic :'
         WRITE(6, 99) STRING1, TDAVID
         WRITE(6, 99) STRING2, TMULT
         IF (IFLAGS2(106) .GE. 2) THEN
            WRITE(6, 99) STRING3, TQUAD
         ENDIF
         WRITE(6,*)
      ENDIF
C
 99   FORMAT(2X, A31, F12.4, ' Seconds')
C
   10  CONTINUE
C
       IF (IFLAGS(18) .EQ. 11 .AND. FREQ_EXIST) THEN
C
        DOMORE = .TRUE.
        IJUNK = 0
        READ(10, *, END=800) FREQ
        READ(10, *, END=899) IJUNK
  899   REALFREQ = IJUNK .EQ. 0 .or. abs(freq) .lt. 1.0d-6
C
      ELSE
        DOMORE = .FALSE.
      ENDIF
C
      IF (DOMORE) THEN
        IF (IUHF .EQ. 0 .AND. IFLAGS2(106)  .GE. 2)  THEN
          CALL SPINAD56(ICORE(I0),MAXCOR,IUHF)
        ENDIF
        CALL ACES_IO_REMOVE(53,'DERINT')
        CALL ACES_IO_REMOVE(54,'DERGAM')
        NFREQ = NFREQ + 1
        XFREQ(NFREQ) = FREQ
        if (.not. realfreq) imag_freq(nfreq) = 1
        GOTO 100
      ENDIF     
C
  800 CONTINUE
C
      write(6,*) ' @vcceh calculation of properties completed '
c
      if ((freq_exist) .and. iflags(18) .eq. 11) then
         if (ieomprop .eq. 0 .or. ieomprop .eq. 3) then
            call write_alpha(80, imultipole, xfreq, imag_freq,
     $           nfreq, nopert, icore(i0))
            close(80, status='KEEP')
         endif
         if (ieomprop .eq. 1 .or. ieomprop .eq. 3) then
            call write_alpha(81, imultipole, xfreq, imag_freq,
     $           nfreq, nopert, icore(i0))
            close(81, status='KEEP')
         endif
         if (ieomprop .eq. 2 .or. ieomprop .eq. 3) then
            call write_alpha(82, imultipole, xfreq, imag_freq,
     $           nfreq, nopert, icore(i0))
            close(82, status='KEEP')
         endif
      endif
c
      CALL MODF(ICORE(I0), MAXCOR, IUHF, -1)
      CALL ACES_AUXCACHE_FLUSH
      CLOSE(UNIT=10, STATUS='KEEP')
C Undo few things
C
#ifdef _DEBUG_LVLM
      CALL CHECKHBAR(ICORE(I0),MAXCOR,IUHF)
#endif 
      IF (.NOT. EXTRAPOL) THEN
         CALL MODIAJK(ICORE(I0), MAXCOR, IUHF, -1.0D0)
         CALL MODAIBC(ICORE(I0), MAXCOR, IUHF, -1.0D0)
         CALL RESET(ICORE(I0), MAXCOR, IUHF)

C When AO LADDERS are used, this does not work and not needed. 
C If needed, this need to be audited. Ajith Perera, 05/2020

      IF (IFLAGS(93) .LT. 2) THEN
          CALL FORMW1(ICORE(I0), MAXCOR, IUHF, .FALSE.)
          CALL HBARABCD(ICORE(I0), MAXCOR, IUHF, .FALSE.)
      ENDIF 
      ENDIF 
c      CALL ZEROLIST(ICORE(I0), MAXCOR, 160)

  911 CONTINUE 

      CALL ACES_FIN
      STOP
      END

