
      SUBROUTINE SECPRP(ICORE, MAXCOR, IUHF,FREQ, REALFREQ, freq_exist)
C
      IMPLICIT INTEGER (A-Z)
      LOGICAL IDENTITY,REALFREQ, freq_exist
C
      DOUBLE PRECISION SDOT,X1,X00,MEXPECT,MBAR00,H0QM,HGGM,X
      DOUBLE PRECISION ONE,TWO,ZILCH,FACT,FREQ,TRACE,ZZ,XPROP,
     &   XPROP2, XVEC, FREQ0
      LOGICAL CCSD, MBPT, PARTEOM
      LOGICAL DOLEFT, DORIGHT, NODAVID, INCORE, IMAG_QUAD
C
      CHARACTER*8 LABEL
      CHARACTER*6 SIDES(2)
      CHARACTER*80 STRING
C
      PARAMETER(MAXPRT=100)
C 
      DIMENSION ICORE(MAXCOR),MPERT(8),LENGTH(8)
      DIMENSION ISYMPRT(MAXPRT),IOFFPRT(MAXPRT)
      DIMENSION MEXPECT(MAXPRT),MBAR00(MAXPRT)
      DIMENSION XPROP(MAXPRT,MAXPRT), XPROP2(MAXPRT,MAXPRT),
     $     XVEC(MAXPRT*MAXPRT)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/FLAGS/IFLAGS(100)
      COMMON/PERT1/NOPERT
      COMMON/PERTTYP/LABEL(100)
      COMMON/FLAGS2/ IFLAGS2(500)
      COMMON/FILES/LUOUT,MOINTS
      COMMON/EOMINFO/CCSD, MBPT, PARTEOM, NODAVID
C
      DATA ONE,TWO,ZILCH /1.0D0,2.0D0,0.0D0/
      DATA SIDES /'right.','left. '/
C
      INDXF(I,J,N)=I+(J-1)*N
C
      ISIDE = 1
      IONE  = 1
      IROW  = 1 
      ICOL  = 1
      IOFF  = 1
      FREQ0 = FREQ
C
C Set flags for the unlinked and quadratic term
C
      IEOMPROP = IFLAGS2(106)
      IF (ISIDE .EQ. 2 .AND. IEOMPROP .GE. 2) THEN
        WRITE(6,*) ' SOLVING FROM LEFT HAND SIDE CANNOT BE USED',
     &     ' IN COMBINATION WITH EOMPROP .GE. 2 '
        CALL ERREX
      ENDIF
C
C The IFIRST and ILAST is introduced to do all the possible options in a
C single calculation if requested. When IEOMPROP = 0, 1 and 2 user is requesting
C to do a only one option in a single run: CILIKE or LINEAR or QUADRATIC respectively.
C When IEOMPROP = 3, user is requesting to do all the possible options in a 
C single calculation. 
C
      IF (IEOMPROP .EQ. 0 .OR. IEOMPROP .EQ. 1 .OR.
     &   IEOMPROP .EQ. 2) THEN
C
        IFIRST = IEOMPROP
        ILAST  = IEOMPROP
C
      ELSE IF (IEOMPROP .EQ. 3) THEN
C
        IFIRST = 0
        ILAST = 2 
C
      ENDIF
C     
      NBAS = NOCCO(1)+NVRTO(1)
      CALL GETREC(20,'JOBARC','NBASTOT ',IONE,NBAST)
C     
      KMAX = MAXEXP - 1
C
      I000 = 1
      I010 = I000 + MAX(NBAST,KMAX*KMAX)*IINTFP
      I020 = I010 + KMAX*IINTFP
C
C Memory checks were added by Ajith 07/95
C 
      IF (I020 .GT. MAXCOR) CALL INSMEM('SECPRP', I020, MAXCOR) 
C
      MXCOR = (MAXCOR - I020 + 1)/IINTFP
C
      WRITE(6, *)
      WRITE(6, 5000) SIDES(ISIDE)
 5000 FORMAT(T3,'@SECPRP-I, Linear equations solved from ',A6)
C
C SET UP LISTS TO HOLD PERTURBED T VECTOR AND LEFT-HAND DIPOLE
C ELEMENTS
C
      CALL IZERO(MPERT, 8)
C
      DO 1 IPERT = 1, NOPERT 
C
        CALL GTPRTSYM(ICORE(I000), ICORE(I010), MAXCOR-I010+1,
     &     NBAST, IRREPX, LABEL(IPERT))
        MPERT(IRREPX)  = MPERT(IRREPX) + 1
        IOFFPRT(IPERT) = MPERT(IRREPX)
        ISYMPRT(IPERT) = IRREPX
C
    1 CONTINUE
C      Print*, (IOFFPRT(I), I=1, NOPERT)
C
C Record symmetry of perturbations to perturb density analysis 08/2007
C Ajith perera
C
      CALL PUTREC (20, 'JOBARC', 'NPERTIRR', 8, MPERT)     
C     
      IF(FREQ .NE. ZILCH .AND. REALFREQ) THEN
        ITOP = 2
        FACT = ONE
      ELSE
        ITOP = 1
        FACT = TWO
        IF(.NOT. REALFREQ) FREQ = FREQ*FREQ
      ENDIF
      imag_quad = freq_exist .and. .not. realfreq .and.
     $     (ieomprop .eq. 2 .or. ieomprop .eq. 3)
      write(6,*) ' Include Imaginary quadratic terms : ', imag_quad
c
c this indicates imaginary frequencies plus a quadratic calculation is requested.
C
C
C  create lists for perturbed amplitudes
C
      WRITE(6,*)
C
c      INEWFIL = 1
      call aces_io_remove(53,'DERINT')
      inewfil = 0
      DO 2 IRREP = 1, NIRREP
C
        WRITE(6, 9000) MPERT(IRREP), IRREP
 9000   FORMAT(T3,'@SECPRP-I, There are ',I2,' perturbations ',
     &     ' of symmetry species ',I2,'.')
C
        LENGTH(IRREP)=0
C
        DO 3 ISPIN =1, 1 + IUHF
          LENGTH(IRREP) = LENGTH(IRREP) + IRPDPD(IRREP, 8+ISPIN)
    3   CONTINUE
C
        DO 4 ISPIN = 3, 3-2*IUHF, -1
C
          LSTTYP = 43 + ISPIN
          LENGTH(IRREP) = LENGTH(IRREP) +
     &       IDSYMSZ(IRREP,ISYTYP(1,LSTTYP),
     &       ISYTYP(2,LSTTYP))
    4   CONTINUE
C
C In John stanton's original implementation these list numbers 
C were 473 and 474. They were changed to 300 lists to accommodate
C 400 list's to the quadratic contributions. 06/1994 by Ajith.
C Note that in the first UPDMOI call to the list 373, it opens
C the 'DERINT' file. For the frequancy dependent calculations
C including quadratic contributions we have to save first-order
C perturb amplitudes for both frequencies +W and -W. The list
C number 373 store the amplitudes corresponding to the -W, while
C the list No. 375 store the amplitudes corresponding to the +W
C frequency.
C
        CALL UPDMOI(MPERT(IRREP),LENGTH(IRREP),IRREP,373,INEWFIL,0)
        INEWFIL = 0
        CALL UPDMOI(MPERT(IRREP),LENGTH(IRREP),IRREP,374,INEWFIL,0)
C
        IF (ITOP .EQ. 2 .or. imag_quad) THEN
          CALL UPDMOI(MPERT(IRREP),LENGTH(IRREP),IRREP,375,INEWFIL,0)
        ENDIF
C
    2 CONTINUE
C
      DO 5 ICOUNT = 1, ITOP 
C
        DO 10 IPERT = 1, NOPERT          
C
          WRITE(6, *)
          WRITE(6, 500)IPERT, LABEL(IPERT)
C
          IRREPX  = ISYMPRT(IPERT)
          IOFFSET = IOFFPRT(IPERT)
          NUMPERT = IOFFSET
C
          IF(IRREPX .NE. 0) THEN
C
            DORIGHT = ISIDE. EQ. 1
            DOLEFT = ISIDE .EQ. 2
            CALL NEWLST(IRREPX, ICORE, MAXCOR, IUHF)
            CALL OPRBAR(ICORE, MAXCOR, IUHF, LABEL(IPERT), IRREPX,
     &         IOFFSET, IPERT, MEXPECT(IPERT), 
     &         MBAR00(IPERT), .TRUE.,
     &         DOLEFT, DORIGHT, IEOMPROP)
C
            if (imag_quad) then

c pick up inhomogenous part and store temporarily on 375                                                                
               CALL GETLST(ICORE(I000),IOFFSET,1,1,IRREPX,372+ISIDE)
               CALL PUTLST(ICORE(I000),IOFFSET,1,1,IRREPX,374+ISIDE)
            endif

            WRITE(6, 501)
C
            INCORE = .TRUE.
C
            IF (PARTEOM .AND. NODAVID .AND. REALFREQ) THEN
C
C  SEE IF IN CORE DAVIDSON PROCEDURE IS POSSIBLE
C

              NDIM = IRPDPD(IRREPX, 9)
              N1AA = NDIM
              IF (IUHF .NE. 0) NDIM = NDIM + IRPDPD(IRREPX,10)
              KMAX = MIN(NDIM, MAXEXP-1)
              IC = I020
              IHXC = IC + KMAX * NDIM * IINTFP
              IY1 = IHXC + KMAX * NDIM * IINTFP
              IX1 = IY1 + NDIM* IINTFP
              IH0 = IX1 + NDIM * IINTFP
              ISCR = IH0 + NDIM*IINTFP
              MXCOR = (MAXCOR - ISCR + 1) / IINTFP
              IMAX = ISCR + 2 * NSIZEC * IINTFP
              IF (IMAX .GE. MAXCOR) THEN
                     CALL INSMEM('SECPRP', IMAX, MXCOR)
                INCORE = .FALSE.
              ENDIF
            ELSE
              INCORE = .FALSE.
            ENDIF
C
            IF (INCORE) THEN
              CALL PLINEQY(ICORE(ISCR),MXCOR,IUHF,ICORE(I000),
     &           ICORE(I010),KMAX,IRREPX,
     &           IOFFSET,IPERT,FREQ,ISIDE,REALFREQ,
     &           ICORE(IC), ICORE(IHXC), ICORE(IY1),
     &           ICORE(IX1), ICORE(IH0), N1AA, NDIM)
              
            ELSE
              CALL LINEQY(ICORE(I020), MXCOR, IUHF, ICORE(I000),
     &           ICORE(I010), KMAX, IRREPX, IOFFSET, IPERT,
     &           FREQ, ISIDE, REALFREQ)
C
            ENDIF
          ENDIF 
C     
C Save the first-order perturb T amplitudes corresponding to +W
C frequency in list 375 if it is a frequency dependent calculation
C      
          IF (ICOUNT .EQ. 1 .AND. ITOP .EQ. 2) THEN
C     
            IA = I000
            IB = IA + NSIZEC*IINTFP
C     
            CALL GETLST(ICORE(IA), IOFFSET, 1, 1, IRREPX, 372+ISIDE)
            CALL PUTLST(ICORE(IA), IOFFSET, 1, 1, IRREPX, 375)
C
          ENDIF
c
          if (imag_quad) then
c
c calculate imaginary part of T-response amplitude (1/freq) * (<n|[Hbar, T_real]|0> - <n|Obar|0>)
c
            IA = I000
            IB = IA + NSIZEC*IINTFP
            CALL GETLST(ICORE(I000),IOFFSET,1,1,IRREPX,372+ISIDE)
            CALL UPDATES(IRREPX,ICORE(I000),444,0,IUHF)
            CALL HBARXC(ICORE(I000),MAXCOR,IUHF,ISIDE,IRREPX)
            CALL FETCHVEC(IRREPX,ICORE(I000),IUHF,2,490,460)
c
c pick up Obar quantities, stored on 375
            CALL GETLST(ICORE(IB),IOFFSET,1,1,IRREPX,374+ISIDE)
            call DAXPY(nsizec, -1.0d0, ICORE(IB), 1, ICORE(I000), 1)
c
            call DSCAL(nsizec, 1.0d0/freq0, ICORE(I000), 1)
            CALL PUTLST(ICORE(I000),IOFFSET,1,1,IRREPX,374+ISIDE)
            if (.false.)  then
            write(6,*) ' Imaginary perturbed amps '
            np = min(nsizec, 100)
            call output(ICORE(I000), 1, 1, 1, np, 1, np, 1)
            endif
         endif
C
  500     FORMAT(T3,'@SECPRP-I, Perturbation ',I2,' is ',A8,'.')
  501     FORMAT(T3,'@SECPRP-I, Solving for derivative T amplitudes.')
C
   10   CONTINUE
C
        FREQ = -FREQ 
C
    5 CONTINUE 
C     
C Calculate properties
C
      DO 20 IEOMPROP = IFIRST, ILAST
C
        IF (ISIDE .EQ. 2 .AND. IEOMPROP .EQ. 2) GOTO 20
C
C
C Print the results to the user with appropriate information
C so the that user knows what he is getting out of the machine.
C
        write(6,*)
        if (ccsd) then
          if (parteom) then
            string ='$Partitioned EOM-CCSD second order properties.'
          else
            string ='$Full EOM-CCSD second order properties.'
          endif
        else
          if (parteom) then
            string ='$Partitioned EOM-MBPT second order properties.'
          else
            string ='$Full EOM-MBPT second order properties.'
          endif
        endif
        Write (LUOUT, '(1X,A)') String(2:Len(string)-1)
        if (ieomprop .eq. 0) then
          string='$CI-like approximation to evaluation of properties.'
        else if (ieomprop .eq. 1) then
          string='$Linear approximation to evaluation of properties.'
        else if (ieomprop .eq. 2) then
          string='$Fully quadratic evaluation of properties.'
        endif
        Write (LUOUT, '(1X,A)') String(2:Len(string)-1)
        WRITE(LUOUT,*)
C
C Initialize the target array (XPROP) to zero for the accumulation
C of final results.
C
        CALL ZERO(XPROP, MAXPRT*MAXPRT)
C
        TRACE = ZILCH
C            
        DO 30 IPERT = 1, NOPERT
          IRREPX  = ISYMPRT(IPERT)
          IOFFSET = IOFFPRT(IPERT)
          NUMPERT = IOFFSET
C               
          IF (IRREPX .NE. 0) THEN
            DOLEFT = ISIDE .EQ. 1
            DORIGHT = ISIDE .EQ. 2
            CALL NEWLST(IRREPX, ICORE, MAXCOR, IUHF)
            CALL OPRBAR(ICORE, MAXCOR, IUHF, LABEL(IPERT), IRREPX,
     &         IOFFSET, IPERT, MEXPECT(IPERT),
     &         MBAR00(IPERT), .TRUE.,
     &         DOLEFT, DORIGHT, IEOMPROP)

          ENDIF
C
   30   CONTINUE
C
        DO 60 ICOUNT = 1, ITOP
C
          IF (ITOP .EQ. 2) THEN
            IF (ICOUNT .EQ. 1) THEN
              WRITE(6,*) ' Positive frequency contribution '
            ELSE
              WRITE(6,*) ' Negative frequency contribution '
            ENDIF
          ENDIF
          DO 40 IPERT = 1, NOPERT
C
            IRREPX = ISYMPRT(IPERT)
C
            IF (IRREPX .NE. 0) THEN
C     
              IOFFSET = IOFFPRT(IPERT)
C
              CALL NEWLST(ISYMPRT(IPERT), ICORE, MAXCOR, IUHF)
C
              IA = I000
              IB = IA + NSIZEC*IINTFP
              MXCOR = MAXCOR - IB + 1
C
              LISTIN = 372 + ISIDE
              IF (ICOUNT .EQ. 1 .AND. ITOP .EQ. 2) 
     &           LISTIN = 375
C
              CALL GETLST(ICORE(IA), IOFFSET, 1, 1, IRREPX,
     &           LISTIN)
#ifdef _DEBUG_LVL0
              call checksum("374",ICORE(IA),NSIZEC)
#endif 

              IF (IUHF .EQ. 0) THEN
                CALL SPNTSING(ISYMPRT(IPERT), ICORE(IA),
     &             ICORE(IB), MXCOR)
              ENDIF
C
              DO 50 JPERT = 1, NOPERT
C
                IF (ISYMPRT(IPERT) .EQ. ISYMPRT(JPERT)) THEN
C
                  JOFFSET = IOFFPRT(JPERT)
C     
                       
     
                  CALL GETLST(ICORE(IB), JOFFSET, 1, 1, IRREPX,
     &               372 + 3 - ISIDE)
#ifdef _DEBUG_LVL0
                  call checksum("376",ICORE(IB),NSIZEC)
#endif 
C
                  X = FACT*SDOT(NSIZEC, ICORE(IA), 1, ICORE(IB),
     &               1)
                  XPROP(IPERT, JPERT) = XPROP(IPERT, JPERT) + X
C
                  IF (IPERT .EQ. JPERT) THEN
                     TRACE = TRACE + X
                  ENDIF
C     
               ENDIF
C     
 50          CONTINUE
C
           ENDIF
C
 40      CONTINUE
#ifdef _DEBUG_LVL0
        Write(6,"(a)") " CI-LIke Pol. matrix"
        call output( XPROP,1,NOPERT,1,NOPERT,MAXPRT,MAXPRT,1)
#endif 

C
         IF (ITOP .EQ. 2) THEN
             
            CALL NEATPRINT(LUOUT, XPROP, NOPERT, NOPERT, MAXPRT,
     &                      MAXPRT)
C
            WRITE(6,*)
            IF (ICOUNT .EQ. 1) THEN
              CALL SCOPY(MAXPRT*MAXPRT,XPROP,1,XPROP2,1)
              CALL ZERO(XPROP, MAXPRT*MAXPRT)
            ELSE
              DO I = 1 , MAXPRT
                DO J=1, MAXPRT
                  XPROP(I,J) = XPROP(I,J) + XPROP2(I,J)
                ENDDO
              ENDDO
            ENDIF
          ENDIF

   60   CONTINUE
C
C  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  &  Connected quadratic contributions are required to make               &
C  &  original EOM-CC (JCP 97, 2299, 1993) size extensive. Approximate     &
C  &  but fully size extensive scheme can also be developed by resorting   &
C  &  to keep only connected contributions in original formulation (Sekino &
C  &  Perera and Nooijen unpublish results).                               &
C  &  Connected quadratic contributions were implemented by Ajith in       &
C  &  06/1994 collaboration with Marcel Nooijen. The current               &
C  &  implementation is benefited from the routines written by J. F.       &
C  &  Stanton, J. Gauss and Marcel Nooijen.                                &
C  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     
C Include the quadratic contribution if requested
C
        IF (IEOMPROP .EQ. 2) THEN
C
           itop2 = itop
           if (imag_quad) itop2 = 2
           DO 70 IFREQ = 1, ITOP2
              CALL CALCQUAD1(ICORE(I000), XPROP, MAXCOR, ISYMPRT,
     &             IOFFPRT, NOPERT, IUHF, FACT, MAXPRT,
     &             IFREQ,ITOP,imag_quad)
              if (.false.) then
              write(6,*) ' Quadratic contribution, ifreq ', ifreq
              CALL NEATPRINT(LUOUT, XPROP, NOPERT, NOPERT, MAXPRT,
     &             MAXPRT)
              WRITE(6,*)
              endif
              
 70        CONTINUE
C
        ENDIF
        CALL SYMMET2(XPROP, MAXPRT)
C
C Write the polarizabilities into a file for Raman Intensity 
C calculations. Ajith & John, August 1998. 
C
        OPEN(UNIT=83,FILE='POLAR',STATUS='UNKNOWN',FORM='FORMATTED')
        REWIND(83)
        WRITE(83,3000) (XPROP(I,1), I=1, 3)     
        WRITE(83,3000) (XPROP(I,2), I=1, 3)     
        WRITE(83,3000) (XPROP(I,3), I=1, 3)
 3000   FORMAT(3F20.10)
        CLOSE(UNIT=83,STATUS='KEEP')
C
C For the moment let's assume it can only do polarizabilities.
C After piotr finished his hyperpolarizability we will incorporate
C his work and set up proper printing flags. Ajith October 1996.
        
        CALL HEADER(' Molecular Polarizabilities (in a.u.)', 0, 6)     
C
        CALL NEATPRINT(LUOUT, XPROP, NOPERT, NOPERT, MAXPRT,
     &     MAXPRT)
        WRITE(6,*)
c
c write elements of polarizability to file. At this point frequencies, line with all properties
c
        if (freq_exist) then
         if (ieomprop .eq. 0) then
            iunit = 80
         elseif(ieomprop .eq. 1) then
            iunit = 81
         elseif(ieomprop .eq. 2) then
            iunit = 82
         endif
         write(iunit, 880) freq0
         icount = 0
         do i = 1, nopert
            do j = 1, i
               icount = icount+1
               xvec(icount) = xprop(i,j)
            enddo
         enddo
 880     format(F18.10)
         write(iunit, 881) (xvec(j), j = 1, icount)
 881     format(200F16.8)
         endif
C
   20 CONTINUE
C
      RETURN
      END
