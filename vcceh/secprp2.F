C     
      SUBROUTINE SECPRP2(ICORE, MAXCOR, IUHF)
C
C Similar to SECPRP, except that it is based on VDINT integrals.
C
      IMPLICIT INTEGER (A-Z)
C
      LOGICAL IDENTITY
      LOGICAL CCSD, MBPT, PARTEOM
      LOGICAL DOLEFT, DORIGHT, NODAVID, INCORE, REALFREQ,
     &        BLOCK
      LOGICAL BUILD_DENS
      LOGICAL SAVE_VECS
C     
      DOUBLE PRECISION SDOT,X1,X00,MEXPECT,MBAR00,H0QM,HGGM,X
      DOUBLE PRECISION ONE,TWO,ZILCH,FACT,FREQ,TRACE,ZZ,XPROP
      DOUBLE PRECISION YPROP
      LOGICAL PRINT,EXTRAPOL
C
      CHARACTER*8 LABEL
      CHARACTER*6 SIDES(2)
      CHARACTER*8 PROJECT(3)
      CHARACTER*80 STRING
C
      PARAMETER(MAXPRT=600)
C
      DIMENSION ICORE(MAXCOR),NPERT(8),LENGTH(8)
      DIMENSION IWORK(1)
      DIMENSION MEXPECT(MAXPRT),MBAR00(MAXPRT)
      DIMENSION XPROP(MAXPRT,MAXPRT), YPROP(MAXPRT*MAXPRT)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT (8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FLAGS2/IFLAGS2(500)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/PERTTYP/LABEL(100)
      COMMON/FILES/LUOUT,MOINTS
      COMMON/EOMINFO/CCSD, MBPT, PARTEOM, NODAVID
C
      DATA ONE,TWO,ZILCH /1.0D0,2.0D0,0.0D0/
      DATA SIDES /'right.','left. '/
      DATA PROJECT /'CORE    ','VALENCE ','INNERVAL'/        
C
      INDXF(I,J,N)=I+(J-1)*N    
C
      PRINT = IFLAGS(1) .GE. 0
      IOFF = 1     
      IROW = 1
      ICOL = 1
      IONE = 1
      FREQ = 0.0d0
      REALFREQ   = .TRUE.
      FACT = TWO
      ISIDE  = 1
C
C Allocate 10,000 integers. This is useful in new amplitude extrapolation
C algorithm.

      CALL GETREC(0,"JOBARC","ETA_4EXP",ILENGTH,IJUNK)
      IF (ILENGTH .GT. 0) EXTRAPOL = .TRUE.
      IF (EXTRAPOL) THEN
          IMAXCOR  = 10000
          IWORK(1) = 0
          CALL ACES_MALLOC(IMAXCOR,IWORK,J000)
          IF (IWORK(1).eq.0) then
              PRINT*, '@-secrprp2.F: Unable to allocate at least ',
     &                 IMAXCOR, ' integers of memory'
             CALL ERREX
          ENDIF 
      ENDIF

      CALL GETREC(20, 'JOBARC', 'NREALATM', 1, NATOMS)
      IF (NATOMS .GT. 100) THEN
          WRITE(6,"(A)") " This program has a 100 atoms limit."
          CALL ERREX
      ENDIF 

      IF (IFLAGS(18) .EQ. 8)  ICALL  = 0
      IF (IFLAGS(18) .EQ. 9)  ICALL  = 2
      IF (IFLAGS(18) .EQ. 10) ICALL  = 3
C
C Set the flags for various methods (CI like, Linear and
C quadratic)
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
     &    IEOMPROP .EQ. 2) THEN
C
         IFIRST = IEOMPROP
         ILAST  = IEOMPROP
C
      ELSE IF (IEOMPROP .EQ. 3) THEN
C
         IFIRST = 0
         ILAST  = 2 
C
      ENDIF
C
C Flag for antisymmetric perturbations (Paramagnetic spin-orbit 
C interaction)
C
      IF (IFLAGS(18) .EQ. 8) THEN
         IANTI = 1
      ELSE
         IANTI = 0
      ENDIF
C
      NBAS = NOCCO(1) + NVRTO(1)
      CALL GETREC(20,'JOBARC','NBASTOT ',IONE,NBAST)
C
      KMAX = MAXEXP - 1
C
      I000  = 1
      I010  = I000 + MAX(NBAST, KMAX*KMAX)*IINTFP
      I020  = I010 + KMAX*IINTFP
C
C Memory checks were added by Ajith 07/95
C 
      IF (I020 .GT. MAXCOR) CALL INSMEM('SECPRP2', I020, MAXCOR) 
      MXCOR = (MAXCOR-I020+1)/IINTFP
C   
      CALL GETREC(20, 'JOBARC', 'NTOTPERT', 1, NTPERT) 
      CALL GETREC(20, 'JOBARC', 'NPERTIRR', 8, NPERT) 
      CALL GETREC(20, 'JOBARC', 'NUCIND  ', 1, NUCIND) 
C
      WRITE(6,*)
      WRITE(6, 5000) SIDES(ISIDE)
 5000 FORMAT(T3, '@SECPRP2-I, Linear equations solved from ', A6)
      WRITE(6,*)
C      
      INEWFIL = 0
C
      DO 2 IRREP = 1, NIRREP
C
         WRITE(6, 9000) NPERT(IRREP),IRREP
 9000    FORMAT(T3,'@SECPRP2-I, There are ',I2,' perturbations ',
     &          ' of symmetry species ',I2,'.')
C
         LENGTH(IRREP) = 0
C
         DO 3 ISPIN = 1, 1 + IUHF
            LENGTH(IRREP) = LENGTH(IRREP)+IRPDPD(IRREP,8+ISPIN)
 3       CONTINUE
C
         DO 4 ISPIN = 3, 3-2*IUHF, -1
            LSTTYP = 43 + ISPIN
            LENGTH(IRREP) = LENGTH(IRREP) +
     &                      IDSYMSZ(IRREP, ISYTYP(1,LSTTYP),
     &                      ISYTYP(2,LSTTYP))
 4       CONTINUE
C
C In John stanton's original implementation these list numbers 
C were 473 and 474. They were changed to 300 lists to accommodate
C 400 list's to the quadratic contributions. 06/1994 by Ajith
C Note that in the first UPDMOI call to the list 373, it opens
C the 'DERINT' file.
C The list number 375 is added to carry out the core, inner-core
C and valence projections Ajith 01/1997.
C
         CALL UPDMOI(NPERT(IRREP), LENGTH(IRREP), IRREP, 373, 
     &               INEWFIL, 0) 
         INEWFIL = 0
         CALL UPDMOI(NPERT(IRREP), LENGTH(IRREP), IRREP, 374,
     &               INEWFIL, 0)
         CALL UPDMOI(3, LENGTH(IRREP), IRREP, 375, INEWFIL, 0)

         call aces_list_memset(irrep,373,0)
         call aces_list_memset(irrep,374,0)
         call aces_list_memset(irrep,375,0)

         IF (EXTRAPOL) THEN
            INULL = 0 
            LENOO = IRPDPD(1,21)
            LENVV = IRPDPD(1,19)
            CALL UPDMOI(NPERT(IRREP),LENOO,IRREP,380,INULL,INULL)
            CALL UPDMOI(NPERT(IRREP),LENVV,IRREP,381,INULL,INULL)
            IF (IUHF .NE. 0) THEN
               LENOO = IRPDPD(1,22)
               LENVV = IRPDPD(1,20)
               CALL UPDMOI(NPERT(IRREP),LENOO,IRREP,382,INULL,INULL)
               CALL UPDMOI(NPERT(IRREP),LENVV,IRREP,383,INULL,INULL)
            ENDIF 
         ENDIF 
C
 2    CONTINUE
C
      WRITE(6,*)
C
      DO 5 IRREPX = 1, NIRREP
C     
         NUMPERT = NPERT(IRREPX)
         IBEGIN  = 1
C
         IF (NUMPERT .GT. 0 ) THEN
C
            NDIM = IRPDPD(IRREPX, 9)
            N1AA = NDIM
            IF (IUHF .NE. 0) NDIM = NDIM + IRPDPD(IRREPX,10)
C
            CALL NEWLST(IRREPX, ICORE, MAXCOR, IUHF)
C
            WRITE(6,2400) IRREPX
 2400       FORMAT (T3,' Solving for perturbations in irrep ', I3)
            WRITE(6,2500) NDIM, NSIZEC
 2500       FORMAT(T3, 'Dimension singles ',I4,' Total Dimension ',I8)
C
C When requsted (IFLAGS(129)=1) to have an estimate of core, inner valence 
C and valence contributions to second-order properties, we need to divide 
C the determinatal space accordingly. The EXCITEPROJ and dependencies 
C do the job.
C 
            IF (IFLAGS2(129) .EQ. 1 .AND. IEOMPROP .EQ. 0) THEN
C
               DO 35 IPRJECT = 1, 3
C
                  CALL EXCITEPROJ(ICORE, MAXCOR/IINTFP, IRREPX, 
     &                            PROJECT(IPRJECT), NSIZEC, IUHF)
C
                  CALL PUTEXCITEPRJ(ICORE, MAXCOR/IINTFP, IRREPX, 
     &                                375, NSIZEC, IBEGIN, IUHF)
                  IBEGIN = IBEGIN + 1 
C
 35            CONTINUE
C
            ENDIF
C
C We need to avoid double counting which is unavoidable in the way
C I have programmed for convenience. Let's take care of the double
C counting now.
C
            CALL FIXDBLCONT(ICORE, MAXCOR/IINTFP, IRREPX, NSIZEC)
C            CALL PROJECTCHECK(ICORE, MAXCOR/IINTFP, IRREPX, NSIZEC)
C
            DO IPERT = 1, NUMPERT
C     
               IOFFSET = IPERT
               DORIGHT = ISIDE. EQ. 1
               DOLEFT  = ISIDE .EQ. 2
C
               CALL OPRBAR(ICORE, MAXCOR, IUHF, LABEL(IPERT),
     &                     IRREPX, IOFFSET, IPERT, MEXPECT(IPERT),
     &                     MBAR00(IPERT), .FALSE.,
     &                     DOLEFT, DORIGHT, IEOMPROP)
C     
            ENDDO
C
            IF (PARTEOM .AND. NODAVID) THEN
C
C CALCULATE MEMORY NEEDED FOR INCORE (BLOCK) DAVIDSON PROCEDURE 
C
               BLOCK = .TRUE.
               KMAX  = MIN(N1AA, MAXEXP-1)
C
               IF (BLOCK) THEN
                  IZ = 1
                  IHZ = IZ + KMAX * NDIM * IINTFP
                  IY1 = IHZ + KMAX * NDIM * IINTFP
                  IZOUT = IY1 + NDIM * NUMPERT * IINTFP
                  IH0 = IZOUT + NDIM * NUMPERT * IINTFP
                  ICONV = IH0 + NDIM*IINTFP
                  ISCR = ICONV + NUMPERT*IINTFP
               ELSE
                  IC = I020
                  IHXC = IC + NDIM * KMAX * IINTFP
                  IY1 =  IHXC + NDIM * KMAX * IINTFP
                  IX1 = IY1 + NDIM * IINTFP
                  IH0 = IX1 + NDIM * IINTFP
                  ISCR = IH0 + NDIM * IINTFP
               ENDIF
C
               INCORE = .TRUE.            
               MXCOR2 = (MAXCOR - ISCR + 1) / IINTFP
               IMAX = ISCR + (2 * NSIZEC + NSIZEC /4)* IINTFP
C
               IF (IMAX .GE. MAXCOR) THEN
                  CALL INSMEM('SECPRP', IMAX, MAXCOR)
                  INCORE = .FALSE.
               ENDIF
C
            ELSE
               INCORE = .FALSE.
            ENDIF
C
            IF (.NOT. INCORE) BLOCK = .FALSE.
            IF (BLOCK) THEN
               WRITE(6,2600) ISCR
 2600          FORMAT(T3, ' Memory used in Blcklineq ', I10)
            ELSE
               IF (INCORE) WRITE(6,2700) ISCR
 2700          FORMAT(T3, ' Memory used in plineqy ', I10)
            ENDIF
C
C Perform partitioned EOM calculation and built P-EOM densities
C and obtain the P-EOM natural orbitals.

            IF (EXTRAPOL) THEN
                BUILD_DENS = .TRUE.
                SAVE_VECS  = .TRUE.
                IF (3*NSIZEC .GT. MAXCOR) THEN
                    Write(6,"(a,ES17.9E2,a)") 
     &                             " The current implementation",
     &                             " of amplitude extrapolation",
     &                             " assumes having at least 3",
     &                             " copies of ",Nsizec," vectors."
                    Call INSMEM("secprp2",NSIZEC,MAXCOR)
                Endif 
                CALL GET_PEOM_NOS(ICORE,MAXCOR,IWORK(J000),
     &                            IMAXCOR,IUHF,NPERT,IPERT,
     &                            NUMPERT,IRREPX,ISIDE,
     &                            BUILD_DENS,SAVE_VECS)
C                 -
C Reconstruct the M vectors (this step is cheap to redo instead of
C storing a copy) 
C
               DO IPERT = 1, NUMPERT
                  IOFFSET = IPERT
                  DORIGHT = ISIDE. EQ. 1
                  DOLEFT  = ISIDE .EQ. 2
                  CALL OPRBAR(ICORE,MAXCOR,IUHF,LABEL(IPERT),
     &                        IRREPX,IOFFSET,IPERT,MEXPECT(IPERT),
     &                        MBAR00(IPERT),.FALSE.,DOLEFT,
     &                        DORIGHT,IEOMPROP)
               ENDDO
               DO I=1,2
                  CALL UPDMOI(MAXEXP,NSIZEC,I,470,0,0)
                  CALL UPDMOI(MAXEXP,NSIZEC,I,471,0,0)
               ENDDO

            ENDIF 
C
C SOLVE LINEAR SYSTEM OF EQUATIONS.
C
            IF (BLOCK) THEN

C The extrapolation needs to tell the blcklineq to save/no save the solution
C vector to list 372+Iside.

               SAVE_VECS = .TRUE.
               CALL BLCKLINEQ(IUHF, IRREPX, ISIDE,
     &                        NUMPERT, KMAX,
     &                        N1AA, NDIM, ICORE(IZOUT), ICORE(IY1),
     &                        ICORE(IZ), ICORE(IHZ), ICORE(IH0),
     &                        ICORE(ICONV), ICORE(ISCR), MXCOR2, 
     &                        FREQ, REALFREQ, SAVE_VECS)
            ELSE
C
               DO 10 IPERT = 1, NUMPERT
C
                  IOFFSET = IPERT
                  WRITE(6,500) IPERT
                  WRITE(6,501)
C
                  IF (INCORE) THEN
                     CALL PLINEQY(ICORE(ISCR),MXCOR2,IUHF,
     &                            ICORE(I000),ICORE(I010),KMAX,
     &                            IRREPX,IOFFSET,IPERT,FREQ,ISIDE,
     &                            REALFREQ,ICORE(IC), ICORE(IHXC),
     &                            ICORE(IY1),ICORE(IX1),ICORE(IH0),
     &                            N1AA, NDIM)
                  ELSE
C
                     I000  = 1
                     I010  = I000 + MAX(NBAST, KMAX*KMAX)*IINTFP
                     I020  = I010 + KMAX*IINTFP
                     MXCOR = (MAXCOR-I020+1)/IINTFP
C
                     CALL LINEQY(ICORE(I020), MXCOR, IUHF, 
     &                           ICORE(I000),ICORE(I010), KMAX,
     &                           IRREPX, IOFFSET, IPERT, FREQ, ISIDE,
     &                           REALFREQ)
                  ENDIF
C
 500              FORMAT(T3,'@SECPRP2-I, Current Perturbation # is '
     &                   ,I2,'.')
 501              FORMAT(T3,'@SECPRP2-I, Solving for derivative T',
     &                   ' amplitudes.')
C
C Transform perturbed CCSD T2 amplitudes to the P-EOM natural orbitals
C  
                  IF (EXTRAPOL) THEN
                     CALL FORM_NOTRANS_PT2(ICORE(I000),MAXCOR,
     &                                     IWORK(J000),IMAXCOR,
     &                                     IUHF,IRREPX,NSIZEC,
     &                                     IOFFSET,NPERT,ISIDE)
                    
                  ENDIF  

 10            CONTINUE
C
            ENDIF
C
         ENDIF
C
 5    CONTINUE
C
C Calculate properties 
C      
      DO 20 IEOMPROP = IFIRST, ILAST
C 
C Begins the the projection loop. Do it only when you do this
C type of calculations.
C
         IF (IFLAGS2(129) .EQ. 1 .AND. IEOMPROP .EQ. 0) THEN
            IPRJFRST = 1
            IPRJLAST = 3
         ELSE
            IPRJFRST = 1
            IPRJLAST = 1
         ENDIF
C         
         DO 45 IPRJECT = IPRJFRST, IPRJLAST
C
            IF (ISIDE .EQ. 2 .AND. IEOMPROP .EQ. 2) GOTO 20
C
            CALL ZERO(XPROP, MAXPRT*MAXPRT)
            CALL ZERO(YPROP, MAXPRT*MAXPRT)
            IOFFPT1 = 0
            IOFFPT2 = 0
C
            DO 30 IRREPX = 1, NIRREP
C     
               NUMPERT = NPERT(IRREPX)
C
               DO 40 IPERT = 1, NUMPERT
C
                  IOFFSET = IPERT
C               
                  IF (IRREPX .NE. 0) THEN
C
                     CALL NEWLST(IRREPX, ICORE, MAXCOR, IUHF)
C     
                     DOLEFT  = ISIDE .EQ. 1
                     DORIGHT = ISIDE .EQ. 2
                     CALL OPRBAR(ICORE, MAXCOR, IUHF, LABEL(IPERT),
     &                           IRREPX, IOFFSET, IPERT,
     &                           MEXPECT(IPERT),MBAR00(IPERT),
     &                           .FALSE. ,DOLEFT, DORIGHT, IEOMPROP)
                  ENDIF
C
 40            CONTINUE
C
               IF (IFLAGS2(129) .EQ. 1 .AND. IEOMPROP .EQ. 0) THEN
                  IP = I000
                  IA = IP + NSIZEC*IINTFP
                  CALL GETLST(ICORE(IP), IPRJECT, 1, 1, IRREPX, 375)
               ENDIF
C               
               DO 50 IPERT = 1, NUMPERT
C     
                  DO 60 JPERT = 1, NUMPERT
C
                     IOFFSET = IPERT
                     JOFFSET = JPERT
C
                     CALL NEWLST(IRREPX, ICORE, MAXCOR, IUHF)
C                     
                     IF (IFLAGS2(129) .EQ. 1 .AND. IEOMPROP .EQ. 0)
     &                   THEN
                        IB = IA + NSIZEC*IINTFP
                     ELSE
                        IA = I000
                        IB = IA + NSIZEC*IINTFP
                     ENDIF
C
                     MXCOR = MAXCOR - IB + 1
C
                     CALL GETLST(ICORE(IA), IOFFSET, 1, 1, IRREPX,
     &                           374)
C
                     IF (IUHF .EQ. 0) THEN
                        CALL SPNTSING(IRREPX, ICORE(IA), ICORE(IB), 
     &                                MXCOR)
                     ENDIF
C
                     IF (IFLAGS2(129) .EQ. 1 .AND. IEOMPROP .EQ. 0)
     &                   THEN 
                        IC     = IB + NSIZEC*IINTFP
                        INEED  = NSIZEC*IINTFP + 1 
                        IAVAIL = MAXCOR - IB 
                        IF (IC .GT. MAXCOR) CALL INSMEM("@-SECPRP2", 
     &                                                 INEED, IAVAIL)
                        CALL VECPRD(ICORE(IP), ICORE(IA), ICORE(IC), 
     &                              NSIZEC)
                        CALL SCOPY(NSIZEC, ICORE(IC), 1, ICORE(IA), 1)
                     ENDIF
C
                     CALL GETLST(ICORE(IB), JOFFSET, 1, 1, IRREPX,
     &                           373)
C
#ifdef _DEBUG_LVLM
      call checksum("secprp2", icore(ib), nsizec)
#endif 
                     X = FACT*SDOT(NSIZEC, ICORE(IA), 1, ICORE(IB),
     &                             1)
C     
                     XPROP(IPERT + IOFFPT1, JPERT + IOFFPT1) = X
C     
 60               CONTINUE
 50            CONTINUE
C
               IOFFPT1 = IOFFPT1 + NUMPERT
C
 30         CONTINUE
             
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
C Print the description of the program, What it does, authors etc.
C to the user.
C
            IF (IEOMPROP .EQ. 2) THEN
               CALL BANNER(LUOUT)
            ENDIF
C
C Create list Hbar(MB,EJ), Hbar(Mb,Ej) and Hbar(Ab,Cd) for RHF
C calculations. These lists are used in the quadratic term, and
C for general comments see below.
C 
            IF (IUHF .EQ. 0 .AND. IEOMPROP .EQ. 2) THEN
               CALL GENRHFLST(ICORE, MAXCOR, IUHF)
            ENDIF
C
C Include the quadratic contribution if requested
C     
            DO 70 IRREPX = 1, NIRREP
C
               NUMPERT = NPERT(IRREPX)
C
               IF (IEOMPROP .EQ. 2 .AND. IRREPX .NE. 0 .AND. 
     &             NUMPERT .NE. 0) THEN
C     
                  CALL CALCQUAD2(ICORE(I000), XPROP, MAXCOR, IRREPX, 
     &                           NUMPERT, IUHF, FACT, IOFFPT2, MAXPRT)
C 
               ENDIF
C     
               IOFFPT2 = IOFFPT2 + NUMPERT
C
 70         CONTINUE
C
C Print the results to the user with appropriate information
C so the that user knows what he is geeting out of the machine.
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
C
            Write (LUOUT, '(1X,A)') String(2:Len(string)-1)
C
         if (ieomprop .eq. 0) then
           string='$CI-like approximation to evaluation of properties.'
           else if (ieomprop .eq. 1) then
           string='$Linear approximation to evaluation of properties.'
           else if (ieomprop .eq. 2) then
           string='$Fully quadratic evaluation of properties.'
         endif
C
            Write (LUOUT, '(1X,A)') String(2:Len(string)-1)
C
         if (iflags2(129) .eq. 1 .and. ieomprop .eq. 0) then
          if (iprject .eq. 1) then 
           string='$Core contribution to the coupling tensor'
          else if (iprject .eq. 2) then 
           string='$Valence contribution to the coupling tensor'
          else if (iprject .eq. 3) then
           string='$Inner-Valence contribution to the coupling tensor'
          endif                
          Write (LUOUT, '(1X,A)') String(2:Len(string)-1)
         endif
C     
            IF (IFLAGS(18) .NE. 9 .AND. IFLAGS(18) .NE. 10 .AND. 
     &          IFLAGS(18) .NE. 8) THEN
C
               CALL NEATPRINT(LUOUT, XPROP, NOPERT, NOPERT, MAXPRT,
     &                        MAXPRT)
               WRITE(6,*)
C     
            ENDIF
C           
            IF (IFLAGS(18) .EQ. 9 .OR. IFLAGS(18) .EQ. 10 .OR. 
     &          IFLAGS(18) .EQ. 8) THEN
C
               JOFF = 1
C     
C Symmetrize the out-put matrix.
C     
               DO 80 IOFF = 1, NTPERT
C             
                  CALL SCOPY(NTPERT, XPROP(1, IOFF), 1, YPROP(JOFF),
     &                       1)
C
                  JOFF = JOFF + NTPERT
C
 80            CONTINUE
C
               IF (IANTI .NE. 0) CALL VMINUS(YPROP, NTPERT*NTPERT)
               CALL SYMMET2(YPROP, NTPERT)
C
               I000 = 1
               I010 = I000 + IINTFP*NTPERT*NTPERT
               I020 = I010 + IINTFP*NTPERT*NTPERT
               IF (I020 .GT. MAXCOR) CALL INSMEM("secprp2",I020,MAXCOR)
C
               CALL JJSPNCUP(YPROP, ICORE(I000), ICORE(I010), NTPERT, 
     &                       NIRREP, IEOMPROP,  ICALL)
            ENDIF  
C
 45      CONTINUE
C     
 20   CONTINUE
C
 900  CONTINUE
C
      CLOSE(UNIT=10, STATUS='KEEP')
C     
      RETURN
      END
