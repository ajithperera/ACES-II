      Subroutine FBlocal(Iuhf, Spntyp, Nalpa, Norbs, Nbasis, Libuf,
     &                   Flag, Dipol, Cdipol, Scfvec, Trans, Scr1, 
     &                   Buf, Ibuf, Itemp, Index)
C
C Description
C
C   Does FB localization (Rev. Mod. Phys. 32,300,(1960)).
C   Basic idea of FB is to maximize the distances between
C   orbital centroids: Maxmize Sum{i=occupied}<i|r|i>**2. The
C   above sum, FB localization sum is a quantitative 
C   measure of degree of localization. In practice this is
C   done by iterative pair wise maximizations until predefined
C   criterion is satisfied.
C
C Variables
C   Lun      : Unit number corresponds to vprop integral file
C   Filnam   : File name for propetry integrals 
C   Dipol    : Area to store dipole integrals from vprop
C   Vplab    : Label for dipole integral in vprop integral file
C   Buf      : Buffer for reading integrals (scratch)
C   Ibuf     : Buffer for reading integrals (scartch)
C   Spntyp   : Distinguish ALPHA and BETA spin cases
C   Norbs    : Number of occupied orbitals
C   Iuhf     : Flag for UHF and RHF runs
C   Nbasis   : Number of basis functions
C   Itemp    : Temporary storge
C   Index    : Keep the ordering of applying 2X2 rotations
C   Cdipol   : Changed dipole integral with rotation
C   CnSSD    : Change in sum of squares of dipole moments with
C              rotations
C   CnSSDM   : Maximum change for particular choice of orbital one
C   IMsave   : First orbital Number corresponds to maximum change
C   JMsave   : Second orbital Number corresponds to maximum change
C   Mxiter  : Maximum number of iterations
C   Fnlsum   : Final localization sum
C   Trans    : Keep the changes of orbital with rotation
C   Nvirt    : Number of virtual orbitals
C   Nalpa    : Number of ALPHA elctrons
C   Flag     : Flag for error return
C Declarations
C
      Implicit Double Precision (A-H, O-Z)
      Dimension Dipol(Nbasis, Nbasis, 3), Cdipol(Nbasis, Nbasis, 3),
     &          Buf(Libuf), Ibuf(Libuf), Scfvec(Nbasis, Nbasis),
     &          Scr1(Nbasis, Nbasis), Trans(Norbs, Norbs),
     &          Index(Norbs), Itemp(Norbs), Sum(3)
      Character*8 Vplab(3)
      Character*5 Filnam, Spntyp
      Logical Flag
      Parameter(Oned10 = 1.0D-10, Oned8 = 1.0D-08, Mxiter = 150,
     &          Maxtry = 3, Fac = 2.541766D00) 
C
C Common Block informations
C
      Common/Machsp/Iintln,Ifltln,Iintfp,Ialone,Ibitwd
      Common/Files/Luout,Moints
C
      Data Vplab /'     X  ', '     Y  ','     Z  '/
C     
      Lun = 7
      Filnam = 'VPOUT'
      Zlich = 0.0D0
      One = 1.0D0
      Two = 2.0D0
      Four = 4.0D0
      Flag = .FALSE.
C     
C     Start Localization trys here
C     
      Itry = 0
C     
 9999 Continue
C     
      Itry = Itry + 1
C     
C     Look for property integral file, open it and read it         
C     
      Call Zero(Dipol, Nbasis*Nbasis*3)
      Do 10 I = 1, 3
         Call Rdvprp(Nbasis, Lun, Filnam, Vplab(I), Dipol(1, 1, I),
     &      Ierr, Buf, Ibuf, Libuf)
         If (Ierr .ne. 0) then
            Flag = .TRUE.
            Return
         Endif
C****Debug
C         Write(Luout,*)
C         Write(Luout,*) 'AO BASIS DIPOLE INTEGRALS'
C         Call Tab(Luout,Dipol(1,1,I),Nbasis,Nbasis,Nbasis,Nbasis)
C*****End
 10   continue
C     
C     AO basis Dipole moment integral are now in DipolX,Y,Z vectros.
C     Transform them to MO Basis by the transformation C(dagger)DC where
C     C is the SCF eigenvectors. They also constitute the orbitals
C     which we want localize.
C     
      If (Spntyp .eq. 'ALPHA') then
         Call Getrec(20, 'JOBARC', 'SCFEVCA0', Nbasis*Nbasis*Iintfp,
     &      Scfvec)
      Else
         Call Getrec(20, 'JOBARC', 'SCFEVCB0', Nbasis*Nbasis*Iintfp,
     &      Scfvec)
      Endif
C     
      Do 20 I = 1, 3
         Call Xgemm('N', 'N', Nbasis, Nbasis, Nbasis, one, 
     &      Dipol(1,1,I), Nbasis, Scfvec, Nbasis, Zlich, 
     &      Scr1, Nbasis)
C     
         Call Xgemm('T', 'N', Nbasis, Nbasis, Nbasis, one, Scfvec,
     &      Nbasis, Scr1, Nbasis, Zlich, Dipol(1,1,I),
     &      Nbasis)
C****Debug
C         Write(Luout,*)
C         Write(Luout,*) 'MO BASIS DIPOLE INTEGRALS'
C         Call Tab(Luout,Dipol(1,1,I),Nbasis,Nbasis,Nbasis,Nbasis)
C*****End
 20   Continue
C     
C     Calculate initial FB localiztion sum
C     
      Bgnlsum = 0.0D0
      Call Zero(Sum, 3)
      Do 30 I = 1, 3
         Call Fbsum( Nbasis, Norbs, Sum(I), Dipol(1,1,I))
         Bgnlsum = Bgnlsum + Sum(I)
 30   Continue

      If (Bgnlsum .Lt. Oned8) Then
          Write(6,"(a,a)") " The oribtal dipole moments are zero.", 
     &                     " This must be an atom or totally symmetric",
     &                     " molecule." 
          Write(6,"(a,a)") " The Foster-Boys localization is not",
     &                     " applicable."
          Call Errex 
      Endif 

C     Print the initial localization sum in Debye**2
C     
      Bgnlsum = Bgnlsum*Fac*Fac
      Write(Luout, 5000) 'Initial localization sum = ', Bgnlsum, 
     &   ' Debye**2'
 5000 Format(2X, A, F12.6, A)
      Write(Luout, *)
      Write(Luout, *)
C     
      Call Zero(Trans, Norbs*Norbs)
      
      Do 40 I = 1, Norbs
         Trans(I, I) = One
 40   Continue
C     
C     Initiate the localization cycles. For each pair of orbitals
C     a 2X2 unitary transformations is performed. The Jacobi type
C     transformation is 
C     Prime(i) = Cos(T)*i + Sin(T)*j
C     Prime(j) = -Sin(T)*i + Cos(T)*j
C     Boys method requires that T should be such as to maximize the sum
C     of squares of orbital dipole moments (SSD, Fancy orbital centroids).
C     Working equations can be obtained from W. L. Lipscomb et al. 
C     JCP 61, 3905 (1974). In practice selecting orbital pairs is random.
C     
C     Seed the random number generator. Also generate an array of
C     random integers of size Norbs. Particular random Genetrator
C     used here is obtained from GAMESS program system. Eventually
C     it will be replaced by an alternative scheme.
C     
      Call Genran(One, Scfvec, Ranum, Nalpa, Nbasis)
      Shift = Atan(One)
      Iter = 0
C
 999  Continue
C     
      Iter = Iter + 1
      Change = Zlich
C
      Do 50 I = 1, Norbs
         Itemp(I) = I
 50   Continue
C     
      NNN = Norbs
      Do 60 I = 1, Norbs
         Call Genran(Zlich, Scfvec, Ranum, Nalpa, Nbasis)
         IXX = Int(Ranum*Real(NNN) + One)
         Index(I) = Itemp(IXX) 
         Itemp(IXX) = Itemp(NNN)
         NNN = NNN - 1
 60   Continue
C     
      Do 70 III = 1, Norbs
         I = Index(III)
         CnSSDM = Zlich
         JMsave = 1
         TMsave = Zlich
         CosTM  = Zlich
         sinTM  = Zlich
C     
         Do 80 J = 1, Norbs
            If ( I .eq. J) Goto 80
C     
C     Compute T which maximize SSD for each I, J pair (I =/ J)
C           
            AIJ = Zlich
            BIJ = Zlich
            Do 90 K = 1, 3
C*****
C               Write(luout,*) 'MO Dipole integrals in loop'
C               Call Tab(Luout, Dipol(1,1,K),Nbasis,Nbasis,Nbasis,Nbasis)
C***** 
               CIJ = Dipol(J, J, K) - Dipol(I, I, K)
               AIJ = AIJ + (Dipol(I, J, K)**2 - CIJ**2/Four)
               BIJ = BIJ + CIJ*Dipol(I, J, K)
 90         Continue
C     
C     Maximum occurs when Tan4(T) = AIJ/BIJ
C     
            If ((Abs(AIJ) .le. Oned10).and.(Abs(BIJ) .le. Oned10))
     &         Goto 80
            T = Atan2(BIJ, AIJ)/Four
C     
C     Impose Sign Convention (For consistancy)
C     
             If ( T .gt. Zlich) then
               T = T - One*Shift
            Else
               T = T + One*Shift
            Endif
C     
C     Update Dipole integrals for this rotation. Also compute the
C     change in SSD for this rotation
C     
            Irot = 0
C     
 99         Continue
C     
            Irot = Irot + 1
            CosT = Cos(T)
            SinT = Sin(T)
C     
            CnSSD = Zlich
            Do 100 K = 1, 3
               Cdipol(I,I,K) = (CosT**2)*Dipol(I,I,K) + (SinT**2)*
     &            Dipol(J,J,K) + Two*SinT*CosT*Dipol(I,J,K)
C     
               Cdipol(J,J,K) = (SinT**2)*Dipol(I,I,K) + (CosT**2)*
     &            Dipol(J,J,K) - Two*SinT*CosT*Dipol(I,J,K)
C     
               CnSSD = CnSSD + (Cdipol(I,I,K))**2 + (Cdipol(J,J,K))**2
     &            - (Dipol(I,I,K))**2 - (Dipol(J,J,K))**2
C
 100        Continue
C     
C     Test to see whether SSD increases with 2X2 rotations. Also 
C     for machine precision
C     
            Tshift = Abs(T) - Shift
            If ((Abs(Tshift) .le. Oned8) .or. (Abs(T) .le. Oned8)
     &         .or. (CnSSD .ge. (-Oned8))) Goto 110
C     
            If ( Irot .le. 1) then
               If ( T .gt. Zlich) then
                  T = T - One*Shift
               Else
                  T = T + One*Shift
               Endif
               Goto 99 
            Endif
C 
C SSD doesn't increase with the rotation. Inform the user and 
C start up a new sequence.
C     
            Write(Luout, 5005)' @Fblocal-I No rotation increases 
     & dipole integrals for the sequance = ', Itry
 5005       Format(2X, A, I1)
            Write(Luout, 5010) I, J, T, CnSSD
 5010       Format(1X, 'Dipole integrals ','I = ',I3, 1X, 'J = ', I3
     &             , 1X, 'Theta = ', F12.6, 1X, 'Change = ', F12.6)
C     
            Return
 110        Continue
C 
            If (CnSSD .gt. CnSSDM) then
               CnSSDM  = CnSSD
               JMsave  = J
               TMsave  = T
               CosTM   = CosT
               SinTM   = SinT
            Endif
C     
 80      Continue
C     
         T    =   TMsave
         J    =   JMsave
         CosT =   CosTM
         SinT =   SinTM
         Change = Change + T**2
C
C Keep track of orbital changes since we need final orbitals
C
         Call Drot(Norbs, Trans(1,I), 1, Trans(1,J), 1, CosT, SinT)
C     
C Now update the dipole integrals corresponds to this
C successful rotation
C     
         Do 120 K = 1, 3
            Cdipol(I,I,K) = (CosT**2)*Dipol(I,I,K) + (SinT**2)
     &         *Dipol(J,J,K) + Two*SinT*CosT*Dipol(I,J,K)
C     
            Cdipol(J,J,K) = (SinT**2)*Dipol(I,I,K) + (CosT**2)
     &         *Dipol(J,J,K) - Two*SinT*CosT*Dipol(I,J,K)
C     
            Cdipol(I,J,K) = (CosT**2 - SinT**2)*Dipol(I,J,K) +
     &         SinT*CosT*(Dipol(J,J,K) - Dipol(I,I,K))
C              
            Cdipol(J,I,K) = Cdipol(I,J,K)
C
C Set other orbital pairs ( Not I or J)
C
            Do 130 MM = 1, Norbs
               If ((I .ne. MM) .and. (J .ne. MM)) Then
C     
                  Cdipol(I,MM,K) = CosT*Dipol(I,MM,K) + SinT
     &               *Dipol(J,MM,K)    
C     
                  Cdipol(MM,J,K) = CosT*Dipol(MM,J,K) - SinT
     &               *Dipol(MM, I, K)
C                           
                  Dipol(I,MM,K) = Cdipol(I,MM,K)
                  Dipol(MM,I,K) = Cdipol(I,MM,K)
                  Dipol(MM,J,K) = Cdipol(MM,J,K)
                  Dipol(J,MM,K) = Cdipol(MM,J,K)
               Endif
 130        Continue
C     
            Dipol(I,I,K) = Cdipol(I,I,K)
            Dipol(J,J,K) = Cdipol(J,J,K)
            Dipol(I,J,K) = Cdipol(I,J,K)
            Dipol(J,I,K) = Cdipol(I,J,K)
C     
 120     Continue
C 
 70   Continue
C
C Calculate the root-mean-square change. Write the current 
C iteration and change to the user
C
      Change = Dsqrt(Two*Change/DFloat(Norbs*(Norbs-1)))
      Write(Luout, 5030) Iter, Change
 5030 Format(15X, 'Boys iteration ', I3, ' Orbital change = ', 
     &   F15.8)
C
C If convergence hasn't been reached start another series of
C two center rotations
C
      If ((Iter .lt. Mxiter) .and. (Change .gt. Oned8)) goto 999
C      
      If (Change .gt. Oned8) then
         If (Itry .lt. Mxtry) then
C     
C Print a warning sign to the user and start a new sequence
C
            Write(Luout, 5040) 
 5040       Format(20X, '!!! Localization Failed. Restarting a new 
     & sequance !!!')
            Goto 9999
         Else
C     
C Localization was not successful, Tell the user about it
C set the flag and return to the caller
C
            Write(Luout, 5050) Maxtry, MXiter, Change
 5050       Format(1X, 'Boys localization failed after ', I3, 'tries 
     & of ', I3, ' iterations each. Last change = ', F15.8)
            Flag = .TRUE.
            Return
         Endif
      Endif
C
C Localization is successful, Tell the user about it
C 
      Write(Luout, *)
      Write(Luout, *)
      Write(Luout, 5060) Iter
 5060 Format(1X, ' Localization Complete after', I3, ' iterations.')
C
C Calculate the new localiztion sum 
C
      Fnlsum = 0.0D0
      Call Zero(Sum, 3)
      Do 150 I = 1, 3
         Call Fbsum( Nbasis, Norbs, Sum(I), Dipol(1, 1, I))
         Fnlsum = Fnlsum + Sum(I)
 150  Continue
C
C  Print the final localiztion sum in Debye**2
      Fnlsum = Fnlsum*Fac*Fac

      Write(Luout, *)
      Write(Luout, *)
      Write(Luout, 5070) 'Final localization sum = ', Fnlsum,
     &   ' Debye**2'
 5070 Format(2X, A, F12.6, A)
C     
#ifdef _DEBUG_LVL0
      Write(Luout,*)
      Write(Luout, 5080) 'The Canonical HF orbitals'
      call output(Scfvec, 1, Nbasis, 1, Nbasis, Nbasis, Nbasis, 1)
#endif
C
C Write out the Transformation matrix and localized MOs
C 
#ifdef _debug_LVL0
      Call Xgemm("T", "N", Norbs, Norbs, Norbs, 1.0D0, Trans,
     &            Norbs, Trans, Norbs, 0.0D0, Scr1, Nbasis)
      Write(6,*) 
      Write(6,"(2x,a)") "Unitary check of FB transformation"
      Call Tab(Luout, scr1, Norbs, Norbs, Nbasis, Nbasis)
#endif

      Write(Luout, *)
      Write(Luout, *)
      Write(Luout, 5080) 'THE LOCALIZTION TRANSFORMATION MATRIX'
 5080 Format(20X, A)
      Call Tab(Luout, Trans, Norbs, Norbs, Norbs, Norbs)
      Write(Luout, *)
      Write(Luout, *)
C
C Transform the canonical to localized sets and copy the virtuals.
C
      Call Zero(Scr1, Nbasis*Nbasis)
      Call Xgemm('N', 'N', Nbasis, Norbs, Norbs, One, Scfvec, Nbasis,
     &           Trans, Norbs, 0.0D0, Scr1, Nbasis)
     
      Nvrt = Nbasis - Norbs 
      Call Dcopy(Nvrt*Nbasis, Scfvec(1,Norbs+1), 1, Scr1(1,Norbs+1), 1)

      Write(Luout, 5090) 'THE LOCALIZED ORBITALS'
 5090 Format(30X, A)
      Call Tab(Luout, Scr1, Nbasis, Nbasis, Nbasis, Nbasis)
C
      Write(Luout, *)
      Write(Luout, *)
C
C Print orbital centroids (Fancy name for orbital dipole moments)
C
      Write(Luout, 5100) 'THE ORBITAL CENTROIDS'
 5100 Format(30X, A)
      Write(Luout, *) 
C      
      Write(Luout, 5110) 'X', 'Y', 'Z'
 5110 Format(18X, A, 22X, A, 22X, A)
      Write(Luout, *)
      Do 160 I = 1, Norbs
         Write(Luout, 5120) Dipol(I,I,1), Dipol(I,I,2),  Dipol(I,I,3)
 160  Continue
C
 5120 Format(8X, F15.8, 8X, F15.8, 8X, F15.8)
      Write(Luout, *)

C Evrey thing is done. So Return
C
      Return
      End

      
      
