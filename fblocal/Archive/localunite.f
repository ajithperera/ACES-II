      Program Local
C
C Description
C
C   Perform Foster-Boys and Edmiston-Ruedenberg
C   Localizations. Works with both UHF and RHF 
C   reference functions. Frozen core and virtual
C   localiztions are available in near future.
C   Current version was coded by Ajith 05/17/93
C
C Local Variables
C   Itemp  : Temporary storge area
C   Nocupy : Keep No. of occupied orbitals for RHF and UHF runs
C   Spntyp : Corresponds to ALPHA  and BETA spin cases
C   Nbasis : Number of basis functions
C   Nalpa  : No. of ALPHA occupied orbitals
C   
C Declarations
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Nocupy(2), Itemp(12)
      Character*5 Spntyp(2)
      Character*8 Dump(2)
      Logical Abort
      Parameter(Lirvpr=600, Oned10 = 1.0D-10, Oned8 = 1.0D-8, 
     &          MXNiter = 150, Mxtry = 3)
C
C Common Block informations
C
      Common//Icore(1)
      Common/Istart/I0
      Common/Machsp/Iintln,Ifltln,Iintfp,Ialone,Ibitwd
      Common/Files/Luout,Moints
      Common/Info/Nocco(2),Nvrto(2)
      Common/Flags/ Iflags(100)
      Common/Syminf/Nstart,Nirrep,Irreps(255,2),Dirprd(8,8)
C
      Data Spntyp /'ALPHA', 'BETA'/
      Data Dump /'LOCALMAO', 'LOCALMBO'/

C Initialize the Memory Allocation
C
      Call Crapsi(Icore, Iuhf, 0)
C
      Ione = 1
      Ithree = 3
      Nbasis = Nocco(1) + Nvrto(1)
C  
      Call Izero(Nocupy, 2)
      Call Getrec(20, 'JOBARC', 'SYMPOPOA', Nirrep, Itemp)
C
      Do 10 I = 1, Nirrep
         Nocupy(1) = Nocupy(1) + Itemp(I)
 10   Continue
      
      If (Iuhf .ne. 0) then 
         Call Getrec(20, 'JOBARC', 'SYMPOPOB', Nirrep, Itemp)
         Do 20 I = 1, Nirrep
            Nocupy(2) = Nocupy(2) + Itemp(I)
 20      Continue
      Endif
C
      Nalpa = Nocupy(1)
C
C Write out some important stuff to the user
C
      Write(Luout,*)
      Write(Luout, 5000)'FOSTER-BOYS ORBITAL LOCALIZATION'
      Write(Luout, 5000)'--------------------------------'
 5000 Format(25X, A)
      Write(Luout,*) 
      Write(Luout,*)'  Maximum number of localization retrys = ', Mxtry
      Write(Luout,*)'  Maximum number of iteration per try = ', MXNiter
      Write(Luout,5010)'Threshold for neglect a rotation = ', Oned8
      Write(Luout,5010)'Threshold for convergence = ', Oned8
      Write(Luout,*)
 5010 Format(2X, A, E7.1)
C      
C Allocate the memory
C      
      I000 = I0
      I010 = I000 + Ithree*Nbasis*Nbasis*Iintfp
      I020 = I010 + Ithree*Nbasis*Nbasis*Iintfp
      I030 = I020 + Nbasis*Nbasis*Iintfp
      I040 = I030 + Nbasis*Nbasis*Iintfp
      I050 = I040 + Nbasis*Nbasis*Iintfp
      I060 = I050 + Lirvpr*Iintfp
      I070 = I060 + Lirvpr
      I080 = I070 + Norbs
      I090 = I080 + Norbs
C
      If (Iuhf .ne. 0) then
         
         Do 30 I = 1, 2
            Write(Luout, 5020)'-UHF WAVE FUNCTION : SPIN TYPE',
     &         Spntyp(I),'-'
 5020       Format(22X, A, 1X, A, A)
            Write(Luout, *)
             Call FBlocal(Iuhf, Spntyp(I), Nalpa, Nocupy(I), Nbasis, 
     &                   Lirvpr, Abort, Icore(I000), Icore(I010),
     &                   Icore(I020), Icore(I030), Icore(I040), 
     &                   Icore(I050), Icore(I060), Icore(I070),
     &                   Icore(I080))
 30      continue

C Write localized MO's into JOBARC file
C
      Call Putrec(20, 'JOBARC', Dump(I), Nbasis*Nbasis*Iintfp,
     &            Icore(I040))
C      
      Else
C
         Write(Luout, 5030)'-RHF WAVE FUNCTION-' 
         Write(Luout, *) 
 5030    Format(30X, A) 
         Call FBlocal(Iuhf, Spntyp(1), Nalpa, Nocupy(1), Nbasis, 
     &                   Lirvpr, Abort, Icore(I000), Icore(I010),
     &                   Icore(I020), Icore(I030), Icore(I040), 
     &                   Icore(I050), Icore(I060), Icore(I070),
     &                   Icore(I080))
C
C Write localized MO's into JOBARC file for RHF runs
C
      Call Putrec(20, 'JOBARC', Dump(I), Nbasis*Nbasis*Iintfp,
     &            Icore(I040))

      Endif
C
C Write out the successful completion message and also 
C the time taken
C 
      If ( Abort) Then
         Call Errex
         Call Crapso
      Else
         Write(Luout,*)
         Call Timer(Tcpu, Tjunk, 1)
         Write(Luout,5040) Tcpu
 5040    FORMAT(T3,'@LOCAL-I, Orbital localization completed in'
     &   ,F8.3,' cpu seconds.')
C
         Write(Luout, *)
         Call Crapso
      Endif
C
      Stop
      End
C
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
C     
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
 5005       Format(2X, A, I)
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
C
C  
      Write(Luout, *)
      Write(Luout, *)
      Write(Luout, 5070) 'Final localization sum = ', Fnlsum,
     &   ' Debye**2'
 5070 Format(2X, A, F12.6, A)
C
      Call Zero(Scr1, Nbasis*Nbasis)
      Call Xgemm('N', 'N', Nbasis, Nbasis, Norbs, One, Scfvec, Nbasis, 
     &           Trans, Norbs, Zlich, Scr1, Nbasis)
C
C Copy the virtuals in Scr1 
C
      Nvirt = Nbasis - Norbs
      Call Dcopy(Nvirt*Nbasis, Scfvec(1, Norbs+1), 1, Scr1(1, Norbs+1),
     &   1)
C
C Write out the Transformation matrix and localized MOs
C 
      Write(Luout, *)
      Write(Luout, *)
      Write(Luout, 5080) 'THE LOCALIZTION TRANSFORMATION MATRIX'
 5080 Format(20X, A)
      Call Tab(Luout, Trans, Norbs, Norbs, Norbs, Norbs)
      Write(Luout, *)
      Write(Luout, *)
C
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
C      
      Subroutine Fbsum(Nbasis, Norbs, Sum, Dipole)
C
C Description
C    Calculate the FB localization sums : Sum{occupy}<i|ri|>**2
C 
C Variables
C    Nbasis    : Number of Basis functions
C    Norbs     : Number of occupied orbitals  
C    Sum       : Final sum returns in sum
C    Dipole    : MO dipole integrals
C
C Declarations
C 
      Implicit Double Precision (A-H, O-Z)
      Dimension Dipole(Nbasis, Nbasis)
C 
C Common block informations
C
      Common/Machsp/Iintln,Ifltln,Iintfp,Ialone,Ibitwd
      Common/Files/Luout,Moints
C 
C     
      Sum = 0.0D0
      Do 10 J = 1, Norbs
         Sum = Sum + Dipole(J, J)**2
 10   Continue
C
      Return
      End
C
      Subroutine Filabl(Lun, Label, Ierr)
C
C Description
C
C   Look for a record label in property integral 
C   program and place the file position pointer 
C   at that label.
C
C Variables
C
C   Lun   = Unit number corresponds to vprop integral file (Input)
C   Label = Record label (input)
C   Ierr  = Error Handler 
C            = 0 Success
C            = 1 unable to find the record label
C Declarations
C
      Implicit Double Precision (A-H, O-Z)
      Character*8 Star, Blank, Check(4), Label
      Parameter (Star  ='********')
      Parameter (Blank ='        ')
C
C Common Block informations
C
      Common/Machsp/Iintln,Ifltln,Iintfp,Ialone,Ibitwd
      Common/Files/Luout,Moints
C
      Ierr = 0
 100  Read(UNIT=Lun, END=1000, ERR=2000, IOSTAT=IOS) Check(1)
      If (Check(1) .ne. Star) Goto 100
C
C Now check the rest to make sure ...
      Backspace (Lun)
      Read(UNIT=Lun, END=1000, ERR=2000, IOSTAT=IOS) Check
      If (Check(4) .eq. Label) Return
      Goto 100
C
C Handle Errors
C
 1000 Ierr = 1
      Return
C     
 2000 Continue

C Get here via an error reading the file
C
      Write(luout,*) ' @Filabl-I, Read error on file attched to
     &                unit ', Lun, '.'

      Return
      End
      Subroutine genran(Test, Vecin, Ranum, Nalpa, Nbasis)
C
C Description 
C   Generate a random number. The present
C   version is obtained from GAMESS program
C   system
C 
C Local variables
C   Vecin  : SCF eigenvectors
C   Nalpa  : Number of ALPHA electrons
C   Nbasis : Number of basis functions
C   Ranum  : Random number (output)
C   U, XY  : Arbitrary labels
C
C Declarations
C  
      Implicit Double Precision (A-H, O-Z)
      Dimension Vecin(Nbasis, Nbasis)
      Parameter(Zero = 0.0D0, One = 1.0D0)
C
      Save U
C
      PI = Acos(-One)
C
      If (Test .ne. Zero) Then
         N = Abs(Nalpa - Nbasis) + 1
         M = N + 5
         XY = Vecin(N, M) * Atan(one)
         U  = (PI + XY)**5
         XY = Real(Int(U))
         U =  U - XY
         Ranum = U
      else
         U = (PI + U)**5
         XY = Real(Int(U))
         U = U - XY
         Ranum = U
      Endif
C
      Return
      End
C         
      subroutine Opnfil(Lun, Filnam, Format, Iopned, Ierr)

C Description
C   Open files attach to a particular unit number
C 
C Variables
C   Lun     : Unit number corresponding to Vprop integral file (input)
C   Filnam  : File name for property integrals
C   Format  : Structure of the file (input)
C   Iopned  : Check to see the file is being used(input) 
C              = 0 then close it and continue
C              = 1 Print a error message and abort
C   Ierr    : Return code
C              = 0 Succesful completion
C              = 1 unable to open/read file
C
C   Ioerr   : Holds return from IOSTAT
C
C Declarations
      Implicit Double Precision(A-H, O-Z)
      Character*(*) Filnam, Format
      Logical Isopen
C
C Common block informations
C
      Common/Machsp/Iintln,Ifltln,Iintfp,Ialone,Ibitwd
      Common/Files/Luout,Moints
     
C See whether particular file exists
C
      Inquire(UNIT=Lun, OPENED=Isopen)
      If (Isopen) then 
         Ierr = 2 + Iopned
         If (Iopned .eq. 1) Return
         Close(UNIT=Lun, IOSTAT=Ioerr)
         If (Ioerr .ne. 0) then
            Ierr = 5
            Write(Luout,*) '@Opnfil-I I/O error on Unit ', Lun,'.'
            Return
         Endif
      Endif
         
C It is not opend, good deal
C
      Open(UNIT=Lun, FILE=Filnam, FORM=Format, Access='SEQUENTIAL',
     &     STATUS='UNKNOWN',IOSTAT=Ioerr) 
C
      If (Ioerr .ne. 0) then
         Ierr = 7
         Write(Luout,*) '@OPnfil-I I/O error on Unit ', Lun,'.'
         Return
      Endif
C
      Rewind (Lun)
C
      Return
      End
C      
      Subroutine Rdvprp(Nbasis, Lun, Filnam, Label, Dipole, Ierr, Buf,
     &                  Ibuf, Libuf)
C
C Description
C
C   Read specified property integral matrix from the disk.
C   Because integrals may break the symmetry of the molecule,
C   they are returned in full square-matrix.
C
C Variables
C
C   Nbasis  : Number of basis function
C   Lun     : Unit number corresponds to vprop integral file (input)
C   Filnam  : File name for property integrals (input)
C   Dipole  : Area to store dipole integrals from vprop (output)
C   Label   : Vprop record label for properties (input)
C   Ierr    : Return code
C             = 0 Succesful completion
C             = 1 unable to open/read file
C   Nbuf    : Number of integral in Vprop integral buffer
C   Buf     : Buffer for reading integrals (scratch)
C   Ibuf    : Buffer for reading integrals (scartch)
C
C Declarations
C
      Implicit Double Precision (A-H, O-Z)
      Character*(*) Filnam, Label
      Dimension Ibuf(Libuf), Buf(Libuf), Dipole(Nbasis,Nbasis)
C
C Common Block informations
C
      Common/Machsp/Iintln,Ifltln,Iintfp,Ialone,Ibitwd
      Common/Files/Luout,Moints

C Open the integral File
 
      Call Opnfil(Lun, Filnam, 'UNFORMATTED', 0, Ierr)
      If (Mod(Ierr, 2) .ne. 0) then
         Write (Luout, *) ' @Rdvprp-I, Unable to open file ',Filnam,'.'
         Return
      endif
C 
C Find the pertinent integrals
C
      Rewind (Lun)
      Call Filabl(Lun, Label, Ierr)
      If (Ierr .eq. 1) then
         Write(Luout, 1000) Label, Filnam, Lun
         Return
      Endif
 1000 Format (1X, ' @Rdvprp-I, Label,',A,', not found on file ', A,
     &        ' (unit ', I2,').')
C
      Call zero(Dipole, Nbasis*Nbasis)

C     Read the vprop integral file

 200  Read (Lun) Buf, Ibuf, Nbuf
      If (Nbuf .ne. -1) then
C
C Map the indices from triangular to square
C
         Do 10 J = 1, Nbuf
C              
C Ibuf carry the index corresponding to elements of the square matrix
C So m*(m-1)/2 is greater than Ibuf(J). 

            m = 1 + (-1 + INT(Dsqrt(8.D0*Ibuf(J)+0.999D0)))/2
C
C Now fairly easily n 
C
            n = Ibuf(J) - (m*(m-1))/2
            
            Dipole(n, m) = Buf(J)
            Dipole(m, n) = Buf(J)
 10      continue
         Goto 200
      Endif
         
C
      Close(Lun)
      Return
      End
C
C $Header: /home/orange/2/qtp/rjb/rp/lib/util/RCS/tab.f,v 1.2 90/09/11 12:22:20 peloquin Exp Locker: peloquin $
C $Log:	tab.f,v $
c Revision 1.2  90/09/11  12:22:20  peloquin
c adjusted format for column headers
c 
c Revision 1.1  89/07/14  10:26:42  sosa
c Initial revision
c 
C....+!..1.........2.........3.........4.........5.........6.........7.*
      SUBROUTINE TAB(NOUT,A,N,M,NNN,MMM)
C  modified by RPM for an 80-character wide display.  un-comment the 
C  appropriate set of parameters and format statements to choose how
C  many columns across the page.
      IMPLICIT REAL*8(A-H,O-Z)
C ***
      DIMENSION A(NNN,MMM)
       integer wid
C WID 120
C      parameter(wid = 10)
C    1 FORMAT(1X,I3,3X,10F12.6)
C   12 FORMAT(4X,10(8X,I4))
C wid 80
C       parameter(wid=5)
C    1 FORMAT(1X,I3,3X,5e14.6)
C12    FORMAT(4X, 5(8x,i4))
C      parameter (wid = 4)
C   1 FORMAT(1X,I3,1X,4e18.10)
C12    FORMAT(4X, 4(8x,i4))
       parameter(wid=6)
 1     FORMAT(1X,I3,2X,6F12.6)
 12    FORMAT(4X, 6(8x,i4))
C ***
   11 FORMAT(/)
      MM=M/wid
      IF(MM.EQ.0) GO TO 6
          DO 4 II=1,MM
             JP=(II-1)*wid+1
             JK=II*wid
             WRITE(nout,11)
             WRITE(nout,12)(I,I=JP,JK)
             DO 4 I=1,N
                WRITE(nout,1)I,(A(I,J),J=JP,JK)
    4     CONTINUE
    6 CONTINUE
      MA=MM*wid+1
      IF(MA .GT. M) RETURN
      WRITE(nout,11)
      WRITE(nout,12)(I,I=MA,M)
      DO 5 I=1,N
          WRITE(nout,1) I,(A(I,J),J=MA,M)
    5 CONTINUE
      RETURN
      END


