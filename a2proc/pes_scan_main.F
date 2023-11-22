#include <flags.h>
      Subroutine Pes_scan_main(Streaming,Stationary,Drive_IRC,Prop)
      Implicit Double Precision (A-H, O-Z) 
C
#include <icore.com>
#include <machsp.com>
#include <mxatms.par>
#include <jodaflags.com>
#include <const.h>
C
      Double Precision Norm_coords, Omega, Lamsearch_Tol, 
     &                 Ln_Intrp_Tol
      Character*80 Fname, Direction
      Character*12 Derivatives 
      Character*80 tmp
      Character*10 Prop
      Character*5 AtmLabel(Mxatms), Constr_Min
      Character*4 SymLabel(3*Mxatms), Units
      Character*12 Grid_file
      Character*11 Vib_Type(3*Mxatms)
      Character*5 Label(Mxatms) 

C
C     Conversion factors from bohr to angstrom and atomic units to cm-1
C
      Parameter (B2ang=CONV_a0_a, Au2Invcm=5.14048D03)
C
      Logical ZMT_Present, GBS_Present, JARC_Present, JNDX_Present,
     &        Can_Do_Aces, Ref_done, User_input, Direct, 
     &        Grid_file_Present, Drive_IRC, Trans_state, 
     &        Get_Hess, More_pts_left, End_point, Begin_IRC,
     &        Forward, Get_grad, Mass_Weigh_nm, Mass_Weigh_gr,
     &        New_IRC, Look4_IRC,Print_header,Do_next,Abort 

      Logical FlgACESGeom,FlgACESElec,FlgACESGrad,FlgACESHess,
     &        FRQARC_EXIST, Streaming, Stationary, Wrt_extrnls
      Integer Iatmchrg(Mxatms)
C 
      Dimension VCoords(Mxatms*3), A2grad(Mxatms*3), Rgrad(Mxatms*3),
     &          A2hess(Mxatms*Mxatms*9), Imap(Mxatms), IOrd(Mxatms),
     &          AtmMass(Mxatms), Rhess(Mxatms*Mxatms*9),
     &          Norm_coords(9*Mxatms*Mxatms), Omega(3*Mxatms),
     &          Vomega(3*Mxatms),Btmp(3*Mxatms), Deltaq(3*Mxatms),
     &          Deltax(3*Mxatms), Coords(Mxatms*3), 
     &          Coords_K0O(Mxatms*3), Coords_K1C(Mxatms*3),  
     &          Coords_K1P(Mxatms*3), Coords_K1CM(Mxatms*3),
     &          Work(50*Mxatms*Mxatms),Grad_K0O(3*Mxatms),
     &          Vec_K1C0(Mxatms*3), Grad_on_K0O(Mxatms*3), 
     &          Grad_on_K1C(Mxatms*3), Vec_K0O(Mxatms*3), 
     &          Steep_Dpath(3*Mxatms), Grad_K1CM(3*Mxatms),
     &          Vec_K1C1(Mxatms*3),Grad_K0OM(3*Mxatms), Grad_stat(6)
      Data Inull,Dnull,Ione /0,0.0D0,1/
c
c First check whether JOBARC/JAINDX present in the current directory.
c
      Call Gfname('ZMAT',   Fname, Length)
      Inquire(FILE=Fname(1:Length), EXIST=ZMT_Present)
      Call Gfname('GENBAS', Fname, Length)
      Inquire(FILE=Fname(1:Length), EXIST=GBS_Present)
  
      Call Gfname('JOBARC', Fname, Length)
      Inquire(FILE=Fname(1:Length), EXIST=JARC_Present) 
      Call Gfname('JAINDX', Fname, Length)
      Inquire(FILE=Fname(1:Length), EXIST=JNDX_Present)
      If (ZMT_Present .And. GBS_Present)   Can_Do_Aces = .True.
      If (JARC_Present .And. JNDX_Present) Ref_done    = .True. 
c 
      FlgACESGeom = .TRUE.
      Wrt_Extrnls = .TRUE.
      FlgACESElec = .FALSE.
      FlgACESGrad = .FALSE.
      FlgACESHess = .FALSE.

      Maxmem_Stack = 3*Mxatms*25* + 9*Mxatms*Mxatms*3 
     &               + Mxatms*Mxatms*50
      Maxmem_Stack = 8*Maxmem_Stack/1.0D06
      Print*, Maxmem_Stack

      Write(6,"(a,ES12.6E2,a)") " This code assumes that there is ",
     &                            Dble(Maxmem_Stack),
     &                          "  of Mb stack memory."
C
      If (Can_Do_Aces .And. .Not. Ref_done) Then
c 
c There is a reference ZMAT and GENBAS (minimum requiremnt for 
c an ACES calculations) Lets run the reference point. For PES scans
c on a user defined grid, the reference run is a single point 
c energy or gradient calculation. 
c
         Call Runit("runaces2a")
         call aces_init_rte
         call aces_com_parallel_aces
         call aces_ja_init
C
      Else if (Ref_done) Then
C
         Call aces_init_rte
         Call aces_com_parallel_aces
         Call aces_ja_init
C
      Endif 
C
      call getrec(1,'JOBARC','IFLAGS',   100,iflags)
      call getrec(1,'JOBARC','IFLAGS2',  500,iflags2)
      call getrec(1,'JOBARC','NREALATM', 1,  Nreals)
      Call Getrec(1,'JOBARC','ZMATATMS', 1,  Natoms)
      Call Getrec(1,'JOBARC','ATOMMASS', Natoms*IINTFP,AtmMass)
      Call Getrec(1,'JOBARC','UHFRHF  ', 1,iUHF)
      Call Getrec(1,'JOBARC','NIRREP  ', 1,NIrRep)
      Call Getrec(1,'JOBARC','MAP2ZMAT', Natoms,Imap)
      Call Getrec(1, 'JOBARC','LINEAR  ', 1, Iamlinear)
      Call Getcrec(1,'JOBARC','ZSYM',5*Natoms,AtmLabel)
      Call Getrec(1, 'JOBARC', 'ATOMCHRG', Natoms, Iatmchrg)
C
      If (Iamlinear .EQ. 1) Then
         Nvibs = 3*Nreals - 5
         Icol  = 5
      Else
         Nvibs = 3*Nreals - 6 
         Icol  = 6
      Endif

      Ivib_level = Iflags(h_IFLAGS_vib) 
      Igrad_calc = Iflags2(h_IFLAGS2_grad_calc)

C Ivib_level = 1 (analytic Hessian)
C Ivib_level = 3 (finite difference) (2 is obsoleted)
C Igrad_calc = 1 (analytical)
C Igrad_calc = 2 (numerical)
c
CSSS      Call a2_reset_jarc()
C
      If (streaming) Then
C
C This is for direct dynamics runs. The reference point has alredy
C been run and the geometry of the point that need to be done is 
C pass from the dynamics program. 
c
         Call a2_reset_jarc()
C
#ifdef _DEBUG_LVL0
         Call Getrec(20, "JOBARC", "REF_GEOM", Natoms*3*IINTFP, 
     &               Coords)
         Write(6, "(3F10.5)") (Coords(i),i=1,3*Natoms)
C
         Ioff = 1
         Do Iatm = 1, 2
            Do Ixyz = 1, 3
               Coords(ioff) = (Coords(ioff)) + 0.01
               ioff = ioff + 1
            Enddo
         Enddo
         Write(6,*)
         Write(6, "(3F10.5)") (Coords(i),i=1,3*Natoms)
         Wrt_Extrnls = .True.
#endif
C
C Write the new incomming geometry. The next calculation is performed
C at this geometry.
C
         Call Putrec(20, "JOBARC", "COORD  ", Natoms*3*IINTFP,
     &               Coords)
         Call aces_ja_fin
C
C Run the point and read the energy and the gradients and Hessians
C if available. 
C
         Call Runit("runaces2b")
         Call aces_ja_init
C
#ifdef _DYNAMIC_MEMORY
         icrsiz = iflags(h_IFLAGS_mem)
         icore(1) = 0
         do while ((icore(1).eq.0).and.(icrsiz.gt.1000000))
            call aces_malloc(icrsiz,icore,i0)
            if (icore(1).eq.0) icrsiz = icrsiz - 1000000
         end do
         if (icore(1).eq.0) then
             print *, '@MAIN: unable to allocate at least ',
     &                 1000000,' integers of memory'
            call aces_exit(1)
         end if
         iatmMass = iocc + 16
         icoord   = iatmMass + 3*Natoms*iintfp
         ivgrad   = icoord   + 3*NAtoms*iintfp
         izgrad   = ivgrad   + 3*Nreals*iintfp
         ihess    = izgrad   + 3*Nreals*iintfp
         ifreq    = ihess    + 9*Natoms*Natoms*iintfp
         inex     = ifreq    + 3*Nreals*iintfp
#endif
C
         Call Get_flags(FlgACESElec, FlgACESGrad,  FlgACESHess,
     &                  FlgACESForc)
C
         Call A2read_Jarc(FlgACESGeom, FlgACESElec, FlgACESGrad,
     &                    FlgACESHess, FlgACESForc, Natoms, Nreals, 
     &                    Imap, Coords, A2Ener, AtmMass, AtmLabel, 
     &                    Rgrad, A2grad, Rhess, A2hess, Vomega, 
     &                    Omega, Iuhf, Nirrep, Nvibs, B2ang, 
     &                    Au2Invcm, Wrt_extrnls)
C
         Call a2_reset_jarc()
C
      Else if (Stationary) Then 
C
C This is primarily reserved for running single points calculations 
C on a grid. The grid can be given in a file or can be genrerated.
C The GIRD file (if not automatic) must be formatted file of Cartesian 
C coordinates. 
C
         Call A2get_gridnmlist(Stride,Derivatives,Grid_file_present,
     &                         Npoints) 

C 1. GIRD can be user defined. If user defined, give the number
C    numbers of atoms (as in the reference ZMAT) whose coordinates
C    are changing (others are kept at the reference values).
C
C 2. Automatically generated on normal modes. The refence point
C    must be a frequency calculations. 
C
         If (Grid_file_present) Then

            Inquire(FILE=Grid_file,EXIST=Grid_file_Present)

            If (Grid_file_Present) Then
                Open(Unit=13, File=Grid_file, Form="FORMATTED",
     &               Status="OLD")
            Else
                Write(6, "(a,1x,a,a)") "The file", Grid_file, 
     &                                     "does not exsist." 
                Call Errex
            Endif 
C
            Call a2_reset_jarc()
            Do Ipts = 1, Nmbrof_gridpts 
C
C Get the point and write the JOBARC record
C
               Read(13, *, End=20)
               Do Iatm = 1, Natoms
                  Ioff = 3*(Iatm - 1)
                  Read(13,*,End=20) (Coords(Ioff + I), I=1, 3)
               Enddo
C       
#ifdef _DEBUG_LVL0
                Write(6,*)
                Write(6, "(a)") "Grid points read from Grid_file"
                Write(6, "(3F10.5)") (Coords(i),i=1,3*Natoms)
#endif
    
               Call Putrec(20, "JOBARC", "COORD  ", Natoms*3*IINTFP, 
     &                    Coords)
               Call aces_ja_fin
C
C Run the point and read the energy and the gradients 
C
               Call Runit("runaces2b")
               Call aces_ja_init
C
               Call Get_flags(FlgACESElec, FlgACESGrad,  FlgACESHess,
     &                        FlgACESForc)
C
C Extract the energy, gradients and Hessians (if they are in JOBARC).
C
               Call A2read_Jarc(FlgACESGeom, FlgACESElec, 
     &                          FlgACESGrad, FlgACESHess, 
     &                          FlgACESForc, Natoms, 
     &                          Nreals, Imap, Coords, A2Ener, 
     &                          AtmMass, AtmLabel, Rgrad, A2grad, 
     &                          Rhess, A2hess, Vomega, Omega, Iuhf, 
     &                          Nirrep, Nvibs, B2ang, Au2Invcm, 
     &                          Wrt_extrnls)
               Call a2_reset_jarc()
            Enddo
  20        Close(13)
C
         Else 
C
C This block handles scan along normal modes. The requirement is that
C refrence job is a vibrational frequency calculation. The displacement
C size and the length must be provided in PES_SCAN namelist. Note: 
C If Mass_weigh_nm False return normal modes transformation in L dimensions
C (instead of M^(1/2)L). 
C 
            Trans_state   = .False.
            Mass_weigh_nm = .False. 
            Mass_weigh_gr = .False. 
            Get_Hess      = .False.
            Get_Grad      = .False.
            Call Get_refvib_data(VCoords, Norm_coords, Vomega, Omega,
     &                           AtmMass, AtmLabel, SymLabel, Coords,  
     &                           Btmp, A2grad, Rgrad, A2Hess, Rhess, 
     &                           Imap, Vib_Type, Nreals, Natoms, 
     &                           Nvibs, B2ang, Au2invcm, Trans_state, 
     &                           Mass_Weigh_nm, Mass_Weigh_gr, 
     &                           Get_hess, Get_Grad, Ivib_level,
     &                           Igrad_calc,Icol)

#ifdef _DEBUG_LVL0
            Write(6,*)
            Write(6,*) "The vibrational frequencies and symmetries"
            Write(6, "(4F12.5,1x)") (Omega(I), I = 1, Nvibs)
            Write(6,*)
            Write(6,*) "The reference geometry (in Angstrom)"
            Write(6, "(3F10.5)") (B2ang*Coords(i),i=1,3*Nreals)
            Write(6,*)
            Write(6, "(4a4)") (SymLabel(I), I = 1, Nvibs)
            Write(6,*) "Normal modes"
            Do I = 1, Nvibs
               Ioff = (I-1)*3*Nreals
                  Write(6,"(3F17.13)") (Norm_coords(Ioff + k), k=1,
     &                                 3*Nreals)
                  Write(6,*)
            Enddo
           Write(6,"(a)") " The atomic charges:"
           Write(6,"(10I3)") (Iatmchrg(I), I=1, Natoms)
           Write(6,"(a)") " The atom lables:"
           Write(6,"(10a)") (AtmLabel(I), I=1, Natoms)
           Write(6,"(a)") " The Vmol<->ZMAT:"
           Write(6,"(10I3)") (Imap(I), I=1, Natoms)
#endif
CSSS           Call a2_reset_jarc_0()
C 
C The vibrational SCF, averaging etc and many other post electronic
C structure calcualtions depends on having PES on normal modes.
C
C V(q_1,q_2, ...,q_n) = E_0(q_1(0), q_2(0)..) + Sum_i(DE/dq_i)q_i +
C Sum_ij(DE^2/dq_idq_j)q_iq_j + Sum_ijk(D^3E/dq_idq_jdq_k) 
C q_iq_jq_k + Sum_ijkl(D^4E/dq_idq_jdq_kdq_l)q_iq_jq_kq_l
C
C Check whether there are imaginary modes, if so warn abd quit.

            Ineg = 0
            Do I = 1, Nvibs
               If (Omega(I) .Le. Dnull) Ineg = Ineg + 1
            Enddo 
            If (Ineg .Gt. Inull) Then
               Write(6,"(2a)") " There are imaginary frequencies and",
     &                         " the scan_type=stationary is allowed"
               Write(6,"(a)") " only for minima."
               Call Errex
            Endif 
 
            Delx = Stride 
            Call Stpsize(Omega, Deltaq, Delx, Nvibs, B2Ang)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " Dimensionless step size will be"
      Write(6, "(4F15.5)") (Deltaq(I), I=1, Nvibs)
      Write(6,"(a,1x,12a)") " The derivative type         :", 
     &                        Derivatives
      Write(6,"(a,1x,I3)") " The number of displacements  :",
     &                      Npoints 
      call check_norm(Norm_coords,Work,Nreals,Nvibs)
#endif
C The dimensionless normal modes (q) are related to the mass-weighted normal
C coordinates Q, q = Sqrt(2 pi*c*omega/hbar)Q. Lets denote the constant 
C factor FRED (2pi*c/hbar) then this relation is q = Fred*sqrt(omega)Q 
C (the Fred*sqrt(omega) has 1/(Sqrt(mass)*L) dimensions). 
C The displacements are taken along the dimensionless coords. Therefore
C Delta Q = (1/(FRED*sqrt(Omega)) * T * Delta q. T is the eigenvectors
C of the mass-weighted Hessian multiplied by the 1/sqrt(mass) (therefore 
C T itself is not unitary).The 1/(FRED*Omega) has (Sqrt(mass)*L and 
C the Sqrt(mass) * T convert  it back to the unitarty transformation 
C matrix (i.e. unitary eigenvectors of the mass-weigthed hessian 
C witout any dimensions. This I assume from many numerical experiments
C see below)). So, the dimensions of final Delta Q is Length
C and we can add it to the equilibrium geometry. 
C
C This is my current understanding of this (after many numerical 
C experiments in which I calculate the numerical Cartesian Hessian based on
C the displacments generated as above, mass-weigh and convert to normal 
C coordinate basis. If everything is correct then, I should have a diagonal
C matrix with  square root of the diagoals being the vibrational 
C frequencies. This test is always performed as a check. The question is
C why Q is in dimension of L when it should be M(1/2)L. If the 
C (1/(FRED*sqrt(Omega)) * T * Delta q is M(1/2)L we should not add this 
C to equilibrium geo. without mass-weighing. If I do that above test 
C fails. We can mass-weigh the final Q. That leads to the same problem.
C The mystery seems to reside in what should T. My opinion it should 
C be unitary and have no dimensions. 

            Call Dcopy(3*Nreals,Coords,1,Btmp,1)
CSSS            Call Get_massw_coords(Btmp,AtmMass,Nreals)
            Write(6,"(a)") " The positive displacements" 
            Do Imode = 1, Nvibs
               Ioff = 3*(Imode-1)*Nreals + Ione 
               Do Icalc = 1, Npoints
                  Call Run_a2(Deltaq,Coords,Btmp,Norm_coords(Ioff),
     &                        Iatmchrg,AtmMass,Nreals,Natoms,Icalc,
     &                        Imode,Prop,"Plus ",Atmlabel,
     &                        Label,Imap)
               Enddo
            Enddo 

            If (Derivatives .EQ. "Double_sided") Then
               Write(6,"(a)") " The negative displacements" 
               Do Imode = 1, Nvibs
                  Ioff = 3*(Imode-1)*Nreals + Ione 
                  Do Icalc = 1, Npoints
                     Call Run_a2(-Deltaq,Coords,Btmp,Norm_coords(Ioff),
     &                            Iatmchrg,AtmMass,Nreals,Natoms,Icalc,
     &                            Imode,Prop,"Minus",Atmlabel,
     &                            Label,Imap)
                   Enddo
               Enddo 
            Endif 
         Endif
C
      Else If (Drive_IRC) Then
C
C If the starting point is the transition state, the fist pivot point 
C is half distance from the steepest decent path. The reference calculation
C is a vibrational frequnecy calculation of the transition state.
C
C **** Nature of the starting point must be read from the IRC namelist***
C
           Call Read_IRCnmlist(Stride, Trans_state, Direction, 
     &                         Max_IRC_Steps,Max_IRC_Searches,
     &                         Do_next)
           
           Forward = .TRUE.
           If (Direction .EQ. "REVERSE") Forward = .False.

C Note: Mass_weigh_nm False return normal modes in L dimensions (instead
C of M^(1/2)L).

           Mass_weigh_nm = .False. 
           Mass_weigh_gr = .True. 
           Get_Hess      = .True.
           Get_Grad      = .True. 
           Call Get_refvib_data(VCoords, Norm_coords, Vomega, Omega,
     &                          AtmMass, AtmLabel, SymLabel, Coords, 
     &                          Btmp, A2grad, Grad_K1CM, A2hess, Rhess, 
     &                          Imap, Vib_Type, Nreals, Natoms, Nvibs,  
     &                          B2ang, Au2Invcm, Trans_state, 
     &                          Mass_weigh_nm, Mass_weigh_gr, 
     &                          Get_hess, Get_grad, Ivib_level,
     &                          Igrad_calc,Icol)
#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,*) "The Hessian @-PES_scan"
      Call output(Rhess, 1, 3*Nreals, 1, 3*Nreals, 3*Nreals,
     &            3*Nreals, 1)
#endif
C#ifdef _DEBUG_LVL0 
            Write(6,*)             
            Write(6,*) "The vibrational frequencies and symmetries"
            Write(6, "(4F12.5,1x)") (Omega(I), I = 1, Nvibs)
            Write(6,*)
            Write(6, "(4a4)") (SymLabel(I), I = 1, Nvibs)
            Write(6,*)
            Write(6,*) "The reference geometry"
            Write(6, "(3F10.5)") (Coords(i),i=1,3*Natoms)
            Write(6,*)
            Write(6, "(4a4)") (SymLabel(I), I = 1, Nvibs)
            Write(6,*) "Normal modes"
            Do I = 1, Nvibs
               Ioff = (I-1)*3*Nreals 
                  Write(6,"(3F17.13)") (Norm_coords(Ioff + k), k=1, 
     &                                 3*Nreals)
                  Write(6,*)
            Enddo 
           Write(6,"(a)") " The atomic charges:"
           Write(6,"(10I3)") (Iatmchrg(I), I=1, Natoms)
C#endif      
C           
           Call Getrec(20,'JOBARC','TOTENERG', IINTFP, E_refrence)
           Call a2_reset_jarc()
C 
C The eigenvector corresponding to the negative eigenvalue (no longer
C mass weighted). Generate the pivot point: that is 1/2stride of 
C the steepest decent eigenvector (when the hessian is availabe). 
C Otherwise choose the largest gradient 
C
           N_Irc_step    = 0
           More_pts_left = .True.
           End_point     = .False.
           Begin_IRC     = .True.
           Opt_tol       = Iflags2(h_IFLAGS_convergence)
           Max_Cycles    = Max_IRC_Searches 
           Stride_total  = 0.0D0
           Stride_Tol    = 1.0D-04
           Ln_Intrp_Tol  = 5.0D0
           Hess_Eval_Tol = 1.0D-09
           Delta         = 5.0D-02
           lamsearch_TOl = 1.0D-16
           Binsearch_Tol = 1.0D-08
           Small         = 1.0D-05
           G_Cutoff      = 1.0D-04*(Stride/0.30D0)
           R_Cutoff      = 1.0D-03*(Stride/0.30D0)
           Units         = "Bohr"
           Halfs         = 0.50D0*Stride
           Print_header  = .True.
C
           Do while (More_pts_left .AND. .NOT. End_point) 

              N_Irc_step = N_Irc_step + 1 

              If (Begin_IRC) Then
                  Abort         = .False.
                  If (Trans_state) Then
                  If (Direction .EQ. "FORWARD") Then
                  Write(6,"(3a)") " Begining a forward direction IRC",
     &                            " search starting from the saddle",
     &                            " point."
                  Write(6,"(2a)") "          *--------------------",
     &                            "-------------------*"
                  Else
                  Write(6,"(3a)") " Begining a reverse direction IRC",
     &                            " search starting from the saddle",
     &                            " point."
                  Write(6,"(2a)") "          *--------------------",
     &                            "-------------------*"
                  Endif 
                  Endif 
                  Write(6,*) 
                  Write(6,"(2a,i3,a)") " Begining to look for",
     &                                 " the IRC point no.", 
     &                                   N_Irc_step,"."
              Else
                  Abort         = .False.
                  Write(6,*) 
                  Write(6,"(2a,i3,a)") " Begining to look for",
     &                                 " the IRC point no.", 
     &                                   N_Irc_step,"."
                  Write(6,*)
              Endif 
C
              If (Begin_IRC) Then
C
                  If (Trans_state) Then 
                     Call Dcopy(3*Nreals, Norm_Coords, 1, 
     &                         Steep_DPath, 1)
C
                     Dpath_Norm = 0.0D0
                     Ioff       = 1
                     Do Iatom = 1, Nreals
                        Sqrtmass = Dsqrt(AtmMass(Iatom))
                        Do Ixyz = 1, 3
                           Dlen       = Steep_DPath(Ioff)*Sqrtmass
                           Dpath_Norm = Dlen*Dlen + Dpath_Norm 
                           Ioff = Ioff + 1 
                        Enddo
                     Enddo

                     Dpath_Norm = DSqrt(Dpath_Norm)
C
                     If (Dpath_Norm .LT. Small) Then
                         Write(6,"(a,a,a)") "Norm of the ",
     &                                   " mass-weighted steepest",
     &                                   " decent path is too small."
                          Call Errex
                     Endif

C The Cartesian steepest descent path with unit length in mass-weighted 
C Cartesian space. 
                     Call Dscal(3*Nreals, 1.0D0/Dpath_Norm, 
     &                          Steep_DPath, 1)

#ifdef _DEBUG_LVL0   
            Write(6,*) 
            Write(6,"(2a)") "@-PES_scan The steepest decent path",
     &                      " in ordinary Cartesian Coordinates",
     &                      " (not mass weighted!)"
            Write(6, "(3F17.13)") (Steep_DPath(i),i=1,3*Nreals)
#endif
C
                  Else
#ifdef _DEBUG_LVL0
            Write(6,*)
            Write(6,"(2a)") "Mass weigthed Gradients for non saddle",
     &                      " starting point"
            Write(6, "(3F17.13)") (Grad_K1CM(i),i=1,3*Nreals)
#endif
                     Grad_Norm = Ddot(3*Nreals, Grad_K1CM, 1, 
     &                                Grad_K1CM, 1)
C
                     Grad_Norm = Dsqrt(Grad_Norm)

                     If (Grad_norm .LT. Small/3.0D00) Then
                        Write(6,"(a,a,a,/,a)") " Gradient norm is too",
     &                                   " small, can not start from",
     &                                   " a stationary point except",
     &                                   " the transition state."
                         Call Errex
                     Endif 
C
                     Call Dscal(3*Nreals, 0.5D0*Stride/Grad_Norm, 
     &                          Grad_K1CM, 1)
#ifdef _DEBUG_LVL0
            Write(6,"(a,a)") "Normalized and scalled mass weighted ",
     &                      "Gradients for non saddle starting point"
            Write(6, "(3F17.13)") (Grad_K1CM(i),i=1,3*Nreals)
#endif
                    Ioff = 1
                     Do Iatom = 1, Nreals
                        Sqrtmass = Dsqrt(AtmMass(Iatom))
                        Do Ixyz = 1, 3
                           Grad_K0O(Ioff) = Grad_K1CM(Ioff)*Sqrtmass
                           Ioff = Ioff + 1
                        Enddo
                     Enddo
#ifdef _DEBUG_LVL0
            Write(6,"(a,a)") "Normalized and scalled Gradients for ",
     &                       "non saddle starting point (no. M.W)"
            Write(6, "(3F17.13)") (Grad_K0O(i),i=1,3*Nreals)
            Grad0_norm = Ddot(3*nreals,Grad_K0O,1,Grad_K0O,1)
            Write(6, "(a,F17.13)") "Grad0_norm (no M.W):",
     &                              Dsqrt(Grad0_norm)
#endif
C
                  Endif 
C                  
                  Step = 0.5D0*Stride
C
                  Imax = iDAMAX(3*Nreals, Steep_Dpath, 1)
C
                  If (Forward .AND. Steep_Dpath(Imax) .LT. 0.0D0) 
     &               Step = -1.0D0*Step
                  If (.NOT. Forward .AND. Steep_Dpath(Imax) .GT.
     &               0.0D0)
     &               Step = -1.0D0*Step 
C
C _K0O is the pevious point, _K1P is the pivot point and _K1C is the 
C current point.
C
                  Call Dcopy(3*Nreals,        Coords,    1, 
     &                       Coords_K0O, 1)
                  Call Dcopy(3*Nreals,        Coords,    1, 
     &                       Coords_K1P, 1)
                  Call Dcopy(3*Nreals,        Coords,    1, 
     &                       Coords_K1C, 1)
C
                  If (Trans_state) Then
                      Call Daxpy(3*Nreals, Step,   Steep_DPath, 1, 
     &                           Coords_K1P, 1)
                      Call Daxpy(3*Nreals, 2.0D0*step, Steep_DPath,  
     &                           1, Coords_K1C, 1)
                      Call Dcopy(3*Nreals,       Coords_K1C,    1, 
     &                          Coords_K1CM, 1)
                  Else
                      Call Daxpy(3*Nreals, 1.0D0, Grad_K0O, 1, 
     &                           Coords_K1P, 1)
                      Call Daxpy(3*Nreals, 2.0D0, Grad_k0O, 1, 
     &                           Coords_K1C, 1)
                      Call Dcopy(3*Nreals,   Coords_K1C,    1, 
     &                           Coords_K1CM, 1)
                  Endif
 
#ifdef _DEBUG_LVL0
            Write(6,"(2a)") "---Current,Pivot,Starting geoms. (not mass"
     &                      ," weighted)---"
            Write(6,*)
            Write(6,"(a,1x,F10.5)") "The initial displacement:",
     &                               Step
            Call Dscal(3*Nreals,Step,Steep_Dpath,1)
            Write(6,*) "The displacement along normal modes"
            Write(6, "(3F17.13)") (Steep_DPath(i),i=1,3*Natoms)
            Write(6,*) "The current geometry:K1C"
            Write(6, "(3F17.13)") (Coords_K1C(i),i=1,3*Natoms)
            Write(6,*) "The current geometry:K1CM"
            Write(6, "(3F17.13)") (Coords_K1C(i),i=1,3*Natoms)
            Write(6,*) "The pivot point geometry:K1P"
            Write(6, "(3F17.13)") (Coords_K1P(i),i=1,3*Natoms)
            Write(6,*) "The previous geometry:K0O"
            Write(6, "(3F17.13)") (Coords_K0O(i),i=1,3*Natoms)
#endif  
C
                  Ncycles    = 1 
                  N_IRC_Step = 1
                  Constr_Min = "Null "

                  Do While (.NOT. (Constr_Min .EQ. "Found") .AND. 
     &                      .NOT. (Abort)) 

                           Write(6,"(3a,I3)") " No. cycles in",
     &                                        " starting constrained",
     &                                        " search: ",Ncycles

                           Call Constr_search(Coords_K0O, Coords_K1C, 
     &                                        Coords_K1P, Vcoords,
     &                                        Coords, Coords_K1CM, 
     &                                        A2grad, Grad_K1CM,  
     &                                        Grad_K0O, Grad_K0OM,
     &                                        A2Hess, Rhess, AtmMAss, 
     &                                        Vec_K1C0, Vec_K1C1,
     &                                        Vec_K0O, Grad_on_K0O, 
     &                                        Grad_on_K1C,
     &                                        Work, B2ang, Units,
     &                                        Stride, Delta, 
     &                                        Begin_IRC, New_IRC,
     &                                        Look4_IRC, Constr_Min,
     &                                        Imap, Ncycles, Nreals, 
     &                                        Natoms, Opt_Tol,
     &                                        Stride_Total, 
     &                                        Stride_Tol,  
     &                                        Ln_Intrp_Tol, 
     &                                        Hess_Eval_Tol, 
     &                                        Lamsearch_Tol,
     &                                        Binsearch_Tol,
     &                                        G_cutoff, R_Cutoff,
     &                                        E_currnt,Atmlabel,
     &                                        Iatmchrg,Direction,
     &                                        Print_header,
     &                                        Trans_state,
     &                                        Max_cycles,
     &                                        Do_next,Abort)

#ifdef _DEBUG_LVLM
            Write(6,"(2a)") "@-PES_scan: After Constr_search for",
     &                      " the intial point"
            Write(6,*) "The previous geometry:K0O<->K0OM"
            Write(6, "(3F17.13)") (Coords_K0O(i),i=1,3*Nreals)
            Write(6,*) "The current geometry:K1C"
            Write(6, "(3F17.13)") (Coords_K1C(i),i=1,3*Nreals)
            Write(6,*) "The current geometry:K1CM"
            Write(6, "(3F17.13)") (Coords_K1CM(i),i=1,3*Nreals)
            Write(6,*) "The pivot point geometry:K1P<->K1PM"
            Write(6, "(3F17.13)") (Coords_K1P(i),i=1,3*Nreals)
            Write(6,*) "The M.W. current gradient "
            Write(6, "(3F17.13)") (Grad_K1CM(i),i=1,3*Nreals)
            Write(6,*) "The N.M.W. current gradient "
            Write(6, "(3F17.13)") (A2grad(i),i=1,3*Nreals)
#endif 
#ifdef _DEBUG_LVL0
            Write(6, "(a,2(1x,ES16.9E2))") "E_ref.,E_cur. =", 
     &                                      E_refrence,E_currnt
#endif
                           Begin_IRC  = .False.
                           New_IRC    = .False.
                           Look4_IRC  = .True.
                           Ncycles    = Ncycles + 1
                  Enddo 
C
                  If (E_currnt .GT. E_refrence) Then
                      Write(6, "(a,a)") "Energy can not rise during a",
     &                              " downhill IRC search."
                      Call Errex
                  Else
                      E_refrence = E_currnt
                  Endif  

                  Begin_IRC  = .False.
                  New_IRC    = .True.
                  Look4_IRC  = .False.
              Else
C
C First check to make sure that we are not at a stationary point.
C
                  Call Vstat(A2grad, Grad_stat, 3*Nreals)

                  If (Grad_stat(1) .LT. Opt_tol .AND. 
     &                Grad_stat(5) .LT. Opt_tol/3.0D00) Then
                      Write(6,"(a,a,a,a)")"The RMS Gradient is below",
     &                                     "the tolerence. The last ",
     &                                     "IRC point is closer to  ",
     &                                     "a stationary point"
                        
                      End_point = .TRUE.
                  End if 
C 
C Generate the new pivot point using the gradient. Note that only the
C first step we have an exact Hessian. 
C 
                  Call Dcopy(3*Nreals, Coords_K1C, 1, Coords_K0O, 1)

                  Grad_length = Ddot(3*Nreals, Grad_K1CM, 1, 
     &                               Grad_K1CM, 1)
                  Call Dcopy(3*Nreals, Grad_K1CM, 1, Work, 1)
                  Call Dscal(3*Nreals, 1.0D0/Dsqrt(Grad_length), 
     &                       Work, 1)
#ifdef _DEBUG_LVL0
           Write(6,*) "@-PES_scan: The normalized M.W. current gradient"
           Write(6, "(3F17.13)") (Work(i),i=1,3*Nreals)
#endif
                  Call Dscal(3*Nreals, -Halfs, Work, 1)

#ifdef _DEBUG_LVL0
           Write(6,"(2a)") "@-PES_scan: The normalized M.W. current",
     &                     " gradient*-1/2s"
           Write(6, "(3F17.13)") (Work(i),i=1,3*Nreals)
#endif
                  Ioff = 1
                  Do Iatom = 1, Nreals
                     Sqrtmass = Dsqrt(AtmMass(Iatom)) 
                     Do Ixyz = 1, 3
                        Work(Ioff) = Work(Ioff)/Sqrtmass
                        Ioff = Ioff + 1
                     Enddo
                  Enddo

                  Call Dcopy(3*Nreals, Coords_K1C,  1, Coords_K1P, 1)
                  Call Daxpy(3*Nreals, 1.0D0, Work, 1, Coords_K1P, 1)
                  Call Dcopy(3*Nreals, Coords_K1P,  1, Coords_K1C, 1)
                  Call Daxpy(3*Nreals, 1.0D0, Work, 1, Coords_K1C, 1)
                  Call Dcopy(3*Nreals, Coords_K1C,  1, Coords_K1CM,1)
C
#ifdef _DEBUG_LVL0
            Write(6,*) "@-PES_scan: The current geometry"
            Write(6, "(3F17.13)") (Coords_K1C(i),i=1,3*Nreals)
            Write(6,*) "@-PES_scan: The pivot point geometry"
            Write(6, "(3F17.13)") (Coords_K1P(i),i=1,3*Nreals)
            Write(6,*) "@-PES_scan: previous geometry"
            Write(6, "(3F17.13)") (Coords_K0O(i),i=1,3*Nreals)
            Write(6,"(a,I2)") "The maximum cycles:", Max_Cycles
#endif
                  Ncycles    = 1
                  Constr_Min = "Null "

                  Do While (.NOT. (Constr_Min .EQ. "Found") .AND.
     &                      .NOT. (Abort)) 

                           Write(6,"(3a,I3)") " No. cycles in",
     &                                        " succeding constrained",
     &                                        " searchs: ", Ncycles
C
                           Call Constr_search(Coords_K0O, Coords_K1C,
     &                                        Coords_K1P, Vcoords,
     &                                        Coords, Coords_K1CM,
     &                                        A2grad, Grad_K1CM,
     &                                        Grad_K0O, Grad_K0OM, 
     &                                        A2Hess, Rhess, AtmMAss,
     &                                        Vec_K1C0, Vec_K1C1, 
     &                                        Vec_K0O, Grad_on_K0O,
     &                                        Grad_on_K1C,
     &                                        Work, B2ang, Units,  
     &                                        Stride, Delta, 
     &                                        Begin_IRC, New_IRC,
     &                                        Look4_IRC, Constr_Min, 
     &                                        Imap, Ncycles, Nreals,
     &                                        Natoms, Opt_Tol,
     &                                        Stride_Total,
     &                                        Stride_Tol,
     &                                        Ln_Intrp_Tol,
     &                                        Hess_Eval_Tol,
     &                                        Lamsearch_Tol,
     &                                        Binsearch_Tol,
     &                                        G_Cutoff, R_Cutoff,
     &                                        E_currnt,Atmlabel,
     &                                        Iatmchrg,Direction,
     &                                        .False.,.False.,
     &                                        Max_cycles,
     &                                        Do_next,Abort)
C
                           Begin_IRC  = .False.
                           New_IRC    = .False.
                           Look4_IRC  = .True.
                           Ncycles    = Ncycles + 1
                  Enddo
#ifdef _DEBUG_LVL0
            Write(6, "(a,2(1x,ES13.6E2))") "E_ref.,E_cur. =", 
     &                                      E_refrence,
     &                                      E_currnt
#endif
C
                  If (E_currnt .GT. E_refrence) Then
                      Write(6, "(a,a)") "Energy can not rise during a",
     &                              " downhill IRC search."
                      Call Errex
CSSS                  Else
CSSS                      E_refrence = E_currnt
                  Endif
C
                  Begin_IRC  = .False.
                  New_IRC    = .True.
                  Look4_IRC  = .False.      
C
                  If (N_IRC_Step .GE. Max_IRC_Steps) 
     &                More_pts_left = .False.
C       
              Endif
C
           Enddo                     
           Write(6,*) 
           Write(6,"(a,I3,a)") " The IRC calculation completed the", 
     &                           N_IRC_Step," requested steps." 
     &                     
           Write(6,"(a)") "          *------------------------*" 
C
      Endif
C
      Return
      End

