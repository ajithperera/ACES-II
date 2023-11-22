#include "flags.h"
      Subroutine Vibavg_driver(Work,Maxcor,Iuhf)

      Implicit Double Precision(A-H,O-Z)
       
      Dimension Work(Maxcor)

#include "machsp.com"
#include "mxatms.par"
#include "syminf.com"
#include "sym.com"
#include "flags.com"
#include "flags2.com"
#include "sympop.com"
#include "info.com"
#include "const.h"

C Conversion factors for from bohr to angstrom and atomic units to cm-1

      Parameter (B2ang=CONV_a0_a,Au2Invcm=5.14048D03)
      Parameter (Mxdata = 100)

      Character*6 AtmLabel(Mxatms)
      Character*6 AtmLabel_vmlord(Mxatms)
      Character*6 AtmLabel_zmtord(Mxatms)
      Character*4 SymLabel(3*Mxatms)
      Character*80 Fname 
      Character*12 Derivatives 
      Character*11 Vib_Type(3*Mxatms)
      Character*5 Tmp1(Mxatms)
      Character*5 Tmp2(Mxatms)
      Character*5 Tmp3(Mxatms)

      Dimension AtmMass(Mxatms)
      Dimension AtmMass_zmtord(Mxatms)
      Dimension Imap(Mxatms)
      Dimension Isotopes(Mxatms)
      Dimension Iatmchrg(Mxatms)
      Dimension Iatmchrg_zmtord(Mxatms)
      Dimension Deltaq (3*Mxatms)
      Dimension Temp(Mxdata)
      Dimension Press(Mxdata)
      Dimension RP(3,3)
      Dimension RM(3,3)

      Logical PJJ_SPNCC_present
      Logical MJJ_SPNCC_present
      Logical EJJ_SPNCC_present
      Logical PNMR_shifts_present
      Logical MNMR_shifts_present
      Logical ENMR_shifts_present
      Logical PFORCE_Present
      Logical MFORCE_Present
      Logical EFORCE_Present
      Logical PGRAD_Present
      Logical MGRAD_Present
      Logical EGRAD_Present
      logical PREORI_Present
      Logical MREORI_Present
      Logical ISOTOPES_Present 
      LOgical STEP_LENGTH_present
      Logical NMR_SPINCC
      Logical NMR_SHIFTS
      Logical FORCE
      Logical GRAD
      Logical Grid_file
      Logical Thermo_exist
      Logical Trans_state
      Logical Mass_weigh_nm
      Logical Mass_weigh_gr
      Logical Get_hess
      Logical Get_grad
      Logical Jcc,Sigma

      Data Ione,Inull,I1,I9,I8,I5,I7,It,Ie,Iw,Ih/1,0,1,9,8,5,7,10,11,
     +                                         12,13/

      Call Getrec(20,'JOBARC','NREALATM',Ione,Nreals)
      Call Getrec(20,'JOBARC','ZMATATMS',Ione,Natoms)
      Call Getrec(20,'JOBARC','ATOMMASS',Natoms*IINTFP,AtmMass)
      Call Getrec(20,'JOBARC','MAP2ZMAT',Natoms,Imap)
      Call Getrec(20,'JOBARC','LINEAR  ',Ione,Iamlinear)
      Call Getrec(20,'JOBARC','ATOMCHRG',Natoms,Iatmchrg)
      Call Getcrec(20,'JOBARC','ZSYM',5*Natoms,AtmLabel)

      If (Iamlinear .EQ. 1) Then
         Nmodes = 3*Nreals - 5
         Icol   = 5
      Else
         Nmodes = 3*Nreals - 6
         Icol   = 6 
      Endif

      Ivib_level = Iflags(h_IFLAGS_vib)
      Igrad_calc = Iflags2(h_IFLAGS2_grad_calc)

      Call Gfname("THERMO  ", FNAME, ILENGTH)
      Inquire(FILE=FNAME(1:ILENGTH), EXIST=Thermo_exist)
      If (Thermo_exist) Then
         Open(UNIT=I_inunit, File=FNAME(1:ILENGTH), FORM='FORMATTED')
         Read (I_inunit, *, END=9) Ndata
         If (Ndata .GT. Mxdata) Then
             Write(6,"(2a,I2)") "The number of currently allowed",
     &                          " data points are: ", Mxdata
             Call Errex
         Endif
         Do Idata = 1, Ndata
            Read (I_inunit,*,END=9) Temp(Idata), Press(Idata)
         Enddo
      Else
         Ndata    = 1
         Temp(1)  = 298.15
         Press(1) = 1
      Endif
   9  Close (I_inunit)

C Read the property storage files for displacements and take the derivatives. 
C The data is generated from ACES II runs. The python scripts reads the 
C individual outputs and prepare two files: one for plus displacement
C and one for minus displacement.

      Nreals2 = Nreals*Nreals

      Iunitst = I1
      Iunitp1 = I9
      Iunitm1 = I8
      Iunitp2 = I5
      Iunitm2 = I7
      Iuniteq = It
      Iunitp3 = Ie
      Iunitm3 = Iw
      Iunitis = Ih

      Call Gfname("STEP_LENGT",Fname,Length)
      Inquire(File=Fname(1:Length),Exist=STEP_LENGTH_present)
      If (STEP_LENGTH_present) Then
         Open(Unit=Iunitst,File=Fname(1:Length),Form="Formatted")
         Read(Iunitst,*) Fname 
         Read(Iunitst,"(5(1x,F15.9))") (Deltaq(I),I=1,Nmodes)
      Else 
         Write(6,"(a,a)") " STEP_LENGTH file, which record curveture",
     +                    " independent step lengths is needed."
CSSS         Call Errex
      Endif 
      Close(Iunitst)

      Call A2get_gridnmlist(Stride,Derivatives,Grid_file,
     +                      Nfdif_points)

      If (Derivatives .Eq. "Single_sided") Then
          Write(6,"(a)") " The single sided numerical derivatives",
     +                   " are not advisable"
          Call Errex
      Endif

      Call Gfname("pJJ_SPNCC",Fname,Length)
      Inquire(File=Fname(1:Length),Exist=PJJ_SPNCC_present)
      If (PJJ_SPNCC_present) Open(Unit=Iunitp1,File=Fname(1:Length),
     +                           Form="Formatted")
      Call Gfname("mJJ_SPNCC",Fname,Length)
      Inquire(File=Fname(1:Length),Exist=MJJ_SPNCC_present)
      If (MJJ_SPNCC_present) Open(Unit=Iunitm1,File=Fname(1:Length),
     +                           Form="Formatted")
      Call Gfname("eJJ_SPNCC",Fname,Length)
      Inquire(File=Fname(1:Length),Exist=EJJ_SPNCC_present)
      If (EJJ_SPNCC_present) Open(Unit=Iuniteq,File=Fname(1:Length),
     +                           Form="Formatted")
      If (PJJ_SPNCC_present .or. MJJ_SPNCC_present .or.
     +    EJJ_SPNCC_present) 
     +    NMR_spincc = .True.

#ifdef _DEBUG_LVL0
      write(6,"(a)") "Curveture independent step-lengths"
      Write(6,"(5(1x,F15.9))") (Deltaq(I),I=1,Nmodes)
      Write(6,*)
#endif 

      Iend = Ione 
      If (NMR_spincc) Then

         Iequi_spncc = Iend 
         I1dif_spncc = Iequi_spncc + Nreals2 
         I2dif_spncc = I1dif_spncc + Nreals2*Nmodes 
         Iend        = I2dif_spncc + Nreals2*Nmodes

         Iplus  = Iend 
         Imins  = Iplus + Nreals2*Nfdif_points
         Ichek  = Imins + Nreals2*Nfdif_points 
         
         Ndims        = Nreals 

         If (EJJ_SPNCC_present) Then
#ifdef _DEBUG_LVL0
         Write(6,"(a)") "Getting NMR spin-spinc coupling at equlibrium"
#endif 
            Call Dzero(Work(Iequi_spncc),Nreals2)
            Call Get_prop_values(Work(Iequi_spncc),NMR_spincc,.False.,
     +                           Ndims,Iuniteq)
         Endif 

         Do Imode = 1, Nmodes
            Iplus_off    = Iplus
            Imins_off    = Imins
            Do J = 1, Nfdif_points
               Call Dzero(Work(Iplus_off),Nreals2)
               Call Dzero(Work(Imins_off),Nreals2)
#ifdef _DEBUG_LVLM
               Write(6,"(a)") "Getting JCCs from plus side"
#endif 
               Call Get_prop_values(Work(Iplus_off),NMR_spincc,
     +                              .False.,Ndims,Iunitp1)
#ifdef _DEBUG_LVLM
               Write(6,"(a)") "Getting JCCs from minus side"
#endif 
               Call Get_prop_values(Work(Imins_off),NMR_spincc,
     +                              .False.,Ndims,Iunitm1)
               Iplus_off  = Iplus_off + Nreals2 
               Imins_off  = Imins_off + Nreals2 
            Enddo 

            Step_size = Deltaq(Imode)
            Itype     = 1
            
#ifdef _DEBUG_LVL0
            Write(6,"(a)") "Obtaining the 1st JCC derivatives"
#endif 
C Obtain the first derivatives

            Call Form_2dprop_deriv(Work(Iplus),Work(Imins), 
     +                             Work(I1dif_spncc),Ndims,Nmodes,
     +                             Step_size,1,Nfdif_points,Imode,
     +                             Itype)

C If needed obtain the second derivatives

            If (Nfdif_points .Eq. 2) Then
               Call Form_2dprop_deriv(Work(Iplus),Work(Imins), 
     +                                Work(I2dif_spncc),Ndims,Nmodes,
     +                                Step_size,2,Nfdif_points,Imode,
     +                                Itype)
            Elseif (Nfdif_points .Eq. 1 .And. EJJ_SPNCC_present) Then
               Call Form_2dprop_deriv_eq(Work(Iplus),Work(Imins),
     +                                   Work(Iequi_spncc),
     +                                   Work(I2dif_spncc),Ndims,Nmodes,
     +                                   Step_size,2,Nfdif_points,Imode,
     +                                   Itype)
            Endif 
         Enddo
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(2a)") "The first and second derivatives of the",
     +               " JCC matrices"
      Ioff = 0
      Do Imode = 1, Nmodes
         Call output(Work(I1dif_spncc+Ioff),1,Ndims,1,Ndims,Ndims,
     +               Ndims,1)
         If (Nfdif_points .Eq. 2 .Or. EJJ_SPNCC_present) then
         Call output(Work(I2dif_spncc+Ioff),1,Ndims,1,Ndims,Ndims,
     +               Ndims,1)
         Ioff = Ioff + Nreals2
         endif 
      Enddo 
      Write(6,*)
#endif 
      Endif 

      Close(Iuniteq)

      Call Gfname("pNMR_SHIFT",Fname,Length)
      Inquire(File=Fname(1:Length),Exist=PNMR_SHIFTS_present)
      If (PNMR_SHIFTS_present) Open(Unit=Iunitp2,File=Fname(1:Length),
     +                              Form="Formatted")
      Call Gfname("mNMR_SHIFT",Fname,Length)
      Inquire(File=Fname(1:Length),Exist=MNMR_SHIFTS_present)
      If (MNMR_SHIFTS_present) Open(Unit=Iunitm2,File=Fname(1:Length),
     +                              Form="Formatted")
      Call Gfname("eNMR_SHIFT",Fname,Length)
      Inquire(File=Fname(1:Length),Exist=ENMR_SHIFTS_present)
      If (ENMR_SHIFTS_present) Open(Unit=Iuniteq,File=Fname(1:Length),
     +                              Form="Formatted")
      If (PNMR_SHIFTS_present .or. MNMR_SHIFTS_present .or.
     +    ENMR_SHIFTS_present) NMR_shifts = .True.
       
      If (NMR_shifts) Then
            
         Iequi_shift = Iend
         I1dif_shift = Iequi_shift + Nreals
         I2dif_shift = I1dif_shift + Nreals*Nmodes
         Iend        = I2dif_shift + Nreals*Nmodes 

         Iplus = Iend 
         Imins = Iplus + Nreals*Nfdif_points
         Ichek = Imins + Nreals*Nfdif_points
         Ndims = Nreals

         If (ENMR_SHIFTS_present) Then
#ifdef _DEBUG_LVL0
            Write(6,"(a)") "Getting NMR chemical shifts at equlibrium"
#endif 
            Call Dzero(Work(Iequi_shift),Nreals)
            Call Get_prop_values(Work(Iequi_shift),.False.,NMR_shifts,
     +                                Ndims,Iuniteq)
         Endif 
         Do Imode = 1, Nmodes
            Iplus_off    = Iplus
            Imins_off    = Imins
            Do J = 1, Nfdif_points
               Call Dzero(Work(Iplus_off),Nreals)
               Call Dzero(Work(Imins_off),Nreals)
#ifdef _DEBUG_LVLM
               Write(6,"(a)") "Getting JCCs from plus side"
#endif 
               Call Get_prop_values(Work(Iplus_off),.False.,
     +                              NMR_shifts,Ndims,Iunitp2)
#ifdef _DEBUG_LVLM
               Write(6,"(a)") "Getting JCCs from minus side"
#endif 
               Call Get_prop_values(Work(Imins_off),.False.,
     +                              NMR_shifts,Ndims,Iunitm2)
               Iplus_off  = Iplus_off + Nreals
               Imins_off  = Imins_off + Nreals
            Enddo

            Step_size = Deltaq(Imode)
            Itype     = 1

C Obtain the first derivatives

            Call Form_1dshift_deriv(Work(Iplus),Work(Imins),
     +                              Work(I1dif_shift),Ndims,Nmodes,
     +                              Step_size,1,Nfdif_points,Imode,
     +                              Itype)

C If needed obtain the second derivatives

            If (Nfdif_points .Eq. 2) Then
               Call Form_1dshift_deriv(Work(Iplus),Work(Imins),
     +                                Work(I2dif_shift),Ndims,Nmodes,
     +                                Step_size,2,Nfdif_points,Imode,
     +                                Itype)
            Elseif (Nfdif_points .Eq. 1 .And. ENMR_SHIFTS_present) Then
               Call Form_1dshift_deriv_eq(Work(Iplus),Work(Imins),
     +                                    Work(Iequi_shift),
     +                                    Work(I2dif_shift),Ndims,
     +                                    Nmodes,Step_size,2,
     +                                    Nfdif_points,Imode,Itype)
            Endif 
         Enddo

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(2a)") "The first and second derivatives of the",
     +               " chemical shift vectors"
      Ioff = 0
      Do Imode = 1, Nmodes
         Write(6,"(5(F15.8))") (Work(I1dif_shift+Ioff+i-1),i=1,Ndims)
         Write(6,"(5(F15.8))") (Work(I2dif_shift+Ioff+i-1),i=1,Ndims)
         Ioff = Ioff + Nreals
         Write(6,*) 
      Enddo
      Write(6,*)
#endif
      Endif


      Close(Iunitp1)
      Close(Iunitm1)
      Close(Iunitp2)
      Close(Iunitm2)
      Close(Iuniteq)

C Obtain the cubic force constants. These are obtained as the finite
C difference of force constant matrices. The force-constant matrices
C are obtained fully analytically (CFOUR) or numerical second 
C dervivative of analytic first derivatives.
C The force constant matrices stored in FORCE file (produced by
C python) stored as (3*Nreals,3*Nreals) quantities for no. displacements
C per node basis.

      Call Gfname("pFCM",Fname,Length)
      Inquire(File=Fname(1:Length),Exist=PFORCE_Present)
      If (PFORCE_Present) Open(Unit=Iunitp1,File=Fname(1:Length),
     +                         Form="Formatted")
      Call Gfname("mFCM",Fname,Length)
      Inquire(File=Fname(1:Length),Exist=MFORCE_Present)
      If (MFORCE_Present) Open(Unit=Iunitm1,File=Fname(1:Length),
     +                         Form="Formatted")
      Call Gfname("eFCM",Fname,Length)
      Inquire(File=Fname(1:Length),Exist=EFORCE_Present)
      If (EFORCE_Present) Open(Unit=Iuniteq,File=Fname(1:Length),
     +                         Form="Formatted")
      Call Gfname("pReori",Fname,Length)
      Inquire(File=Fname(1:Length),Exist=PREORI_Present)
      If (PREORI_Present) Open(Unit=Iunitp3,File=Fname(1:Length),
     +                         Form="Formatted")
      Call Gfname("mReori",Fname,Length)
      Inquire(File=Fname(1:Length),Exist=MREORI_Present)
      If (MREORI_Present) Open(Unit=Iunitm3,File=Fname(1:Length),
     +                         Form="Formatted")
      
      If (PFORCE_present .or. MFORCE_present .or. EFORCE_PRESENT) 
     +    FORCE = .True.

      If (FORCE) Then

C Currently We only do first derivatives of the Force connstant matrix.

         Nfdif_points = 1
         Iequi_force  = Iend
         I1dif_force  = Iequi_force + 9*Nreals2
         Iend         = I1dif_force + 9*Nreals2*Nmodes
 
         Iplus  = Iend 
         Imins  = Iplus + 9*Nreals2
         Itemp  = Imins + 9*Nreals2
         Iend   = Itemp + 3*Nreals
         Maxcor = Maxcor - Iend 

         Ndims        = 3*Nreals
         Iplus_off    = Iplus
         Imins_off    = Imins

         If (EFORCE_Present) Then
            Call Dzero(Work(Iequi_force),9*Nreals2)
            Call Get_force_values(Work(Iequi_force),Work(Itemp),Ndims,
     +                            Iuniteq)
         Endif 
         Do Imode = 1, Nmodes
            Do J = 1, Nfdif_points
               Call Dzero(Work(Iplus_off),9*Nreals2)
               Call Dzero(Work(Imins_off),9*Nreals2)

               If (PREORI_Present .OR. MREORI_Present) Then
                  Call Get_reorient(RP,Iunitp3)
                  Call Get_reorient(RM,Iunitm3)
               Endif 

#ifdef _DEBUG_LVLM
               Write(6,"(a)") "Getting force constant from plus side"
#endif 
               Call Get_force_values(Work(Iplus_off),Work(Itemp),
     +                               Ndims,Iunitp1)
               Call Reorient_fc(Work(Iplus_off),RP,Work(Iend),Maxcor,
     +                          Ndims,Nreals)
#ifdef _DEBUG_LVLM
               Write(6,"(a)") "Getting force constant from minus side"
#endif 
               Call Get_force_values(Work(Imins_off),Work(Itemp),
     +                               Ndims,Iunitm1)
               Call Reorient_fc(Work(Imins_off),RM,Work(Iend),Maxcor,
     +                          Ndims,Nreals)
            Enddo

            Step_size = Deltaq(Imode)
            Itype     = 2

            Call Form_2dprop_deriv(Work(Iplus),Work(Imins),
     +                             Work(I1dif_force),Ndims,Nmodes,
     +                             Step_size,1,Nfdif_points,Imode,
     +                             Itype)
            
         Enddo 
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(2a)") "The first derivative of the force constant",
     +                " matrix"
      Ioff = 0
      Do Imode = 1, Nmodes
         Call output(Work(I1dif_force+Ioff),1,Ndims,1,Ndims,Ndims,
     +               Ndims,1)
         Ioff = Ioff + 9*Nreals2 
      Enddo
      Write(6,*)
#endif
      Endif
C
C For debugging purposed, obtain the Hessian numrically. If this 
C matches with the analytica Hessian, then we know I am doing the
C finite difference correctly.

      Call Gfname("pGRAD",Fname,Length)
      Inquire(File=Fname(1:Length),Exist=PGRAD_Present)
      If (PGRAD_Present) Open(Unit=Iunitp1,File=Fname(1:Length),
     +                         Form="Formatted")
      Call Gfname("mGRAD",Fname,Length)
      Inquire(File=Fname(1:Length),Exist=MGRAD_Present)
      If (MGRAD_Present) Open(Unit=Iunitm1,File=Fname(1:Length),
     +                         Form="Formatted")
      Call Gfname("eGRAD",Fname,Lengthj)
      Inquire(File=Fname(1:Length),Exist=EGRAD_Present)
      If (EGRAD_Present) Open(Unit=Iuniteq,File=Fname(1:Length),
     +                         Form="Formatted")
      If (PGRAD_present .or. MGRAD_present .or. EGRAD_present) 
     +    GRAD = .True.

      If (GRAD) Then

         Nfdif_points = 1
         Iequi_grad  = Iend
         I1dif_grad  = Iequi_grad + 3*Nreals
         Iend        = I1dif_grad + 3*Nreals*Nmodes

         Iplus  = Iend 
         Imins  = Iplus + 3*Nreals
         Iend   = Imins + 3*Nreals
         Maxcor = Maxcor - Iend 
   
         Ndims        = Nreals
         Iplus_off    = Iplus
         Imins_off    = Imins

         If (EGRAD_Present) Then
            Call Dzero(Work(Iequi_grad),3*Nreals)
#ifdef _DEBUG_LVLM
               Write(6,"(a)") "Getting gradient from equlibrium"
#endif 
            Call Get_grad_values(Work(Iequi_grad),Ndims,Iuniteq)
         Endif 

         Do Imode = 1, Nmodes
            Do J = 1, Nfdif_points
               Call Dzero(Work(Iplus_off),3*Nreals)
               Call Dzero(Work(Imins_off),3*Nreals)

CSSS               If (PREORI_Present .OR. MREORI_Present) Then
CSSS                  Rewind(Iunitp3)
CSSS                  Rewind(Iunitm3)
CSSS                  Call Get_reorient(RP,Iunitp3)
CSSS                  Call Get_reorient(RM,Iunitm3)
CSSS               Endif 
#ifdef _DEBUG_LVLM
               Write(6,"(a)") "Getting gradient from plus side"
#endif 
               Call Get_grad_values(Work(Iplus_off),Ndims,Iunitp1)
CSSS               Call Reorient_grd(Work(Iplus_off),RP,Work(Iend),Maxcor,
CSSS     +                           Nreals)
#ifdef _DEBUG_LVLM
               Write(6,"(a)") "Getting gradient from minus side"
#endif 
               Call Get_grad_values(Work(Imins_off),Ndims,Iunitm1)
CSSS               Call Reorient_grd(Work(Imins_off),RP,Work(Iend),Maxcor,
CSSS     +                           Nreals)
            Enddo

            Step_size = Deltaq(Imode)
            Itype     = 2

c Obtain the first derivatives (Hessian)

            Call Form_1dprop_deriv(Work(Iplus),Work(Imins),
     +                             Work(I1dif_grad),Nreals,Nmodes,
     +                             Step_size,1,Nfdif_points,Imode,
     +                             Itype)

         Enddo
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(2a)") "Numerically obtained Hessian matrix"
      Call output(Work(I1dif_grad),1,3*Ndims,1,Nmodes,3*Ndims,Nmode,1)
      Write(6,*)
#endif
      Endif 

      Close(Iunitp1)
      Close(Iunitm1)
      Close(Iuniteq)
      Close(Iunitp3)
      Close(Iunitm3)

C Lets make the assumption that the files set from the equilibrium
C virational frequency calculation is saved. Read the equilibrium
C force data. The normal modes returned from the following call
C are not mass weighted. We reconstruct them in get_nomdes.F. So
C what is being read is irrelevant. 

      Trans_state   = .False.
      Mass_weigh_nm = .False.
      Mass_weigh_gr = .False.
      Get_Hess      = .True.
      Get_Grad      = .False.

      IHess_t = Iend 
      IHess_f = IHess_t + 9*Natoms*Natoms 
      I010    = Ihess_f + 9*Natoms*Natoms
      I020    = I010    + 3*Natoms 
      I030    = I020    + 9*Nreals*Nreals 
      I040    = I030    + 3*Nreals
      I050    = I040    + 3*Nreals
      I060    = I050    + 3*Nreals
      I070    = I060    + Natoms 
      I080    = I070    + 3*Nreals 
      Iend    = I080    + 3*Nreals 
   
      If (Iend .Gt. Maxcor) Call Insmem("@-vibavg_driver",Iend,
     +                                   Maxcor)

      Call Get_refvib_data(Work(I010),Work(I020),Work(I030),
     +                     Work(I040),Atmmass,Atmlabel,Symlabel,
     +                     Work(I050),Work(I060),Work(I070),Work(I080),
     +                     Work(Ihess_t),Work(Ihess_f),Imap,Vib_type,
     +                     Nreals,Natoms,Nmodes,B2ang,Au2invcm,
     +                     Trans_state,Mass_Weigh_nm,Mass_Weigh_gr,
     +                     Get_hess,Get_Grad,Ivib_level,Igrad_calc,
     +                     Icol)
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "The vibrational frequencies and symmetries"
      Write(6, "(5F12.5,1x)") (Work(I040+I-1), I = 1, Nmodes)
      Write(6,*)
      Write(6, "(6(a4))") (SymLabel(I), I = 1, Nmodes)
      Write(6,*) "The internal<->external ordering matrix"
      Write(6,"(10(1x,i2))") (Imap(i), i=1,Natoms)
      
#endif 
#ifdef _REORDER     
      Call Reord_2zmat(Atmmass,Atmlabel,Iatmchrg,Imap,Atmmass_zmtord,
     +                 Atmlabel_zmtord,Iatmchrg_zmtord,Natoms,Nreals)
#endif 

      Call Izero(Isotopes,Nreals)
      Call Gfname("iSOTOPES",Fname,Length)
      Inquire(File=Fname(1:Length),Exist=ISOTOPES_Present)
      If (ISOTOPES_present) Then
         Open(Unit=Iunitis,File=Fname(1:Length),Form="Formatted")
         If (Natoms .Ne. Nreals) Then
            Ncount = 0
            Do Iatmvml = 1, Natoms
               Iatmzmt = Imap(Iatmvml)
               If (Iatmzmt .Ne. Inull) Then
                  Nocunt = Ncount + 1
                  Atmmass_zmtord(Ncount)  = Atmmass(Iatmzmt) 
                  Iatmchrg_zmtord(Ncount) = Iatmchrg(Iatmzmt) 
               Endif
            Enddo 
            Do Iatm = 1, Ntoms 
               Atmmass(Iatm)  = Atmmass_zmtord(Iatm)
               Iatmchrg(Iatm) = Iatmchrg_zmtord(Iatm)
            Enddo 
            Call Fixdummy(Atmlabel_vmlord,Tmp1,Tmp2,Tmp3,Imap,
     +                    Iatmchrg,Natoms,Nreals)
            Do Iatmvml = 1, Nreals 
               Iatmzmt = Imap(Iatmvml)
               Atmmass_zmtord(Iatmzmt)  = Atmmass(Iatmvml)
               Iatmchrg_zmtord(Iatmzmt) = Iatmchrg(Iatmvml) 
               Atmlabel_zmtord(Iatmzmt) = Atmlabel_vmlord(Iatmvml)
            Enddo 
         Else
            Do Iatm = 1, Nreals
               Atmmass_zmtord(Iatm)  = Atmmass(Iatm)
               Iatmchrg_zmtord(Iatm) = Iatmchrg(Iatm)
               Atmlabel_zmtord(Iatm) = Atmlabel(Iatm)
            Enddo 
         Endif 
         Call Get_isodata(Atmmass_zmtord,Atmlabel_zmtord,
     +                    Iatmchrg_zmtord,Isotopes,Atmlabel,
     +                    Iunitis,Nreals)
         Do Iatm = 1, Nreals
            If (Isotopes(Iatm) .EQ. Ione) Then
               Iatmchrg(Iatm) = Iatmchrg_zmtord(Iatm)
            Endif
            Atmmass(Iatm)  = Atmmass_zmtord(Iatm)
         Enddo 
         Close(Iunitis)
      Endif 
C
C The ZMAT ordered Hessian in Cartesian cooordinates are returned 
C in Work(Ihess_t) and Work(Ihess_f). The Hessian in Work(Ihess_t) 
C has no dummy atoms. Obtain the normal modes as eigenvectors of
C the mass-weighted Hessian.

      Inmrlm = Ihess_f
      Ivfreq = Inmrlm + 9*Nreals*Nreals 
      Ivproj = Ivfreq + 3*Nreals 
      Icfc   = Ivproj + 9*Nreals*Nreals 
      Iend   = Icfc   + Nmodes*Nmodes*Nmodes 
      Maxcor = Maxcor - Iend 

#ifdef _DEBUG_LVL0
      write(6,"(a)") " The reference geometry"   
      Write(6,"(3(1x,F12.7))") (work(i050-1+j),j=1,3*Nreals)
#endif
      Call Hproject_main(Work(Ivproj),Work(I050),Atmmass,Nreals)

      Call Get_nmodes(Work(Ihess_t),Work(Inmrlm),Work(Ivfreq),
     +                Work(Iend),Maxcor,Atmmass,Nreals)
C
C This is for debugging. Compute the Hessian from numerically obtained 
C gradients (computed in the above block specified by GRAD).
    
      If (Grad) Then
         Call Get_fc(Work(I1dif_grad),Work(Inmrlm),Work(Ivproj),
     +               Work(Ivfreq),Work(Iend),Maxcor,Atmmass,Nreals,
     +               Nmodes,Iamlinear)
      Endif 

C The cubic force field stored as (Cart,Cart,Normal). Convert this
C to (Normal,Normal,Normal) using the eigenvectors of the equilibrium
C Hessian.

      If (Force) Then
          Call Get_cfc(Work(I1dif_force),Work(Icfc),Work(Ivfreq),
     +                 Work(Inmrlm),Work(Iend),Maxcor,Atmmass,
     +                 Nreals,Nmodes,Iamlinear)
      Endif 

C This will give the very basic vibrationally/thermally averaged 
C properties. 

      Maxcor = Maxcor - Iend 
      If (NMR_spincc) Then
          Ndim = Nreals*Nreals
          Call Vibavg_lvl0(Work(I1dif_spncc),Work(I2dif_spncc),
     +                     Work(Icfc),Work(Ivfreq),
     +                     Work(Iend),Imap,Temp,Maxcor,Ndata,
     +                     Nmodes,Ndim,Nreals,B2ang,Au2Invcm,
     +                     NMR_spincc,.False.,Iamlinear)

C Starting from Work(Ispncc), we have the corrections due to vibrational
C We need to add them to the equlibrium values and print them out.

          Call Finalize(Work(Iend),Iatmchrg,Isotopes,Atmlabel,Temp,
     +                  Maxcor,Nreals,Natoms,Ndata,NMR_spincc,.False.) 
      Endif 

      If (NMR_shifts) Then
          Ndim = Nreals
          Call Vibavg_lvl0(Work(I1dif_shift),Work(I2dif_shift),
     +                     Work(Icfc),Work(Ivfreq),
     +                     Work(Iend),Imap,Temp,Maxcor,Ndata,
     +                     Nmodes,Ndim,Nreals,B2ang,Au2Invcm,
     +                     .False.,NMR_shifts,Iamlinear)

C Starting from Work(Ishifts), we have the corrections due to vibrational
C We need to add them to the equlibrium values and print them out.

          Call Finalize(Work(Iend),Iatmchrg,Isotopes,Atmlabel,Temp,
     +                  Maxcor,Nreals,Natoms,Ndata,.False.,NMR_shifts) 
      Endif 
    
      Return
      End 



