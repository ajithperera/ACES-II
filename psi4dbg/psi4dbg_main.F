#include "flags.h"

      Program  psi4dbg
      Implicit none
      Integer iuhf
      Integer IGrad_calc
      Integer IHess_calc
      Integer Psi4dbg_refr,Itest
      Integer OOmicroItr,OOmacroItr,OOtotalItr
      Integer Length
      Integer Junk,Ncycle
      Double Precision Scale 
      Character*80 Fname
      Logical OOexist 
      Logical Macro_iter
      DOUBLE PRECISION GLOBNORM
      INTEGER MACRONUM,MICRONUM,TOTALITR
      LOGICAL MACCONVG,MICCONVG



c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "sympop.com"
#include "machsp.com"
       
c ----------------------------------------------------------------------
      Call Aces_init(icore, i0,icrsiz, iuhf, .true.)
      If (Iuhf .Eq. 1) Then
         Write(6,*)
         Write(6,"(2a)") " UHF reference pCCD calculations are not",
     +                   " supported changed the refrence to RHF"
         Write(6,"(a)")  " and rerun to proceed."
         Call Errex()
      Endif 
      CALL PUTREC(20,"JOBARC","PCCD_RUN",1,1)

      IGrad_calc = 1
      IHess_calc = 1
      
      Call Gfname("ooinfo.dat",Fname,Length)
      Inquire(file=Fname(1:Length),exist=OOexist)
      If (OOexist) Then
         Open(Unit=1011,File=Fname(1:Length),Status="old")
         Read(1011,*) MACRONUM,MICRONUM,TOTALITR,MACCONVG,
     +                         MICCONVG,GLOBNORM
         call putrec(20,"JOBARC","MACRONUM",1,MACRONUM)
         call putrec(20,"JOBARC","MICRONUM",1,MICRONUM)
         call putrec(20,"JOBARC","TOTALITR",1,TOTALITR)
         CALL putrec(20,"JOBARC","MACCONVG",1,.False.)
         CALL putrec(20,"JOBARC","MICCONVG",1,.False.)
         call putrec(20,"JOBARC","GLOBNORM",1,0.0d0)
         close(1011, status='delete')
      Endif
      call getrec(-20,"JOBARC","MICRONUM",1,OOmicroItr)
      call getrec(-20,"JOBARC","MACRONUM",1,OOmacroItr)
      call getrec(-20,"JOBARC","TOTALITR",1,OOtotalItr)

      If (Iflags2(h_IFLAGS2_grad_calc) .Eq. 2) Then
          IGrad_calc = 0
          IHess_calc = 0
          Scale     = Iflags(h_IFLAGS_fd_stepsize)
      Endif
C
      If (Igrad_calc .eq. 1 .and. IHess_calc .eq.1) Then

C If Bruekner optimization is sought iterate until Brueckner convergence
C is achieved before computing orbital gradients. 

        Call getrec(0,"JOBARC","ORBOPITR",Ncycle,junk)
        If (Ncycle .lt. 0) then
           Macro_iter = .True.
        Else
            call getrec(20,"JOBARC","ORBOPITR",1,ncycle)
            Macro_iter = .True. 
        Endif

         Call Psi4dbg_rdriver(Icore(i0),icrsiz/iintfp,Iuhf,
     +                        .Not.Macro_iter)
         Call Psi4dbg_ldriver(Icore(i0),icrsiz/iintfp,Iuhf,
     +                        .Not.Macro_iter)
#ifdef _DEBUG_LVL0
         call psi4dbg_check_t2(icore(i0),icrsiz,iuhf)
#endif   

         Psi4dbg_refr = Iflags2(h_IFLAGS2_pccd_reference)
         Call Getrec(-20, 'JOBARC', 'BRUKTEST', 1, Itest)

         If (Psi4dbg_refr .Eq. 2 .And. Itest .Ne. 1) Go to 10

         Call Psi4dbg_init_2pdens_lists(Icore(i0),icrsiz/iintfp,Iuhf)

         Call Psi4dbg_dens(Icore(i0),icrsiz/iintfp,Iuhf)
         Call Psi4dbg_do_orbrots(Icore(I0),icrsiz/iintfp,Iuhf,
     +                           Igrad_calc,IHess_calc,Scale,
     +                           OOmicroItr,OOmacroItr,
     +                           OOtotalItr)
      Else
         Call Pccd_form_nmrcl_grad_hess(Icore(I0),icrsiz/iintfp,Iuhf,
     +                                  Igrad_calc,IHess_calc,Scale)
      Endif

 10   Continue

      Call aces_fin
C
c ----------------------------------------------------------------------
      Stop
      End

