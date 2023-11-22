#include "flags.h"

      Program  pccd 
      Implicit none
      Integer iuhf
      Integer IGrad_calc
      Integer IHess_calc
      Integer Pccd_refr,Itest
      Integer OOmicroItr,OOmacroItr,OOtotalItr
      Integer Length,Idrop
      Double Precision Scale 
      Character*80 Fname
      Character*6 Status
      Logical OOexist 
      DOUBLE PRECISION GLOBNORM
      INTEGER MACRONUM,MICRONUM,TOTALITR
      LOGICAL MACCONVG,MICCONVG
      Logical Opt_orbs, OO_constr, OV_constr, All_mos

c COMMON BLOCKS
#include "syminf.com"
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "sympop.com"
#include "machsp.com"

      Common /OO_info/Opt_orbs,OO_constr,OV_constr
c ----------------------------------------------------------------------
      Call Aces_init(icore, i0,icrsiz, iuhf, .true.)
#ifdef _NOSKIP
      If (Iuhf .Eq. 1) Then
         Write(6,*)
         Write(6,"(2a)") " UHF reference pCCD calculations are not",
     +                   " supported changed the refrence to RHF"
         Write(6,"(a)")  " and rerun to proceed."
         Call Errex()
      Endif 
#endif 
      CALL PUTREC(20,"JOBARC","PCCD_RUN",1,1)

      IGrad_calc = 1
      IHess_calc = 1
      Status     = "Begin "
      
      Call Gfname("ooinfo.dat",Fname,Length)
      Inquire(file=Fname(1:Length),exist=OOexist)
      If (OOexist) Then
         Open(Unit=1011,File=Fname(1:Length),Status="old")
         Read(1011,*) MACRONUM,MICRONUM,TOTALITR,MACCONVG,
     +                          MICCONVG,GLOBNORM
         call putrec(-20,"JOBARC","MACRONUM",1,MACRONUM)
         call putrec(-20,"JOBARC","MICRONUM",1,MICRONUM)
         call putrec(-20,"JOBARC","TOTALITR",1,TOTALITR)
         CALL putrec(-20,"JOBARC","MACCONVG",1,.False.)
         CALL putrec(-20,"JOBARC","MICCONVG",1,.False.)
         call putrec(-20,"JOBARC","GLOBNORM",1,0.0d0)
         close(1011, status='delete')
      Endif
      call getrec(-20,"JOBARC","MICRONUM",1,OOmicroItr)
      call getrec(-20,"JOBARC","MACRONUM",1,OOmacroItr)
      call getrec(-20,"JOBARC","TOTALITR",1,OOtotalItr)
      Call Getrec(20,"JOBARC",'NDROPGEO',1,Idrop)

      If (Iflags2(h_IFLAGS2_grad_calc) .Eq. 2) Then
          IGrad_calc = 0
          IHess_calc = 0
          Scale     = Iflags(h_IFLAGS_fd_stepsize)
      Endif

      All_mos   = (Idrop .Eq. 0)
      Opt_orbs  = (Iflags2(h_IFLAGS2_oo_method) .Gt. 0)
      Oo_constr = (Iflags2(h_IFLAGS2_oo_constrain) .Eq. 1)
      Ov_constr = (Iflags2(h_IFLAGS2_ov_constrain) .Eq. 1)
C
      If (Igrad_calc .eq. 1 .and. IHess_calc .eq.1) Then

C If Bruekner optimization is sought iterate until Brueckner convergence
C is achieved before computing orbital gradients. 

         If (Iuhf .Eq. 0 .AND. All_mos) Then 
            Call Pccd_banner(6,Status)
            Call Pccd_rdriver(Icore(i0),icrsiz/iintfp,Iuhf)

            Pccd_refr = Iflags2(h_IFLAGS2_pccd_reference)
            Call Getrec(-20, 'JOBARC', 'BRUKTEST', 1, Itest)

            If ((Pccd_refr .Eq. 2) .And. 
     +           Itest .Ne. 1) Go to 10

            Call Pccd_ldriver(Icore(i0),icrsiz/iintfp,Iuhf)

            If (Opt_orbs) Then
               Call Pccd_init_2pdens_lists(Icore(i0),icrsiz/iintfp,Iuhf)
               Call Pccd_init_orbhes_lists(Icore(i0),icrsiz/iintfp,Iuhf)

               Call Pccd_dens(Icore(i0),icrsiz/iintfp,Iuhf)
               Call Pccd_do_orbrots(Icore(I0),icrsiz/iintfp,Iuhf,
     +                              Igrad_calc,IHess_calc,Scale,
     +                              OOmicroItr,OOmacroItr,OOtotalItr)
            Endif 
         Else
            Call Pccd_banner(6,Status)
            Call Pccd_urdriver(Icore(i0),icrsiz/iintfp,Iuhf)

            If ((Pccd_refr .Eq. 2) .And. 
     +           Itest .Ne. 1) Go to 10
            Pccd_refr = Iflags2(h_IFLAGS2_pccd_reference)
            Call Getrec(-20, 'JOBARC', 'BRUKTEST', 1, Itest)

            Call Pccd_uldriver(Icore(i0),icrsiz/iintfp,Iuhf)

            If (Opt_orbs) Then
               Call Pccd_do_orbrots(Icore(I0),icrsiz/iintfp,Iuhf,
     +                              Igrad_calc,IHess_calc,Scale,
     +                              OOmicroItr,OOmacroItr,OOtotalItr)
            Endif 
         Endif 
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

