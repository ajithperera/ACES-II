#include "flags.h"

      Program  pccd_drpmo
      Implicit none
      Integer iuhf
      Integer IGrad_calc
      Integer IHess_calc
      Integer Pccd_refr,Itest
      Integer OOmicroItr,OOmacroItr,OOtotalItr
      Integer Length, Junk
      Double Precision Scale 
      Character*80 Fname
      Character*6 Status
      Logical OOexist 
      DOUBLE PRECISION GLOBNORM
      INTEGER MACRONUM,MICRONUM,TOTALITR
      LOGICAL MACCONVG,MICCONVG
      Logical Opt_orbs, OO_constr, OV_constr 
      Logical Pccd,Ccd,Lccd,Symmetry 

c COMMON BLOCKS
#include "syminf.com"
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "sympop.com"
#include "machsp.com"

      Common /OO_info/Opt_orbs,OO_constr,OV_constr
      Common /CALC/PCCD,CCD,LCCD
      Common /Symm/Symmetry

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
      Pccd  = .FALSE.
      Ccd   = .FALSE.
      Lccd  = .FALSE.
      CALL Getrec(0,"JOBARC","PCCD_RUN",Length,Junk)
      If (Length .Gt. 0) pCCD = .True.
      If (IFLAGS(h_IFLAGS_calclevel).EQ.53) Lccd = .True.
      If (IFLAGS(h_IFLAGS_calclevel).EQ.54) Ccd  = .True. 

      IF(IFLAGS(h_IFLAGS_symmetry).EQ.2) Symmetry  = .TRUE.

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

      If (Iflags2(h_IFLAGS2_grad_calc) .Eq. 2) Then
          IGrad_calc = 0
          IHess_calc = 0
          Scale     = Iflags(h_IFLAGS_fd_stepsize)
      Endif
      Opt_orbs  = (Iflags2(h_IFLAGS2_oo_method) .Gt. 0)
      Oo_constr = (Iflags2(h_IFLAGS2_oo_constrain) .Eq. 1)
      Ov_constr = (Iflags2(h_IFLAGS2_ov_constrain) .Eq. 1)
C
      If (Igrad_calc .eq. 1 .and. IHess_calc .eq.1) Then

C If Bruekner optimization is sought iterate until Brueckner convergence
C is achieved before computing orbital gradients. 

         Call Pccd_drpmo_do_orbrots(Icore(I0),Icrsiz,Iuhf,IGrad_calc,
     +                              IHess_calc,Scale,OOmicroItr,
     +                              OOmacroItr,OOtotalItr)
      Endif

 10   Continue

      Call aces_fin
C
c ----------------------------------------------------------------------
      Stop
      End

