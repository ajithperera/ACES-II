#include "flags.h"
      Subroutine Pccd_do_orbrots(Work,Maxcor,Iuhf,IGrad_calc,IHess_calc,
     +                           Scale,OOmicroItr,OOmacroItr,OOtotalItr)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxcor)
      Logical Grdcnv
      Logical pCCD,CCD,LCCD
      Logical OneP_ONLY

      Common /CALC/pCCD,CCD,LCCD
      Common /ORBR_HESS/ONEP_ONLY

#include "info.com"
#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "flags2.com"

      Data Ione,Itwo,Izero,One,Half /1,2,0,1.0D0,0.50D0/

      Call Getrec(20,"JOBARC",'NDROPGEO',1,Idrop)
       
      If (Iuhf .EQ. 0 .AND. Idrop .EQ. 0) Then
         If(.Not. ONEP_only) Call Pccd_form_htau_pqrs(Work,Maxcor,Iuhf)
         Call Pccd_form_htau_pq(Work,Maxcor,Iuhf,Igrad_calc,
     +                          IHess_calc,Scale)
      Else

C Uhf only or UHF and RHF for forzen core orbitals.

         Call Pccd_form_uhtau(Work,Maxcor,Iuhf)
      Endif 

      Call GetSDInfo(nocca,noccb,nvrta,nvrtb,nbas,nvecDim,Naobfns)

      If (Iflags2(h_IFLAGS2_oo_method) .EQ. 1) then ! BB Steep. Desc.
         call Drive_BB_SD(Work,Maxcor,Iuhf,OOmicroItr,OOmacroItr,
     +                    nocca,nvrta,nvecDim)
      else if (Iflags2(h_IFLAGS2_oo_method) .EQ. 2) then !BFGS
        print*,'Have not include L-BFGS yet; to do'
      else if (Iflags2(h_IFLAGS2_oo_method) .EQ. 3) then !FullNR
        print*,'Have not included FullNR yet; to do'
      else if (Iflags2(h_IFLAGS2_oo_method) .EQ. 4) then !AMSgrad
          call Drive_AMSgrad(Work,Maxcor,Iuhf,OOmicroItr,OOmacroItr,
     +                        nocca,nvrta,nvecDim)
      else if (Iflags2(h_IFLAGS2_oo_method) .EQ. 5) then !AMSgrad
          Call Getrec(-20,'JOBARC','ORBOPITR',Ione,Ncycle)
          Call Pccd_rotgrd(Work,Maxcor,Nbas,Nocca,Noccb,Nvrta,Nvrtb,
     +                     Iuhf,Grdcnv,Tol,Grd_max,Grd_rms,Ncycle)
      endif
C 
      If (Grdcnv .And. .Not. Ncycle .Eq. 0) Call Pccd_finalize(Work,
     +                                      Maxcor,Iuhf,Tol,Grd_max,
     +                                      Grd_rms)
      If (IFLAGS(h_IFLAGS_reference).GT.0) 
     +IFLAGS(h_IFLAGS_reference) = Ione
      IFLAGS(h_IFLAGS_scf_maxcyc)= Izero
      IFLAGS(h_IFLAGS_nonhf)     = Ione
C Do not turn on non-HF key-word sinply because this indirectly add singles 
      IFLAGS(h_IFLAGS_qrhf_orb)  =Izero
      IFLAGS(h_IFLAGS_guess)     =Izero

      Call Putrec(20,'JOBARC','IFLAGS  ',100,IFLAGS)

      Return
      End 
 
