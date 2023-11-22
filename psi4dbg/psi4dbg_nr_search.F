#include "flags.h"
      Subroutine Psi4dbg_nr_search(Grd,Grd_oo,Grd_vv,Grd_vo,
     +                             Grd_ov,Grd_stat,Lenoo,Lenvv,Lenvo,
     +                             Nbas,Nocc,Nvrt,Work,Maxcor)

      Implicit Double Precision(A-H,O-Z)

      
      Dimension Grd_stat(6)
      Dimension Grd(Nbas,Nbas)
      Dimension Grd_oo(Lenoo)
      Dimension Grd_vv(Lenvv)
      Dimension Grd_vo(Lenvo)
      Dimension Grd_ov(Lenvo)
      Dimension Work(Maxcor)

#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "flags2.com"

      Data Ione,Izero,One,Half,Dnull,Two/1,0,1.0D0,0.5D0,0.0D0,
     +                                   2.0D0/
      Data Onem/-1.0D0/

      I000 = Ione
      I010 = I000 + Nbas
      I020 = I010 + Lenvo
      I030 = I020 + Lenvo
      Iend = I030 + Nbas*Nbas
      Memleft = Maxcor - Iend 

      Call Getrec(20,"JOBARC","SCFEVALA",Nbas*Iintfp,Work(I000))

      Indi   = Izero
      Inda   = Izero 
      Indv   = Izero 
      Indo   = Izero 
      Indv0  = Izero 
      Indo0  = Izero 
      Ijunk  = Izero
      Irrepx = Ione 

      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)

         Nocci = Pop(Irrepr,1)
         Nvrti = Vrt(Irrepl,1)

         Do indxi = 1, Nocci
            Indi = Indi + 1
            Indv = Indv0
            Do indxa = 1, Nvrti
               Inda = Inda + 1
               Indv = Indv + 1
               Work(I010+Inda-1) = (Work(I000+Indv-1+Nocc) - 
     +                             Work(I000+Indi-1))
#ifdef _DEBUG_LVLM
               write(6,"(3(1x,F15.10))") Work(I000+Indv-1+Nocc),
     +                                   Work(I000+Indi-1),
     +                                   Work(I010+Inda-1) 
#endif 
            Enddo
         Enddo
         Indv0 = Indv0 + Nvrti
      Enddo 

      Nsize = Inda
      Ioff = I010
      Joff = I020 
      Do Irrepr= 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)
         Nocci = Pop(Irrepr,1)
         Nvrti = Vrt(Irrepl,1)
         Call Transp(Work(Ioff),Work(Joff),Nocci,Nvrti)
         Ioff = Ioff + Nocci*Nvrti
         Joff = Joff + Nocci*Nvrti
      Enddo 

      If (Lenvo .Ne. Nsize) Then
         Write(6,"(2a)") " There is an internal inconsistency in OV/VO",
     +                   " length."
         Call Errex
      Endif

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The VO/OC orbital rotation hessians"
      Ioff = I010
      Joff = I020 
      Do Irrepr= 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)
         Nocci = Pop(Irrepr,1)
         Nvrti = Vrt(Irrepl,1)
         Call output(work(ioff),1,nvrti,1,nocci,nvrti,nocci,1)
         Call output(work(joff),1,nocci,1,nvrti,nocci,nvrti,1)
         Ioff = Ioff + Nocci*Nvrti
         Joff = Joff + Nocci*Nvrti
      Enddo 
#endif
#ifdef _DEBUG
C ---- This block was added to use psi4 gradients for debugging----

       open(unit=911,file="psi4_grad",status="old")
       call getrec(0,"JOBARC","ORBOPITR",ncycle,junk)
       if (ncycle .lt. 0) then 
           ncycle = 1
       else
           call getrec(20,"JOBARC","ORBOPITR",1,ncycle)
       endif 
       write(6,*)
       write(6,"(a,i4)") " The optimization cycle:", ncycle 
       Do i=1,ncycle
          read(911,*) grad_vo_psi4
       enddo 
       write(6,"(a,F15.10)") " The orb. rot. grad from psi4:",
     +                       grad_vo_psi4
       Grd_vo(2) = Grad_vo_psi4
       Grd_ov(2) = Grad_vo_psi4
       Write(6,"(a)") "The VO-MO gradient matrices (psi4)"
       write(6,"(6(1x,F12.6))") (Grd_vo(i),i=1,Lenvo)
       Write(6,"(a)") "The OV-MO gradient matrices (psi4)"
       write(6,"(6(1x,F12.6))") (Grd_ov(i),i=1,Lenvo)

C -----------------debug end-------------------------------------
#endif 

      Call Psi4dbg_vstat(Grd,Grd_stat,Nsize)
      Write(6,*)
      Write(6,"(2a)")  " The statistics of the orbital rotation",
     +                 " gradient matrix"
      Write(6,"(1x,2a)") "--------------------------------------------",
     +                   "-------"
      Write(6,"(5x,a,5xa,5xa)") "Minimum grad.", "Maximum grad.", 
     +                           "RMS grad"
      Write(6,*)
      Write(6,"(3(5x,E12.6))") Grd_stat(3), Grd_stat(4), Grd_stat(5)
      Write(6,*)
      Write(6,"(1x,2a)") "--------------------------------------------",
     +                   "-------"
      Write(6,*)

      If (Grd_stat(5) .Lt. Conv_tol) Then
          Write(6,"(2a)") " The micro-iterations of orbital rotation",
     +                    " gradients reached convergence" 
          Call Putrec("20","JOBARC","MICRO_CV",Iintfp,Ione)
      Endif 

C Note that I am using 1/(ea-ei) instead of 1/2(ea-ei) since I work with only
C the alpha (or beta block).  The main loop that rotates orbitals goes over 
C only alpha block only (in the RHF context). See pccd_rotg.F.
      
      Do I = 1, Nsize
         D = -One/Work(I010-1+I)
         Grd_vo(I) = Grd_vo(I)*D
         D = -One/Work(I020-1+I)
         Grd_ov(I) = Grd_ov(I)*D
      Enddo

      Call Dzero(Grd_oo,Lenoo)
      Call Dzero(Grd_vv,Lenvv)
      Call Pccd_frmful(Grd_oo,Grd_vv,Grd,Work(Iend),Memleft,Nbas,0)
      Call Pccd_frmful_ov(Grd_ov,Grd_vo,Grd,Work(Iend),Memleft,
     +                    Nbas,"COPY",0)
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Scalled gradients (G=-H^-1g)" 
      Call output(Grd,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif

C#ifdef _NOSKIP
C Form  G-G(^t). This appears double dipping but help convergence tremendoulsy. 

      Call Transp(Grd,Work(I030),Nbas,Nbas)
      Call Daxpy(Nbas*Nbas,Onem,Work(I030),1,Grd,1)

#ifdef _DEBUG_LVL0
      write(6,*)
      Write(6,"(a)") "Antisymmetrized Scalled gradients (G-G(^t))" 
      Call output(Grd,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif
      Call Dgemm("N","T",Nbas,Nbas,Nbas,Half,Grd,Nbas,Grd,Nbas,
     +            Dnull,Work(I030),Nbas)

      Do I = 1, Nbas
         Grd(I,I) = Grd(I,I) + One
      Enddo

      Call Daxpy(Nbas*Nbas,One,Work(I030),1,Grd,1)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "U=(1+K+1/2K*K)" 
      Call output(Grd,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif 
      Call Pccd_gramschmidt(Grd,Nbas,Nbas)
C#endif 

#ifdef _NOSKIP
C For G - G^(t)

      Call Transp(Grd,Work(I030),Nbas,Nbas)
      Call Daxpy(Nbas*Nbas,Onem,Work(I030),1,Grd,1)
      Do I = 1, Nbas
         Grd(I,I) = Grd(I,I) + One
      Enddo
      Call Pccd_gramschmidt(Grd,Nbas,Nbas)
#endif 

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Unitary check of Kappa" 
      Call Dgemm("N","T",Nbas,Nbas,Nbas,One,Grd,Nbas,Grd,Nbas,
     +           Dnull,Work(I030),Nbas)
      Call output(Work(i030),1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif
   
      Return
      End

