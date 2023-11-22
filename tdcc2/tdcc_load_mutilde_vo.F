      Subroutine Tdcc_load_mutilde_vo(Work,Memleft,Memleft_modf,Iuhf,
     +                                Ioff_vo,Irrepx)

      Implicit Integer (A-Z)

      Dimension Work(Memleft), Ov_start(2)
      Dimension Ioff_vo(8,2)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "symloc.com"
#include "sym.com"

C The Mu_vo vector reside at the end of the Work array ordered as
C UHF Mu_vo(AA)Mu_vo(BB)
C RHF Mu_vo(AA)
C The vo_off array give offset to the right-hand irrep. 

#ifdef _DEBUG_LVLM
        Write(6,*)
        Write(6,"(a)") "---Entered  tdcc_load_mu_vo---"
        Write(6,*)
#endif
      List_mu0 = 392
      Memleft_modf = Memleft 
      Do Ispin = 2, 2-Iuhf, -1

         If (Iuhf .EQ. 0) Then
            Ndim_vo = Irpdpd(Irrepx,9) 
         Else
           If (Ispin .EQ. 2) Then
            Ndim_vo = Irpdpd(Irrepx,10) 
           Else
            Ndim_vo = Irpdpd(Irrepx,9) 
           Endif 
         Endif

         Ov_start(Ispin) = Memleft_modf-Ndim_vo+1
         If (Iuhf .EQ. 0) Ov_start(1) = Ov_start(2)
         Memleft_modf = Memleft_modf - Ndim_vo
         Ioff = Ov_start(Ispin)

         Do Irrep_i = 1, Nirrep
            Irrep_a = Dirprd(Irrep_i,Irrepx)
            Ioff_vo(Irrep_i,Ispin) = Ioff 
            If (Iuhf .EQ. 0) Ioff_vo(Irrep_i,1) = Ioff
            Ndim = Vrt(Irrep_a,Ispin) * Pop(Irrep_i,Ispin)
            Ioff = Ioff + Ndim 
         Enddo 

         If (Iuhf .EQ. 0) then
            Call Getlst(Work(Ov_start(Ispin)),1,1,1,1,List_mu0)
         Else 
            Call Getlst(Work(Ov_start(Ispin)),1,1,1,ISpin,List_mu0)
         Endif 

#ifdef _DEBUG_LVLM
      Write(6,*)
      call Checksum("Tdcc_load_mu_vo,Mu(a,i):",Work(Ov_start(Ispin)),
     +          Ndim_vo,S)
#endif 
      Enddo 

#ifdef _DEBUG_LVLM
       write(6,"(a)") "The vo offset array"
       write(6,"(6(1x,i8))") ((ioff_vo(i,j),i=1,Nirrep),j=1,2)
#endif 
      Return
      End
