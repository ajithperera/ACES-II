#include "flags.h"
      Subroutine Psi4dbg_dens(Work,Maxcor,Iuhf)

      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Maxcor)
      Logical Non_hf

#include "machsp.com"
#include "sym.com"
#include "sympop.com"
#include "syminf.com"
#include "flags.com"

      Data Ione /1/

      I0 = Ione
      Non_hf = (Iflags(h_IFLAGS_nonhf) .gt.0)
      Call Psi4dbg_init2pdm(Work(I0),Maxcor)
      Call Psi4dbg_gamdrv(work(I0),Maxcor,Iuhf)

      IDOO   = I0 + MAXCOR-(Nfmi(1)+IUHF*Nfmi(2))
      MAXCOR = MAXCOR - (Nfmi(1)+IUHF*Nfmi(2))
C
      IDVV   = IDOO-(Nfea(1)+IUHF*Nfea(2))
      MAXCOR = MAXCOR-(Nfea(1)+IUHF*Nfea(2))
C
      IDOV   = IDVV-(Nt(1)+IUHF*Nt(2))
      MAXCOR = MAXCOR-(Nt(1)+IUHF*NT(2))
C
      Call Dzero(Work(IDOO),Nfmi(1)+IUHF*Nfmi(2))
      Call Dzero(Work(IDVV),Nfea(1)+IUHF*Nfea(2))
      CALL Psi4dbg_densoo(Work(IDOO),Work(I0),Maxcor,Iuhf)
      CALL Psi4dbg_densvv(Work(IDVV),Work(I0),Maxcor,Iuhf)
      If (Non_hf) CALL Psi4dbg_densvo(Work(IDOV),Work(I0),
     +                                Maxcor,Iuhf)

C G(ij,ab) + G(ab,ij)
#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(a)") "Checksums of G(ij,ab)"
#endif
       Call Pccd_Gamma1(Work(I0),Maxcor,Iuhf)
      
       Call Psi4dbg_setg(Work(I0),Maxcor)

#ifdef _DEBUG_LVL0
       Write(6,*)
       Write(6,*) "The checksums of the Density blocks"
       Call checksum("DENSOO ", Work(IDOO),Nfmi(1)+IUHF*Nfmi(2))
       Call checksum("DENSVV ", Work(IDVV),Nfea(1)+IUHF*Nfea(2))
       If (non_hf) Call checksum("DENSOV ", Work(IDOV),
     +                            Nt(1)+IUHF*Nt(2))
#endif

      Return
      End 

