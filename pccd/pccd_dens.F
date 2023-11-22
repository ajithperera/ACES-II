      Subroutine Pccd_dens(Work,Maxcor,Iuhf)

      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Maxcor)
      LOGICAL PCCD,CCD,LCCD

#include "machsp.com"
#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      COMMON/CALC/PCCD,CCD,LCCD
      COMMON/NHFREF/NONHF_REF
 
      Data Ione,Onem /1,-1.0D0/

      I0 = Ione
#ifdef _DEBUG_LVL0
       call pccd_checkintms(work(I0),maxcor,iuhf,0)
#endif
C G(ij,kl),G(ab,cd),G(ai,bj) is constructed in gamdrv. The other
C Gamma elements are zero for pCCD.

      Call Pccd_gamdrv(Work(I0),Maxcor,Iuhf)

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
      Call Dzero(Work(IDOV),Nt(1)+IUHF*Nt(2))
      CALL Pccd_densoo(Work(IDOO),Work(I0),Maxcor,Iuhf)
      CALL Pccd_densvv(Work(IDVV),Work(I0),Maxcor,Iuhf)
      If (Nonhf_ref) Call Pccd_densvo(Work(IDOV),Work(I0),Maxcor,
     +                                Iuhf)

C G(ij,ab) + G(ab,ij)
#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(a)") "Checksums of G(ij,ab)"
#endif
      Call Pccd_gamma1(Work(I0),Maxcor,Iuhf)
CSSS      If (Lccd) Call Pccd_setg(Work(I0),Maxcor)

#ifdef _DEBUG_LVL0
       Write(6,*)
       Write(6,"(a)") "The checksums of the Density blocks"
CSSS       Call Dscal(Nfmi(1)+IUHF*Nfmi(2),Onem, Work(IDOO),1)
CSSS       Call Dscal(Nfea(1)+IUHF*Nfea(2),Onem, Work(IDVV),1)
       Call checksum("DENSOO   :", Work(IDOO),Nfmi(1)+IUHF*Nfmi(2))
       Call checksum("DENSVV   :", Work(IDVV),Nfea(1)+IUHF*Nfea(2))
       Call checksum("DENSVO   :", Work(IDOV),Nt(1)+IUHF*Nt(2))
#endif

      Return
      End 

