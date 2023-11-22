#include "flags.h"
      Subroutine Drive_lhbar(Work, Length, Iuhf)

      Implicit Double Precision (A-H,O-Z)
      Dimension Work(Length)
      Dimension Ecorr(3)
      Logical Sing 
      Logical MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,CC2
      Integer AAAA_LENGTH_IJKA

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "parcc.com"
#include "files.com"

       Common /METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,CC2
       Common /Energy/Energy(500,2),Ixtrle(500)

       Sing   = .False.
       Icontl = Iflags(h_IFLAGS_cc_conv) 
       Iconvg = 1
       Ncycle = 0
       Sing   = (Iflags(h_IFLAGS_calc) .gt. 9)
       Luout  = 6
       Moints = 50

       Coulomb = .False.
       Call Parread(iuhf)
           write(6,2010) paralpha
           write(6,2011) parbeta
           write(6,2012) pargamma
           write(6,2013) pardelta
           write(6,2014) parepsilon
 2010      format(' PCCSD   alpha parameter : ', F14.6)
 2011      format(' PCCSD    beta parameter : ', F14.6)
 2012      format(' PCCSD   gamma parameter : ', F14.6)
 2013      format(' PCCSD   delta parameter : ', F14.6)
 2014      format(' PCCSD epsilon parameter : ', F14.6)
       If (ispar) Then
          write(6,"(a,a)") ' Perform a parameterized CC HBAR',
     &                    ' calculations'
          Write(6,*)
          if (coulomb) write(6,"(a,a)") " The Coulomb integrals are ",
     $                    "used in W(mbej) intermediate."
          write(6,*)
          Fae_scale    = (Paralpha - 1.0D0)
          Fmi_scale    = (Parbeta  - 1.0D0)
          Wmnij_scale  = Pargamma
          Wmbej_scale  = Pardelta
          Gae_scale    = Paralpha  
          Gmi_scale    = Parbeta  
      Else
          write(6,*) '   Perform a regular CC HBAR calculations'
          write(6,*)
          Fae_scale    = 0.0D0
          Fmi_scale    = 0.0D0
          Wmnij_scale  = 1.0D0
          Wmbej_scale  = 1.0D0
          Gae_scale    = 1.0D0
          Gmi_scale    = 1.0D0
      Endif
C
C Unfortunately I have to use the METH common block. At the moment 
C We will set to CCSD uing calc key-word. 
C
       CCSD = (Iflags(h_IFLAGS_calc) .eq. 10)
C
C Copy the converged T1 and T2 amplitudes to L1 and L2 lists.
C These are our starting L1 and L2 amplitudes.
C
       Call init2_leom(Iuhf, Sing)
       Call Drmove(Work, Length, Iuhf, 100, Sing)
       Call Inilam(Iuhf)
       Call Rnabij(Work,Length,iuhf,"L")
       Call Diislst_leom(1, Iuhf, Sing)

       Irrepx = 1
#ifdef _DEBUG_LVL0
      Write(6,*) "The Hbar lists"
      Call checkhbar(Work,length,Iuhf)
#endif 
       Do While (Iconvg .eq. 1)

          Ncycle = Ncycle + 1
     
          Call Lguess(Work, Length, Iuhf, Sing)

#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals before multiplication"
      Call check_leom(Work,length,Iuhf)
#endif 
          Call Hbar_mult_l(Work, Length, Iuhf)

          Call Newl1_leom(Work,Length,Iuhf)
          Call Newt2(Work,Length,Iuhf)

#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after  multiplication"
      Call check_leom(Work,length,Iuhf)
#endif 
          Call Drtsts(Work, Length, Ncycle, Iuhf, Iconvg, Icontl,
     &                Sing, 100, "L") 
          Call Dodiis0_l(Work, Length, Iuhf, 1, Ncycle, Iconvg,
     &                   Icontl, Sing, 144, 61, 190, 0, 90, 
     &                   2, 70, "       ")
          Call Drmove(Work, Length, Iuhf, 100, Sing)
          Call Rnabij(Work,Length,iuhf,"L")

          Call Cmpeng2(Work,Length,60,2,Ecorr,Energy(Ncycle+1,1),
     &                 Energy(Ncycle+1,2),Iuhf,1)

       Enddo 

       Return
       End
