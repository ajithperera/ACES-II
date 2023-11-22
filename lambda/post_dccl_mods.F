#include "flags.h"
      Subroutine post_dccl_mods(Work, Length, Iuhf)

      Implicit Double Precision (A-H, O-Z)

      Dimension Work(Length)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "flags.com"

      Imode  = 0
      Irrepx = 1
C
C Needs the G-intermediates scales by Gae_ and G_mi scales.
C
      IRREPX      = 1
      IHHA_LENGTH = IRPDPD(IRREPX,21)
      IHHB_LENGTH = IRPDPD(IRREPX,22)
      IPPA_LENGTH = IRPDPD(IRREPX,19)
      IPPB_LENGTH = IRPDPD(IRREPX,20)

      Call Updmoi(1,IHHA_LENGTH,9,191, 0, 0)
      Call Updmoi(1,IPPA_LENGTH,9,192, 0, 0)
      If (Iuhf .ne.0) Then
         Call Updmoi(1,IHHB_LENGTH,10,191, 0, 0)
         Call Updmoi(1,IPPB_LENGTH,10,192, 0, 0)
      Endif

      Call zerlst(Work,IHHA_LENGTH,1,1,9,191)
      Call zerlst(Work,IPPA_LENGTH,1,1,9,192)
      If (iuhf .ne. 0) then
            call zerlst(Work,IHHB_LENGTH,1,1,10,191)
            call zerlst(Work,IPPB_LENGTH,1,1,10,192)
      Endif
CSSS      If (Coulomb) Then
CSSS         Call Pdcc_formg1(Work,Length,Iuhf,1.0D0,8)
CSSS      Else 
CSSS        Call Pdcc_formg1(Work,Length,Iuhf,Gae_scale,8)
CSSS      Endif 
      Call Pdcc_formg1(Work,Length,Iuhf,Gae_scale,8)
      Call Pdcc_formg2(Work,Length,Iuhf,Gmi_scale,8)

C Store a flag to indicate that lambda has been done. This flag
C is used in EOM code to test whether lambad is done. 

      Call Putrec(20,"JOBARC","LAMBDA  ",1,Length)

#ifdef _DEBUG_LVLM
CSSS      Write(6,*)
CSSS      Write(6,"(a)") " Antisymmetry check of T2 and L2 Amps."
CSSS      Call Modf_T2L2_inL(Work,Length,Iuhf)
#endif 

      Return
      End
 
