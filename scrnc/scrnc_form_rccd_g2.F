      Subroutine Scrnc_form_rccd_g2(Work,Maxmem,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxmem)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "sym.com"
#include "flags2.com"
#include "flags.com"
#include "info.com"

      Data One /1.0D0/

C Form G(a,e)=-1/2 T2(mn,ef)L2(af,mn) and G(i,j)=1/2 T2(mn,ef)L2(ef,in)
C For rCCD and DrCCD, these were not created in ringcl code. In principle 
C these could be moved there (and add capablities for rCCD and drCCD 
C property calculations).

      Call formg1(Work,Maxmem,Iuhf,One)
      Call formg2(Work,Maxmem,Iuhf,One)

      Return
      End 
