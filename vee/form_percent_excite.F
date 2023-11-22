
#include "flags.h"

      SUBROUTINE FORM_PERCENT_EXCITE(SCR, MAXCOR, IRREPX, IUHF,
     &                               ROOT, ECC, IROOT) 
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DOUBLE PRECISION SCR(MAXCOR)

#include "flags2.com"
#include "flags.com"
#include "syminf.com"
#include "info.com"
#include "machsp.com"

      Data Cen /1.0D+02/

      Nbas = Nocco(1) + Nvrto(1)
  
      Iunit  = 556

C The unit 556 is opened in  vee.F and correspond to "EOM_excite"

      Itop   = 1
      Imap_a = Itop 
      Imap_b = Imap_a + Nbas
      Iend   = Imap_b + Nbas
      If (Iend .GE. Maxcor) CALL INSMEM("@-FORM_PERCENT_EXCITE",Iend,
     &                                   Maxcor)
      
      Call Getrec(20, "JOBARC", "ORBMAP_A", Nbas, Scr(Imap_a))
      Call Getrec(20, "JOBARC", "ORBMAP_A", Nbas, Scr(Imap_b))
      If (Iuhf .Ne. 0) Call Getrec(20, "JOBARC", "ORBMAP_B", Nbas, 
     &                             Scr(Imap_b))

      Write(Iunit,*)
      Write(Iunit,881) iroot, irrepx
 881  Format(' Root # ', i4, 3x, ' [', i1,']')
      Write(iunit,882) Root*27.2113957d0, Root, Root+Ecc
 882  Format( ' Excitation energy (eV): ', F12.6, ' in a.u. ', F14.8,
     &      /,  ' Total energy : ', F16.8)

      Call get_topamps(Scr(Iend),Maxcor-Iend,Scr(Imap_a),Scr(Imap_b),
     &                 Iuhf,Nbas,Irrepx,Iunit,R1,R2)

C Compute the %of R1 and R2 
      R = R1 + R2
      R1 = (R1/R)*Cen
      R2 = (R2/R)*Cen
      Write(6,"(7x,a)") "-------Percent R1 and R2-----------------"
      Write(6,"(16x,a,F5.1)") "%R1 = ", R1
      Write(6,"(16x,a,F5.1)") "%R2 = ", R2
      Write(6,"(7x,a)") "-----------------------------------------"
      Write(6,*)

      Return
      End
     
