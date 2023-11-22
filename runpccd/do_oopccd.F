#include "flags.h"
      Subroutine Do_oopccd(Iuhf)

      integer, intent(in)::Iuhf
      Logical Ca,Fc,Bca,Bfc,Occnums_file,Move
      Character*80 Fname

#include "flags.com"
#include "flags2.com"
#include "sympop.com"

      Data Izz /100/

      print*,'inside job control',Iflags2((h_IFLAGS2_oo_method))
      print*,'Iflags for oo:', Iflags2(h_IFLAGS2_oo_method)
      print*,'true statement:',Iflags2(h_IFLAGS2_oo_method).eq.3

      If (Iflags2(h_IFLAGS2_oo_method) .EQ. 1 .or.
     +         Iflags2(h_IFLAGS2_oo_method) .EQ. 2  .or.
     +         Iflags2(h_IFLAGS2_oo_method).EQ.3    .or.
!SteepestDescent
     +         Iflags2(h_IFLAGS2_oo_method).EQ.4) then

         print*,'Running steepest descent (SD) algo to optimize...'
         call do_FullNR(Work,Maxcor,Iuhf)
         !Call do_steepestDescent(Work,Maxcor,Iuhf)
!      Else if (Iflags(h_IFLAGS2_oo_method) .EQ. 2) then !L-BFGS
!          print*,'Inside L-BFGS'
!         Call do_LBFGS()
!      Else if (Iflags2(h_IFLAGS2_oo_method).EQ.3) then !NR step
!          print*,'Running full NR algo to optimize...'
!         call do_FullNR(Work,Maxcor,Iuhf)
      Else
         print*,'zero off-diag T2'
         Call Runit("xvmol")
         Call Runit("xvmol2ja")
         Call Runit("xvscf")
         Call Runit("xvtran")
         Call Runit("xintprc")
         Call Runit("xpccd")
      Endif 

      Return
      End 


     
