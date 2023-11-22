      Subroutine Run_freq(Calc_level, Analytical_grads)
c
      Implicit None
c
#include "flags.h"
c
      Character*4 Calc_level
      Character*79 szGExtrap  
      Logical Analytical_grads
      Integer Istat, Is, Ius, Irank
c
      iRank=0
      istat=0
      Do while (istat.eq.0)
         Call Runit("xvmol")
         Call Runit("xvmol2ja")
         Call Runit("xvscf")
         Call Runit("xvtran")
         Call Runit("xintprc")
         Call Runit("xvcc")
         Call Runit("xlambda")
         Call Runit("xdens")
         Call Runit("xanti")
         Call Runit("xbcktrn")
         Call Runit("xvdint")
c
         Call a2getrec(1,'JOBARC','LASTGEOM',1,istat)
   
         If (istat.eq.0) then
            Call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
            Call runit("xjoda")
         Else
           If (irank.eq.0) then
              Call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
              Call runit("xjoda")
           Endif
        Endif  
C
        Call rmfiles
c
      End do
c
      Return
      End
