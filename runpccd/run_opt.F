      Subroutine Run_Opt()
#include "flags.h"

      integer icore, is, ius
      logical Do_derint

      ndrgeo=0
      szGExtrap = 'xa2proc grad_extrp'
      call a2getrec(20,'JOBARC','NDROPGEO',1,NDRGEO)

      icycle=0
 1000 continue
      icycle=icycle+1
c
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
      call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
      call runit('xjoda')
      call a2getrec(1,'JOBARC','JODADONE',1,istat)
      if (istat.ne.0) return
      call rmfiles
      goto 1000

      Return
      End

