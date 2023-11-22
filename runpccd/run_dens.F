      Subroutine Run_Dens()

      Implicit Integer(A-Z)

#include "flags.h"
 
      Call Runit("xvmol")
      Call Runit("xvmol2ja")
      Call Runit("xvprops")
      Call Runit("xvscf")
      Call Runit("xvtran")
      Call Runit("xintprc")
      Call Runit("xvcc")
      Call Runit("xlambda")
      Call Runit("xdens")
      Call Runit("xprops")

      Return
      End 
