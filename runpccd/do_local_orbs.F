      Subroutine do_local_orbs()

      Implicit integer(A-Z)

      Call Runit("xvmol")
      Call Runit("xvmol2ja")
      Call Runit("xvscf")
      Write(6,*)
      Write(6,"(2a)") "-------------Entering Natural Local Orbital",
     &               " Generation---------------"
      Call Runit("xnlorb")
      Call Runit("xvtran")
      Call Runit("xintprc")
      Call Runit("xvcc")

      Return 
      End 
