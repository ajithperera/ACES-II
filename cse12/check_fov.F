      Subroutine check_fov(Fockov_a,Fockov_b,Nocc_a,Nocc_b,Nvrt_a,
     +                     Nvrt_b)

      Implicit Double Precision (A-H,O-Z)
  
      Dimension Fockov_a(Nocc_a,Nvrt_a),Fockov_b(Nocc_b,Nvrt_b)

      Write(6,*) "The Fock-OV-A"
      call output(Fockov_a,1,Nocc_a,1,Nvrt_a,Nocc_a,Nvrt_a,1)
      Write(6,*) "The Fock-OV-b"
      call output(Fockov_b,1,Nocc_b,1,Nvrt_b,Nocc_b,Nvrt_b,1)

      Return
      End 


