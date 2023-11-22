      Subroutine check_fod(Fockoo_od_a,Fockoo_od_b,Fockvv_od_a,
     +                     Fockvv_od_b,Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision (A-H,O-Z)
  
      Dimension Fockoo_od_a(Nocc_a,Nocc_a),Fockoo_od_b(Nocc_b,Nocc_b)
      Dimension Fockvv_od_a(Nvrt_a,Nvrt_a),Fockcvv_od_b(Nvrt_b,Nvrt_b)

      Write(6,*) "The OO-od"
      call output(Fockoo_od_a,1,Nocc_a,1,Nocc_a,Nocc_a,Nocc_a,1)
      call output(Fockoo_od_b,1,Nocc_b,1,Nocc_b,Nocc_b,Nocc_b,1)
      Write(6,*) "The vv-od"
      call output(Fockvv_od_a,1,Nvrt_a,1,Nvrt_a,Nvrt_a,Nvrt_a,1)
      call output(Fockvv_od_b,1,Nvrt_b,1,Nvrt_b,Nvrt_b,Nvrt_b,1)

      Return
      End 


