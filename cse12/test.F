
      Subroutine test(
     +                Fock_ov_a,Fock_ov_b,
     +                Fockoo_od_a,Fockoo_od_b,
     +                Fockvv_od_a,Fockvv_od_b,
     +                Nocc_a,Nocc_b,Nvrt_a,Nvrt_b,
     +                Work,Maxcor)

      Implicit Double Precision (A-H,O-Z)

#include "cse12_vars.com"

      Dimension Work(Maxcor)

      Dimension fockov_a(Nocc_a,Nvrt_a)
      Dimension fockov_b(Nocc_b,Nvrt_b)

      Dimension fockoo_od_a(Nocc_a,Nocc_a)
      Dimension fockoo_od_b(Nocc_b,Nocc_b)
      Dimension fockvv_od_a(Nvrt_a,Nvrt_a)
      Dimension fockvv_od_b(Nvrt_b,Nvrt_b)

      call check_fov(Fockov_a,Fockov_b,Nocc_a,Nocc_b,Nvrt_a,
     +                 Nvrt_b)


      Return
      End

