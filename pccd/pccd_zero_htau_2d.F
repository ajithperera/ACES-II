      Subroutine Pccd_zero_htau_2d(Hoo_qp,Hoo_pq,Hvv_qp,Hvv_pq,Hov,Hvo,
     +                             Nocc,Nvrt,Nbas)

      Implicit Double Precision(A-H,O-Z)

      Dimension Hoo_pq(Nocc*Nocc)
      Dimension Hoo_qp(Nocc*Nocc)
      Dimension Hvv_pq(Nvrt*Nvrt)
      Dimension Hvv_qp(Nvrt*Nvrt)
      Dimension Hvo(Nvrt*Nocc)
      Dimension Hov(Nocc*Nvrt)


      Call Dzero(Hoo_pq,Nocc*Nocc)
      Call Dzero(Hoo_qp,Nocc*Nocc)
      Call Dzero(Hvv_pq,Nvrt*Nvrt)
      Call Dzero(Hvv_qp,Nvrt*Nvrt)
      Call Dzero(Hvo,Nocc*Nvrt)
      Call Dzero(Hov,Nocc*Nvrt)


      Return 
      End
