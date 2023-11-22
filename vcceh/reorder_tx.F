      Subroutine Reorder_tx(Tx,Work,Maxcor,T1ln,T2ln,T1lnn_aa,T1ln_bb,
     +                      T2ln_aa,T2ln_bb,T2ln_ab,Nsizec,Iuhf)
      
      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Dimension Tx(Nsizec)

      Integer T2ln_aa,T2ln_bb,T2ln_ab,T2ln
      Integer T1ln,T1ln_aa,T1ln_bb

      Data Ione /1/

C Incoming Tx is ordered ad T1(AA)-T1(BB)-T2(ABAB)-T2(AAAA)-T2(BBBB)
C and it is reordered to T1(AA)-T1(BB)-T2(AAAA)-T2(BBBB)-T2(ABAB)


      Len = T2ln_aa + Iuhf*T2ln_bb + T2ln_ab
      I000 = Ione 
      Iend = I000 + Len

      If (Iend .Gt. Maxcor) Call Insmem("reorder_tx",Iend,Maxcor)
     
      Ioff = T1ln_aa + Iuhf*T1n_bb + Ione

      Call Dcopy(Len,Tx(Ioff),1,Work(I000),1)

      Joff = T2ln_ab
      Ioff = T1ln_aa + Iuhf*T1ln_bb + Ione
      Len  = T2ln_aa + Iuhf*T2ln_bb 

      Call Dcopy(Len,Work(I000+Joff),1,Tx(Ioff),1)
      Ioff = T1ln_aa + Iuhf*T1ln_bb + T2ln_aa + Iuhf*T2ln_bb + Ione 
      
      Call Dcopy(Tln_ab,Work(I000),1,Tx(Ioff),1)

      Return
      End 
