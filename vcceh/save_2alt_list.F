      Subroutine Save_2alt_list(Work,Maxcor,T1ln,T2ln,T1ln_aa,
     +                          T1ln_bb,T2ln_aa,T2ln_bb,T2ln_ab,
     +                          Irrepx,Iuhf)

      Implicit Double Precision(A-H,O-Z)
      Integer T1ln,T2ln
      Integer T1ln_aa,T1ln_bb,T2ln_aa,T2ln_bb,T2ln_ab

      Dimension Work(Maxcor)
 
      Data Ione /1/

      List1_h = 490
      List1_t = 493

      I000 = Ione
      Call Getlst(Work(I000),1,1,1,1,List1_h)
      Call Putlst(Work(I000),1,1,1,1,List1_t)

      If (Iuhf .Ne. 0) Then
         Call Getlst(Work(I000),1,1,1,2,List1_h)
         Call Putlst(Work(I000),1,1,1,2,List1_t)
      Endif

      List2_h = 444
      List2_t = 454

      Call Getall(Work(I000),T2ln_aa,Irrepx,List2_h)
      Call Putall(Work(I000),T2ln_aa,Irrepx,List2_t)

      If (Iuhf .Ne. 0) Then
         Call Getall(Work(I000),T2ln_bb,Irrepx,List2_h+1)
         Call Putall(Work(I000),T2ln_bb,Irrepx,List2_t+1)
      Endif 

      Call Getall(Work(I000),T2ln_ab,Irrepx,List2_h+2)
      Call Putall(Work(I000),T2ln_ab,Irrepx,List2_t+2)

#ifdef _DEBUG
      List1 = 493
      List2 = 453

      Call Getlst(Work(I000),1,1,1,1,List1)
      call checksum("T1aa :",Work(I000),T1ln_aa)
      Ioff = T1ln_aa
      If (Iuhf .Ne. 0) Call Getlst(Work(I000+Ioff),1,1,1,2,List1)
      call checksum("T1bb :",Work(I000+ioff),T1ln_bb)

      Ioff = T1ln_aa + Iuhf*T1ln_bb
      Call Getall(Work(I000+Ioff),T2ln_aa,Irrepx,List2+1)
      call checksum("T2aa :",Work(I000+ioff),T2ln_aa)
      Ioff = T1ln_aa + Iuhf*T1ln_bb + T2ln_aa
      If (Iuhf .Ne. 0) Call Getall(Work(I000+Ioff),T2ln_bb,
     +                             Irrepx,List2+2)
      call checksum("T2bb :",Work(I000+ioff),T2ln_bb)
      Ioff =  T1ln_aa + Iuhf*T1ln_bb + T2ln_aa + Iuhf*T2ln_bb
      Call Getall(Work(I000+Ioff),T2ln_ab,Irrepx,List2+3)
      call checksum("T2ab :",Work(I000+ioff),T2ln_ab)
#endif


      Return
      End
      

