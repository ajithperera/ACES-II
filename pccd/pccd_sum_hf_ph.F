      Subroutine Pccd_sum_hf_ph(Work,Hvo,Nk,Ni,Na,Nocc,Nvrt)

      Implicit Double Precision(A-H,O-Z)
      Integer A 

      Dimension Hvo(Nvrt,Nocc)
      Dimension Work(Nk,Nk,Ni,Na)
      Data Dnull /0.0D0/

      Do A = 1, Na
      Do I = 1, Ni
        Sum = DNull
        Do N = 1, Nk
           Print*, Work(N,N,I,A)
           Sum = Sum + Work(N,N,I,A)
        Enddo
        Hvo(A,I) = Hvo(A,I) + Sum*2.0D0
      Enddo 
      Enddo

      Return
      ENd 
  
     
