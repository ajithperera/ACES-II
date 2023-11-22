      Subroutine Dump_amps(T1aa,T1bb,T2aa,T2bb,T2ab,Nocc_a,
     +                     Nocc_b,Nvrt_a,Nvrt_b)

      Implicit Double Precision(A-H,O-Z)

      Dimension T2aa(Nvrt_a,Nvrt_a,Nocc_a,Nocc_a)
      Dimension T2bb(Nvrt_b,Nvrt_b,Nocc_b,Nocc_b)
      Dimension T2ab(Nvrt_a,Nvrt_b,Nocc_a,Nocc_b)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Integer I,J,A,B
      Logical T_guess_exist
      Character*80 Fullname 

C Write the coverged amplitudes to TGUESS file to be used in 
C another calculation

      Call Gfname("TGUESS",Fullname,Namelen)
      Inquire(File=Fullname(1:Namelen),Exist=T_guess_exist)

      If (T_guess_exist) Then
         Open(Unit=94,File=Fullname(1:Namelen),Status="Old",
     +        Form="Unformatted",Access="Sequential")
         Close(Unit=94,Status="Delete") 
      Endif
      Open(Unit=94,File=Fullname(1:Namelen),Status="New",
     +     Form="Unformatted",Access="Sequential")
 
      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
         Write(94) T1aa(A,I)
      Enddo
      Enddo

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
         Write(94) T1bb(a,i)
      Enddo
      Enddo

      Do J = 1, Nocc_a
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_a
      Do A = 1, Nvrt_a
         Write(94) T2aa(A,B,I,J)
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do I = 1, Nocc_b
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_b
         Write(94) T2bb(a,b,i,j)
      Enddo
      Enddo
      Enddo
      Enddo

      Do J = 1, Nocc_b
      Do I = 1, Nocc_a
      Do B = 1, Nvrt_b
      Do A = 1, Nvrt_a
         Write(94) T2ab(A,b,I,j)
      Enddo
      Enddo
      Enddo
      Enddo

      Close(Unit=94,Status="Keep")

      Return
      End
