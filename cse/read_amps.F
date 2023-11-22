      Subroutine Read_amps(Iunit,T1aa,T1bb,Nocc_a,Nocc_b,Nvrt_a,
     +                     Nvrt_b,Work,Maxcor)

      Implicit Double Precision(A-H,O-Z)

      Dimension T1aa(Nvrt_a,Nocc_a)
      Dimension T1bb(Nvrt_b,Nocc_b)

      Dimension Work(Maxcor)

      Integer I,J,A,B

      Write(6,"(2a)") "The starting amplitudes are read from the",
     +                " TGUESS file."
 
      Do I = 1, Nocc_a
      Do A = 1, Nvrt_a
         Read(Iunit) T1aa(A,I)
      Enddo
      Enddo

      Do I = 1, Nocc_b
      Do A = 1, Nvrt_b
         Read(Iunit) T1bb(a,i)
      Enddo
      Enddo

#ifdef _DEBUG_LVLM
C This is an invarience check. What happens here is that we take
C converged T1 and T2 vectors and then add a diagonal matrix to it. The
C Diagonals of this matrix can be any value > 0 and less than 1. 
C However they should be very close to the original denisty matrices 
C that we used to obtain converged T1 and T2. Otherwise it might not
C Converge. 
 
      Id1  = 1
      Id2  = Id1 + Nvrt_a*Nocc_a
      Iend = Id2 + Nvrt_b*Nocc_b
      Call Augment_ts(T1aa,T1bb,Nocc_a,Nocc_b,Nvrt_a,vrt_b,Work,
     +                Work(Id1),Work(Id2))
#endif 
#ifdef _DEBUG_LVL0
      do i=1,Nocc_a
      do a=1,Nvrt_a
      Write(6,"(2I2,1x,F15.10)") a,i,T1aa(a,i)
      enddo
      enddo
#endif 

      Return
      End
