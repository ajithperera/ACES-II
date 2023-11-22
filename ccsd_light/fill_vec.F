      Subroutine Fill_vec(Vec,Nocc,Nvrt,Iloc,Ioff)

      Implicit Double Precision (A-H,O-Z)

      Dimension Vec(Nvrt*Nvrt*Nocc*Nocc)
      Integer I,J,A,B

      Data Dzero,Done,Izero/0.0D0,1.0D0,0/

      Write(6,"(a)") " @-Fill_vec-in"
      call output(Vec,1,Nvrt*Nvrt,1,Nocc*Nocc,Nvrt*Nvrt,Nocc*Nocc,1)
      Idiag = Ioff
      Icount = Izero 
      Do I = 1, Nocc
      Do J = 1, Nocc
      Do A = 1, Nvrt
      Do B = 1, Nvrt
         Write(6,"(4(1x,I2))") I,J,A,B
         If (I .Eq. J .OR. A .Eq. B) Then
         Icount = Icount + 1
         Idiag  = Idiag  + 1
         Write(6,"(3(1x,I2))") Idiag,Icount,Iloc
         If (Idiag .Eq. Iloc) Vec(Icount) = Dzero
         Else
             Idiag = Idiag + 1
             Icount = Icount + 1
         Endif 
      Enddo 
      Enddo 
      Enddo 
      Enddo 
 10   Continue
      Write(6,"(a)") " @-Fill_vec-ot"
      call output(Vec,1,Nvrt*Nvrt,1,Nocc*Nocc,Nvrt*Nvrt,Nocc*Nocc,1)

      Return 
      End 
