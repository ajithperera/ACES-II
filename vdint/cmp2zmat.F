      Subroutine Cmp2zmat(CMP_GRAD,ZMT_GRAD,IMAP,NAMES_REORD,NATOMS)

      Implicit Double Precision(A-H,O-Z)
     
      Dimension CMP_GRAD(3,NATOMS)
      Dimension ZMT_GRAD(3,NATOMS)
      Dimension IMAP(NATOMS)
      Character*6 NAMES_REORD(NATOMS)
      Character*5 AtmLabel(Natoms)
      Character*5 String
      Character*1 V1
      Character*2 V2
      Character*2 V3

      Call Zero(ZMT_GRAD,3*Natoms)
      Do Iatmvml = 1, Natoms
         Iatmzmt = Imap(Iatmvml)
         If (Iatmzmt .Ne. 0) Then
            ZMT_GRAD(1,Iatmzmt)  = Cmp_Grad(1,Iatmvml)
            ZMT_GRAD(2,Iatmzmt)  = Cmp_Grad(2,Iatmvml)
            ZMT_GRAD(3,Iatmzmt)  = Cmp_Grad(3,Iatmvml)
         Endif 
      Enddo

      Call Getcrec(1,'JOBARC','ZSYM',5*Natoms,AtmLabel)
      Do I = 1, Natoms 
         If (I .Lt. 10) Then
             Write(V1,"(I1)") I
             V2 = "0"//V1
             String = AtmLabel(I)
             Names_reord(I) = String(1:2)//V2
         Else if (I .Lt. 99) Then
             Write(V3,"(I2)") I
             String = AtmLabel(I)
             Names_reord(I) = String(1:2)//V3
         Else
             Write(6,"(a)") "@-prigrd->cmp2zmat, More than 100",
     +                      " atoms is not allowed"
             Call Errex
         Endif
      Enddo

      Return
      End 


