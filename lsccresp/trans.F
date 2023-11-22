      Subroutine Trans(T_in,T_out,Coo,Cvv,Nrow,Ncol,Type)

      Implicit Double Precision(A-H,O-Z)
      Character*1 A,B,C,D
      Character*3 Type

      Dimension T_in(Nrow,Ncol)
      Dimension T_out(Nrow,Ncol)
      Dimension Coo(Ncol,Ncol) 
      Dimension Cvv(Nrow,Nrow)
       
      Data Dnull,Done /0.0D0,1.0D0/

      If (Type .Eq. "C2N") Then
          A = "T"
          B = "N"
          C = "N"
          D = "N"
      elseif (Type .Eq. "N2C") Then
          A = "N"
          B = "N"
          C = "N"
          D = "T"
      Endif

#ifdef _DEBUG_LVL0
CSSS      Write(6,"(a)") " At entry T1 amplitudes"
CSSS      call output(T_in,1,Nrow,1,Ncol,Nrow,Ncol,1)
      call checksum("T1-in:",T_in,nrow*Ncol)
CSSS      Write(6,"(a)") " The occ-occ block of natural orbitals"
CSSS      call output(Coo,1,Ncol,1,Ncol,Ncol,Ncol,1)
      call checksum("COO  :",Coo,ncol*ncol)
CSSS      Write(6,"(a)") " The vrt-vrt block of natural orbitals"
CSSS      call output(Cvv,1,Nrow,1,Nrow,Nrow,Nrow,1)
      call checksum("Cvv  :",Cvv,nrow*nrow)
#endif
C T1(a,i) -> T1(a',i)

      Call Xgemm(A,B,Nrow,Ncol,Nrow,Done,Cvv,Nrow,T_in,Nrow,
     +            Dnull,T_out,Nrow)

C T1(a',i) -> T1(a',i')

      Call Xgemm(C,D,Nrow,Ncol,Ncol,Done,T_out,Nrow,Coo,Ncol,
     +            Dnull,T_in,Nrow)

#ifdef _DEBUG_LVL0
CSSS      Write(6,"(a)") " At exit NO transformed T1 amplitudes"
CSSS      call output(T_out,1,Nrow,1,Ncol,Nrow,Ncol,1)
      call checksum("T1-ot:",T_in,nrow*Ncol)
#endif 

      Return
      End 
