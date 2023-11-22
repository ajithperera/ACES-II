      Subroutine Pccd_form_hinv(H,Work,Maxcor,Nbas,Len,Nocc,Nvrt,
     +                          Tag)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Dimension H(Len,Len)
      Character*2 Tag

      Data One,Dnull,Ione/1.0D0,0.0D0,1/

      Lensq = Len*Len
      I000  = Ione
      I010  = I000 + Lensq
      Iend  = I010 + Lensq

      If (Iend .Gt. Maxcor) Call Insmem("pccd_form_hinv",Iend,Maxcor)

      Call Dcopy(Lensq,H,1,Work(I000),1)
      Call Minv(H,Len,Len,Work(I010),Det,1.0D-8,0,1)

#ifdef _DEBUG_LVL0
      If (Tag .Eq. "OV") Write(6,"(a)") " The OV/VO inverse check"
      If (Tag .Eq. "OO") Write(6,"(a)") " The OO inverse check"
      If (Tag .Eq. "VV") Write(6,"(a)") " The VV inverse check"
      Call Dgemm("N","N",Len,Len,Len,One,H,Len,Work(I000),Len,Dnull,
     +            Work(Iend),Len)
      Call output(Work(Iend),1,Len,1,Len,Len,Len,1)
#endif

      Return
      End 
          
 
