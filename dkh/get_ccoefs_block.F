      Subroutine get_ccoefs_block(Ccoefs,Ccoef_block,Ishell,
     +                            Npcoef,Ccoef_off,Length,Nshells)

      Implicit Double Precision (A-H, O-Z)

      Dimension Ccoefs(Npcoef)
      Dimension Ccoef_block(Length)
      Integer Ccoef_off(Nshells)

      ioff = Ccoef_off(Ishell) 

#ifdef _DEBUG_LVLM
      write(6,*)
      write(6,"(a,(1x,i4))")"Ishell",Ishell
      write(6,"(a,(1x,i4))")"Ioff",ioff
      write(6,*) "The contraction coefs."
      write(6,"(6(1x,f10.6))") (CCoefs(i),i=1,npcoef)
      Write(6,*)
#endif 
      Call Dcopy(Length,Ccoefs(ioff),1,Ccoef_block,1)
      

#ifdef _DEBUG_LVLM
      Write(6,*) "Length", Length
      Write(6,"(a,2i4)") "The contraction coefs. of",ishell
      write(6,"(6(1x,E12.5))") (CCoef_block(i),i=1, Length)
#endif 
      
      Return
      End 
     
