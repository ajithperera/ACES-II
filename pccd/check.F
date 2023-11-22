      Subroutine Check(vv,oo,eval,nocc,nvrt,Nbas)

      Implicit Double Precision(a-h,o-z)
      Dimension VV(nvrt), OO(nocc),eval(nbas)

      write(6,*) (vv(i),i=1,nvrt) 
      write(6,*) (oo(i),i=1,nocc) 
      write(6,*) (eval(i),i=1,nbas) 

      Return
      End 
