      subroutine gen_mis_grad(Igrad,Fgrad,natom,nt)
      implicit none
      integer natom,i,nt,offset 
      double precision Igrad(nt-3),Fgrad(nt)

      do i =1,natom-1
        offset=3*i
        Fgrad(1)=Fgrad(1)+Igrad(offset-2)
        Fgrad(2)=Fgrad(2)+Igrad(offset-1)
        Fgrad(3)=Fgrad(3)+Igrad(offset)
      end do     
      Fgrad(1)=-Fgrad(1)
      Fgrad(2)=-Fgrad(2)
      Fgrad(3)=-Fgrad(3) 
      return
      end  
