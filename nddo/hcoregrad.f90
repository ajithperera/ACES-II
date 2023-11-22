subroutine hcoregrad(loop)
   use constants
use tables
use indices
implicit double precision (a-h,o-z)
integer,dimension(2),intent(in)::loop
double precision,dimension(:),allocatable::htemp

ir=0
jc=0
ir2=0
jc2=  0
!nelectrons=0
! loopf 

do m=1,2
i=loop(m)

   ni=nbas(species(i))
lmax=10
if(ni>4)lmax=45
      istart=ifirst(i)
      iend=ilast(i)
      istart2=ifirst2(i)
      iend2=ilast2(i)           
if(ni==1)then
      t1=zero
      do n=1,2
         j=loop(n)
         if(j==i)goto 1
call pack2(ifirst2(i),ifirst2(j),ij) ! returns ij composite index b/c here i may be greater than j
t1=t1+eff_core(zeff(j))*twoe(ij)
1        continue
      end do
H(istart+offset1(istart))=-t1



   else
       
      allocate(htemp(lmax))  
         htemp=zero
         do kk=1,2
          j=loop(kk)
               if(j==i)goto 2
         charge=eff_core(zeff(j))
         
           do k=0,lmax-1
call pack2(ifirst2(i)+k,ifirst2(j),ij)
 htemp(k+1)=htemp(k+1)+charge*twoe(ij)
          end do
 
2      continue
      end do
H(istart+offset1(istart))=-htemp(1)
H(offset1(istart+1)+istart)=-htemp(2)
H(offset1(istart+2)+istart)=-htemp(3)
H(offset1(istart+3)+istart)=-htemp(4)
if(ni>4)then

H(offset1(istart+4)+istart)=-htemp(11)
H(offset1(istart+5)+istart)=-htemp(22)
H(offset1(istart+6)+istart)=-htemp(37)
H(offset1(istart+7)+istart)=-htemp(16)
H(offset1(istart+8)+istart)=-htemp(29)
end if

i1=istart+1

H(i1+offset1(i1))=-htemp(5)
H(offset1(istart+2)+i1)=-htemp(6)
H(offset1(istart+3)+i1)=-htemp(7)


if(ni>4)then

H(offset1(istart+4)+i1)=-htemp(12)
H(offset1(istart+5)+i1)=-htemp(23)
H(offset1(istart+6)+i1)=-htemp(38)
H(offset1(istart+7)+i1)=-htemp(17)
H(offset1(istart+8)+i1)=-htemp(30)

end if
i1=istart+2
H(offset1(i1)+i1)=-htemp(8)
H(offset1(istart+3)+i1)=-htemp(9)
if(ni>4)then

H(offset1(istart+4)+i1)=-htemp(13)
H(offset1(istart+5)+i1)=-htemp(24)
H(offset1(istart+6)+i1)=-htemp(39)
H(offset1(istart+7)+i1)=-htemp(18)
H(offset1(istart+8)+i1)=-htemp(31)
end if
i1=istart+3

H(i1+offset1(i1))=-htemp(10)

if(ni>4)then

H(offset1(istart+4)+i1)=-htemp(14)
H(offset1(istart+5)+i1)=-htemp(25)
H(offset1(istart+6)+i1)=-htemp(40)
H(offset1(istart+7)+i1)=-htemp(19)
H(offset1(istart+8)+i1)=-htemp(32)

i1=istart+4


H(i1+offset1(i1))=-htemp(15)
H(offset1(istart+5)+i1)=-htemp(26)
H(offset1(istart+6)+i1)=-htemp(41)
H(offset1(istart+7)+i1)=-htemp(20)
H(offset1(istart+8)+i1)=-htemp(33)


i1=istart+5
H(i1+offset1(i1))=-htemp(28)
H(offset1(istart+6)+i1)=-htemp(43)
H(offset1(istart+7)+i1)=-htemp(27)
H(offset1(istart+8)+i1)=-htemp(35)

i1=istart+6

H(i1+offset1(i1))=-htemp(45)
H(offset1(istart+7)+i1)=-htemp(42)
H(offset1(istart+8)+i1)=-htemp(44)




i1=istart+7

H(i1+offset1(i1))=-htemp(21)
H(offset1(istart+8)+i1)=-htemp(34)


i1=istart+8
H(i1+offset1(i1))=-htemp(36)
end if













deallocate(htemp)
end if


end do





end subroutine hcoregrad
