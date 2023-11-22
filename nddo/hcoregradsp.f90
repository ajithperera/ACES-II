subroutine hcoregradsp(i)
   use constants
use tables
use indices
implicit double precision (a-h,o-z)
integer,intent(in)::i
double precision,dimension(:),allocatable::htemp
double precision,dimension(:,:),allocatable::wout
interface
 subroutine twoe_sparkle(i,j,wout)
  double precision,dimension(:,:),intent(inout)::wout
integer,intent(in)::i,j
end subroutine twoe_sparkle
end interface


   ni=nbas(species(i))
lmax=10
if(ni>4)lmax=45
      istart=ifirst(i)
      iend=ilast(i)
      istart2=ifirst2(i)
      iend2=ilast2(i)           
if(ni==1)then
H(istart+offset1(istart))=zero
   allocate(wout(pairs(species(i)),1))
   do j=1,nsparkle
   call twoe_sparkle(i,j,wout)
   H(istart+offset1(istart))=H(istart+offset1(istart))-spcharge(zeffsp(j))*wout(1,1)
   end do
   deallocate(wout)





   else
       
      allocate(htemp(lmax))  
         htemp=zero

            allocate(wout(pairs(species(i)),1))
   do j=1,nsparkle
   call twoe_sparkle(i,j,wout)

   do k=0,lmax-1
 htemp(k+1)=htemp(k+1)+spcharge(zeffsp(j))*wout(k+1,1)
          end do
   end do
deallocate(wout)




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








end subroutine hcoregradsp
