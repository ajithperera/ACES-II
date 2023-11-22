subroutine mkd(rhor,x1,y1,z1)
use constants
use tables
use scratch_array
use indices
implicit double precision (a-h,o-z)
!double precision,dimension(:,:),intent(in):: density
double precision,intent(inout)::rhor
double precision,intent(in)::x1,y1,z1
double precision,dimension(num_basis)::amplit(num_basis)
parameter(pi3=3.0d0*3.141592653589793d0)
rhor=0.0d0
icount=0
do i=1,numat
ni=nbas(species(i))
x_relative=(x1-x(i))
y_relative=(y1-y(i))
z_relative=(z1-z(i))
rsq=x_relative*x_relative+y_relative*y_relative+z_relative*z_relative
r=dsqrt(rsq)
 do j=1,ni
icount=icount+1
 nl=j
! load s data
 if(j.eq.1)then
 ei=zs(species(i))
 eir=-ei*r
 nq=nqs(zeff(i))
! load p data
 elseif(j>1 .and. j<5)then
 ei=zp(species(i))
 eir=-ei*r
 nq=nqp(zeff(i))
! load d data
 else
 ei=zetad(species(i))
 eir=-ei*r
 nq=nqd(zeff(i))
 end if
! 1s  
      if(nq.eq.1)then
      amplit(icount)=dsqrt(ei*ei*ei/pi)*exp(eir)
! 2s
      elseif(nq.eq.2.and.nl.eq.1)then
      amplit(icount)=dsqrt(ei*ei*ei*ei*ei/(pi3))*r*exp(eir)
! 2px
      elseif(nq.eq.2.and.nl.eq.2)then
      amplit(icount)=x_relative*dsqrt(ei*ei*ei*ei*ei/pi)*exp(eir)
! 2py
      elseif(nq.eq.2.and.nl.eq.3)then
      amplit(icount)=y_relative*dsqrt(ei*ei*ei*ei*ei/pi)*exp(eir)
! 2pz
      elseif(nq.eq.2.and.nl.eq.4)then
      amplit(icount)=z_relative*dsqrt(ei*ei*ei*ei*ei/pi)*exp(eir)
! 3s
      elseif(nq.eq.3.and.nl.eq.1)then
      amplit(icount)=.118941607d0*dsqrt(ei*ei*ei*ei*ei*ei*ei)*exp(eir)*r*r
! 3px
      elseif(nq.eq.3.and.nl.eq.2)then
       amplit(icount)=.206012907d0*dsqrt(ei*ei*ei*ei*ei*ei*ei)*exp(eir)*r*x_relative
! 3py
      elseif(nq.eq.3.and.nl.eq.3)then
       amplit(icount)=.206012907d0*dsqrt(ei*ei*ei*ei*ei*ei*ei)*exp(eir)*r*y_relative
! 3pz
      elseif(nq.eq.3.and.nl.eq.4)then
       amplit(icount)=.206012907d0*dsqrt(ei*ei*ei*ei*ei*ei*ei)*exp(eir)*r*z_relative
       end if
         





!         rhor=rhor+amplit(icount)*amplit(icount)*s3(icount,icount)
       rhor=rhor+amplit(icount)*amplit(icount)*s3(icount+offset1(icount))
      end do
      end do
    do k=1,numat
       if(nbas(species(k)).gt.1)then 
          istart=ifirst(k)

          iend=ilast(k)
          do i=istart+1,iend
          do j=istart,i-1
!          rhor=rhor+amplit(i)*amplit(j)*s3(j,i)*two
             call pack1(i,j,ij)
             rhor=rhor+amplit(i)*amplit(j)*s3(ij)*two
          end do
          end do
        end if
        end do

!      do i=2,num_basis
!      do j=1,i-1
!          rhor=rhor+amplit(i)*amplit(j)*s3(j,i)*two
!         end do
!         end do
   end subroutine mkd
