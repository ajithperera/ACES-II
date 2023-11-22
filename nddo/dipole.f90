subroutine dipole
use tables
use constants
use indices
use scratch_array
use control
implicit double precision(a-h,o-z)
double precision,dimension(3)::center
double precision,dimension(:),allocatable::charge
factor=.577350269d0
total=zero
center=zero
allocate(charge(numat))
do i=1,numat
    ! center of gravity of coordinates
    weight=atmass(zeff(i))
    total=total+weight
    center(1)=center(1)+weight*x(i)
    center(2)=center(2)+weight*y(i)
    center(3)=center(3)+weight*z(i)
    ! compute atomic charges
    pop=zero
    k=ifirst(i)
    l=ilast(i)
    do j=k,l
!    pop=pop+s3(j,j)
       pop=pop+s3(j+offset1(j))
    end do
    charge(i)=eff_core(zeff(i))-pop
end do

center=center/total
elecx=zero
elecy=zero
elecz=zero
ptchgx=zero
ptchgy=zero
ptchgz=zero


do i=1, numat

! compute pt chg contribution to dipole
ptchgx=ptchgx+4.80324D0*charge(i)*(x(i)-center(1))
ptchgy=ptchgy+4.80324D0*charge(i)*(y(i)-center(2))
ptchgz=ptchgz+4.80324D0*charge(i)*(z(i)-center(3))
if(nbas(species(i))==1)goto 10
!compute electronic contribution to dipole - sp terms
elecx=elecx-hyfsp(species(i))*s3(ifirst(i)+offset1(ifirst(i)+1))
elecy=elecy-hyfsp(species(i))*s3(ifirst(i)+offset1(ifirst(i)+2))
elecz=elecz-hyfsp(species(i))*s3(ifirst(i)+offset1(ifirst(i)+3))


       if(nbas(species(i))>4)then
elecx= elecx - ( s3(offset1(ifirst(i)+7)+ifirst(i)+3) + s3(offset1(ifirst(i)+4)+ifirst(i)+1) & 
+ s3(offset1(ifirst(i)+6)+ifirst(i)+2) - factor * s3(offset1(ifirst(i)+5)+ifirst(i)+1) ) * hyfpd(species(i))
 
elecy= elecy - ( s3(offset1(ifirst(i)+8)+ifirst(i)+3) - s3(offset1(ifirst(i)+4)+ifirst(i)+2) &
+ s3(offset1(ifirst(i)+6)+ifirst(i)+1) - factor * s3(offset1(ifirst(i)+5)+ifirst(i)+2) ) * hyfpd(species(i))

elecz= elecz - ( s3(offset1(ifirst(i)+7)+ifirst(i)+1) + s3(offset1(ifirst(i)+8)+ifirst(i)+2) &
                    + two*factor * s3(offset1(ifirst(i)+5)+ifirst(i)+3) ) * hyfpd(species(i))
       end if

10 continue
end do
if(sparkles)then
do i=1,nsparkle
ptchgx=ptchgx+4.80324D0*spcharge(zeffsp(i))*(xsparkle(i)-center(1))
ptchgy=ptchgy+4.80324D0*spcharge(zeffsp(i))*(ysparkle(i)-center(2))
ptchgz=ptchgz+4.80324D0*spcharge(zeffsp(i))*(zsparkle(i)-center(3))
end do
end if



dptchg=dsqrt(ptchgx * ptchgx + ptchgy * ptchgy + ptchgz * ptchgz)
delec=dsqrt(elecx * elecx + elecy * elecy + elecz * elecz)
dx = ptchgx + elecx
dy = ptchgy + elecy
dz = ptchgz + elecz
dmoment=dsqrt(dx * dx + dy * dy + dz * dz)



write(*,*)''
write(*,*)'TABLE OF POINT CHARGE AND HYBRIDIZATION CONTRIBUTIONS TO DIPOLE MOMENT'
write(*,*)'----------------------------------------------------------------------'
write(*,20)'X','Y','Z'
write(*,21)'Point Charge',ptchgx,ptchgy,ptchgz
write(*,21)'Electronic  ',elecx,elecy,elecz
write(*,21)'Total       ',dx,dy,dz
write(*,*)''
write(*,22)'Dipole Moment = ',dmoment



20 format(A19,A10,A10)
21 format(A12,F10.5,F10.5,F10.5)
22 format(A15,F10.5)
end subroutine dipole


