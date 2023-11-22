subroutine geometry
use constants
use tables
implicit double precision (a-h,o-z)
double precision,dimension(4,3)::coor
double precision,dimension(4)::xout,yout,zout
double precision,dimension(3,3)::matrix
double precision,dimension(:),allocatable::oldx,oldy,oldz

interface
subroutine rotcor(coor,phi,theta,xout,yout,zout,matrix)
double precision,dimension(:,:),intent(inout)::coor,matrix
double precision,intent(in)::phi,theta
double precision,dimension(:),intent(inout)::xout,yout,zout
end subroutine rotcor
end interface

if(allocated(x))deallocate(x)
if(allocated(y))deallocate(y)
if(allocated(z))deallocate(z)
allocate(x(ninput))
allocate(y(ninput))
allocate(z(ninput))
x=zero
y=zero
z=zero
  
!place first atom at origin and place second atom down x axis
x(2)=q(2,1)  !whew, that was hard. on to the third atom
! convert all angles to radians
!factor=pi/180.0d0
!do i=1,numat
!q(i,2)=factor*q(i,2)
!end do
! use law of cosines to find angle of atom 3 - 1 - x axis. (compute phi in polar coords)
r1=q(2,1)
r2=q(3,1) ! r1 and r2 define lengths of side of triangle, now we need third side
r1sq = r1 * r1
r2sq = r2 * r2
r3sq = r1sq + r2sq - two * r1 * r2 * dcos(q(3,2))
r3=dsqrt(r3sq)
! now use law of cosines again to compute angle we seek

angle = r2sq - r1sq - r3sq 
angle = - angle / (two * r1 * r3)
x(3)=r3 * angle
angle=acos(angle)
y(3)= r3 * dsin(angle)


! now add the fourth atom in a loop
! so, we are trying to insert 4 - 3 - 2 - 1 which is labeling of atoms across row in ZMAT
! lets put atoms 3 and 2 on z axis (3 at origin)

! loop over remaining atoms

do i=4,ninput

i1=ref(i,1)
i2=ref(i,2)
i3=ref(i,3)
xold=x(i2)
yold=y(i2)
zold=z(i2)

x(i1)=x(i1)-xold
y(i1)=y(i1)-yold
z(i1)=z(i1)-zold
x(i3)=x(i3)-xold
y(i3)=y(i3)-yold
z(i3)=z(i3)-zold
x(i2)=zero
y(i2)=zero
z(i2)=zero

coor(1,1)=x(i3)
coor(2,1)=y(i3)
coor(3,1)=z(i3)
coor(1,2)=x(i2)
coor(2,2)=y(i2)
coor(3,2)=z(i2)
coor(1,3)=x(i1)
coor(2,3)=y(i1)
coor(3,3)=z(i1)


! now, rotate so that atoms 2 and 3 are on z axis ! this code is taken directly from sub two_electron

rsq=x(i1)*x(i1) + y(i1)*y(i1) + z(i1)*z(i1)
   r=dsqrt(rsq)
x1=x(i1)
y1=y(i1)
z1=z(i1)

 costh=z1/r

if(abs(costh).gt.one)costh=dsign(one,costh)
 theta=dacos(costh)
 sinth=dsin(theta)
 if(abs(sinth)<1D-10)sinth=zero
 if(abs(sinth)==zero)then  ! there is a problem here with machine precision. 1d-20 effectively turns
! this check off
     phi=zero
     if(z1.gt.zero)then
     theta=zero
     else
     theta=pi
     end if
! if not on z axis then...
else
         cosphi=x1/(r*sinth)
         if(abs(cosphi).gt.one)cosphi=dsign(one,cosphi)
         phi=dacos(cosphi)
         cosphi=dcos(phi)
!         print*,'so phi is',phi*180.0d0/pi

! if on y axis or in plane of y axis then
!        rdiff=abs(phi)-(pi/two)
!    if(abs(phi).eq.(pi/two))then
!        if(abs(rdiff).lt.1.0d-6)then !BBBBBBBBBBBBBb
         if(abs(cosphi)<1D-10)cosphi=zero
         if(cosphi==zero)then
              if(y1.gt.zero)then
              phi=pi/two
              else
              phi=three*pi/two
              end if
!  elseif(abs(cosphi).lt.0.99999999999999999999d0)then
         elseif(y1==zero.and.z1==zero)then
               if(x1<zero)then
                  phi=pi
                  else
                  phi=zero
                  end if
         else
!            print*,'call quad',x1,y1,z1,phi*180.0d0/pi
         call quadrant(x1,y1,z1,phi,sinth)
!         print*,'after',phi*180.0d0/pi
          end if


end if

xout=zero
yout=zero
zout=zero

call rotcor(coor,phi,theta,xout,yout,zout,matrix)


r1=zout(3)
r2=q(i,1)
r1sq=r1*r1
r2sq=r2*r2
r3sq=r1sq + r2sq - two *r1 * r2 * dcos(q(i,2))
r3=dsqrt(r3sq)

angle = r2sq - r1sq - r3sq 
angle = - angle / (two * r1 * r3)
zout(4)=r3 * angle
angle=acos(angle)
xout(4)= r3 * dsin(angle)

! so now we have atom 4 in at proper r and theta but now we need to rotate to proper dihedral orientation
! first, compute phi for atom 1
rsq=xout(1)*xout(1) + yout(1)*yout(1) + zout(1)*zout(1)
   r=dsqrt(rsq)
x1=xout(1)
y1=yout(1)
z1=zout(1)

 costh=z1/r

if(abs(costh).gt.one)costh=dsign(one,costh)
 theta=dacos(costh)
 sinth=dsin(theta)
 if(abs(sinth)<1D-10)sinth=zero
 if(abs(sinth)==zero)then  ! there is a problem here with machine precision. 1d-20 effectively turns
! this check off
     phi=zero
     if(z1.gt.zero)then
     theta=zero
     else
     theta=pi
     end if
! if not on z axis then...
else
         cosphi=x1/(r*sinth)
         if(abs(cosphi).gt.one)cosphi=dsign(one,cosphi)
         phi=dacos(cosphi)
         cosphi=dcos(phi)
!         print*,'so phi is',phi*180.0d0/pi

! if on y axis or in plane of y axis then
!        rdiff=abs(phi)-(pi/two)
!    if(abs(phi).eq.(pi/two))then
!        if(abs(rdiff).lt.1.0d-6)then !BBBBBBBBBBBBBb
         if(abs(cosphi)<1D-10)cosphi=zero
         if(cosphi==zero)then
              if(y1.gt.zero)then
              phi=pi/two
              else
              phi=three*pi/two
              end if
!  elseif(abs(cosphi).lt.0.99999999999999999999d0)then
         elseif(y1==zero.and.z1==zero)then
               if(x1<zero)then
                  phi=pi
                  else
                  phi=zero
                  end if
         else
!            print*,'call quad',x1,y1,z1,phi*180.0d0/pi
         call quadrant(x1,y1,z1,phi,sinth)
!         print*,'after',phi*180.0d0/pi
          end if


end if


!angle=phi+q(i,3)*factor
angle=phi+q(i,3)
! now, about z axis of new atom
xx=xout(4)
yy=yout(4)
xout(4)=dcos(angle)*xx-dsin(angle)*yy
yout(4)=dsin(angle)*xx+dcos(angle)*yy
do ii=1,4
coor(ii,1)=xout(ii)
coor(ii,2)=yout(ii)
coor(ii,3)=zout(ii)

end do

! put it back in molecular frame

do k=4,4
xx=xout(4)
yy=yout(4)
zz=zout(4)
xout(k)=matrix(1,1)*xx + matrix(1,2)*yy+matrix(1,3)*zz
yout(k)=matrix(2,1)*xx + matrix(2,2)*yy+matrix(2,3)*zz
zout(k)=matrix(3,1)*xx + matrix(3,2)*yy+matrix(3,3)*zz
xout(k)=xout(k)+xold
yout(k)=yout(k)+yold
zout(k)=zout(k)+zold
end do

x(i1)=x(i1)+xold
y(i1)=y(i1)+yold
z(i1)=z(i1)+zold
x(i3)=x(i3)+xold
y(i3)=y(i3)+yold
z(i3)=z(i3)+zold
x(i2)=xold
y(i2)=yold
z(i2)=zold
x(i)=xout(4)
y(i)=yout(4)
z(i)=zout(4)














end do


! put angles back in degrees for geometry optimization
!factor=180.0d0/pi
!do i=1,numat
!q(i,2)=factor*q(i,2)
!end do






end subroutine geometry
  








