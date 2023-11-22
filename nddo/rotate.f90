subroutine rotate(twoemol,theta,phi,twoelcl,overlcl,overmol)
use constants
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(inout) ::twoemol,overmol
double precision,intent(in)::theta,phi
double precision,dimension(:),intent(in)::twoelcl,overlcl
double precision,dimension(3,3)::rotmat
!**********************************************
!compute rotation matrices. 
!the px,py,and pz orbitals transform exactly as the 
!cartesian x,y,and z axes. (no sh*t)
! this matrix corresponds to (RyRz)^T(pi pi' o) => (x y z)
! so, the first row is x (molecular coordinates) expanded in terms of local coordinates
! i.e. x(molecular) = rotmat(1,1)*pi + rotmat(1,2)*pi' + rotmat(1,3)*o
! therefore, to derive it, take the standard textbook rotation matrix product 
!Ry*Rz and transpose it.  The matrix is as follows:
!
!   [  cos phi * cos theta      -sin phi        cos phi * sin theta   ]     >x across row
!   |                                                                 |
!   |                                                                 |
!   |  sin phi * cos theta       cos phi        sin phi * sin theta   |     >y across row
!   |                                                                 |
!   |                                                                 |
!   [  -sin theta                  0                cos theta         ]     >z across row

cosphi=dcos(phi)
sinphi=dsin(phi)
costh=dcos(theta)
sinth=dsin(theta)
rotmat(1,1)=cosphi*costh
rotmat(1,2)=-sinphi
rotmat(1,3)=cosphi*sinth
rotmat(2,1)=sinphi*costh
rotmat(2,2)=cosphi
rotmat(2,3)=sinphi*sinth
rotmat(3,1)=-sinth
rotmat(3,2)=zero
rotmat(3,3)=costh
fac=180.0d0/pi
!print*,'theta phi',theta*fac,phi*fac
!call matprt(rotmat,3,3,3,3)
!***************************************************

!****************************************************





icol=1
! 1=(ss|ss)
twoemol(1,icol)=twoelcl(1)
!print*,'(ss|ss) =',twoemol(1,1)


! (ss|ps) 
do i=1,3
icol=icol+1
twoemol(1,icol)=twoelcl(4)*rotmat(i,3)
end do
!print*,'(ss|ps) order xyz',twoemol(1,2),twoemol(1,3),twoemol(1,4)



!! (ss|pp)
!do i=1,3
!icol=icol+1
!twoemol(1,icol)=( rotmat(i,1)*rotmat(i,1) +  rotmat(i,2)*rotmat(i,2) )*twoelcl(2)  + rotmat(i,3)*rotmat(i,3) * twoelcl(3)
!end do

!print*,'(ss|pp) order xyz',twoemol(1,5),twoemol(1,6),twoemol(1,7)

! (ss|p p)
do i=1,3
do j=i,3
icol=icol+1
twoemol(1,icol)=(rotmat(i,1)*rotmat(j,1) + rotmat(i,2)*rotmat(j,2)) * twoelcl(2) + &
rotmat(i,3)*rotmat(j,3)*twoelcl(3)
end do 
end do
!print*,'(ss|pp'') order xy xz yz',twoemol(1,8),twoemol(1,9),twoemol(1,10)


! now do the second row <xs|...

 irow=2

! (ps|ss) 
do i=1,3
twoemol(irow,1)=twoelcl(7)*rotmat(i,3)
irow=irow+1
end do
!print*,'(ps|ss) order xyz',twoemol(2,1),twoemol(3,1),twoemol(4,1)

! (ps | ps)
irow=2
icol=1
do i=1,3
do j=1,3
icol=icol+1
twoemol(irow,icol)=(rotmat(i,1)*rotmat(j,1)+rotmat(i,2)*rotmat(j,2))*twoelcl(12)  + &
  rotmat(i,3)*rotmat(j,3)*twoelcl(13)
end do 
irow=irow+1
icol=1
end do

!print*,'[ps|ps]',twoemol(2,2),twoemol(2,3),twoemol(2,4)
!print*,'[ps|ps]',twoemol(3,2),twoemol(3,3),twoemol(3,4)
!print*,'[ps|ps]',twoemol(4,2),twoemol(4,3),twoemol(4,4)

! (ps | pp)
irow=1
icol=5
do i=1,3
irow=irow+1
do j=1,3
do k=j,3
twoemol(irow,icol)=rotmat(i,1)*(rotmat(j,1)*rotmat(k,3) + rotmat(j,3)*rotmat(k,1))*twoelcl(14) &
+ rotmat(i,2)*(rotmat(j,2)*rotmat(k,3) + rotmat(j,3)*rotmat(k,2))*twoelcl(14) &
+ rotmat(i,3)*(rotmat(j,1)*rotmat(k,1)*twoelcl(8) + rotmat(j,2)*rotmat(k,2)*twoelcl(8) + rotmat(j,3)*rotmat(k,3)* &
twoelcl(10))
icol=icol+1
end do
end do
icol=5
end do

!print*,'ps|pp',twoemol(2,5),twoemol(2,6),twoemol(2,7),twoemol(2,8),twoemol(2,9),twoemol(2,10)
!print*,'ps|pp',twoemol(3,5),twoemol(3,6),twoemol(3,7),twoemol(3,8),twoemol(3,9),twoemol(3,10)
!print*,'ps|pp',twoemol(4,5),twoemol(4,6),twoemol(4,7),twoemol(4,8),twoemol(4,9),twoemol(4,10)

irow=4
! (pp|ss)
do i=1,3
do j=i,3
irow=irow+1
twoemol(irow,1)=( rotmat(i,1)*rotmat(j,1) +  rotmat(i,2)*rotmat(j,2) )*twoelcl(5)  + rotmat(i,3)* &
rotmat(j,3)* twoelcl(6)
end do
end do
!print*,'pp|ss',twoemol(5,1)
!print*,'pp|ss',twoemol(6,1)
!print*,'pp|ss',twoemol(7,1)
!print*,'pp|ss',twoemol(8,1)
!print*,'pp|ss',twoemol(9,1)
!print*,'pp|ss',twoemol(10,1)



! (pp | ps)
irow=4
icol=1
do i=1,3
do j=i,3
irow=irow+1
do k=1,3
icol=icol+1
twoemol(irow,icol)=rotmat(k,1)*(rotmat(i,1)*rotmat(j,3) + rotmat(i,3)*rotmat(j,1))*twoelcl(15)  + &
rotmat(k,2)*(rotmat(i,2)*rotmat(j,3) + rotmat(i,3)*rotmat(j,2))*twoelcl(15) + &
rotmat(k,3)*(rotmat(i,1)*rotmat(j,1) + rotmat(i,2)*rotmat(j,2))*twoelcl(9) + &
rotmat(k,3)*(rotmat(i,3)*rotmat(j,3))*twoelcl(11)
end do
icol=1
end do
end do

!print*,'pp|ps',twoemol(5,2),twoemol(5,3),twoemol(5,4)
!print*,'pp|ps',twoemol(6,2),twoemol(6,3),twoemol(6,4)
!print*,'pp|ps',twoemol(7,2),twoemol(7,3),twoemol(7,4)
!print*,'pp|ps',twoemol(8,2),twoemol(8,3),twoemol(8,4)
!print*,'pp|ps',twoemol(9,2),twoemol(9,3),twoemol(9,4)
!print*,'pp|ps',twoemol(10,2),twoemol(10,3),twoemol(10,4)

irow=4
do i=1,3
do j=i,3
irow=irow+1
icol=4
do k=1,3
do l=k,3
icol=icol+1
term1 = rotmat(k,1)*rotmat(l,1)*twoelcl(17) + rotmat(k,2)*rotmat(l,2) * twoelcl(18) +  rotmat(k,3)*rotmat(l,3) * twoelcl(19) 
term2 = rotmat(k,1)*rotmat(l,1)*twoelcl(18) + rotmat(k,2)*rotmat(l,2) * twoelcl(17) +  rotmat(k,3)*rotmat(l,3) * twoelcl(19)
term1=rotmat(i,1)*rotmat(j,1)*term1 
term2=rotmat(i,2)*rotmat(j,2)*term2 

term3 = (rotmat(k,1)*rotmat(l,1)+ rotmat(k,2)*rotmat(l,2)) * twoelcl(20) +  rotmat(k,3)*rotmat(l,3) * twoelcl(21) 
term3=term3*rotmat(i,3)*rotmat(j,3)

term4=rotmat(i,1)*rotmat(j,2)*(rotmat(k,1)*rotmat(l,2)+rotmat(k,2)*rotmat(l,1))*twoelcl(22)
term5=rotmat(i,2)*rotmat(j,1)*(rotmat(k,2)*rotmat(l,1)+rotmat(k,1)*rotmat(l,2))*twoelcl(22)

term6=rotmat(i,1)*rotmat(j,3)*(rotmat(k,1)*rotmat(l,3)+rotmat(k,3)*rotmat(l,1))*twoelcl(16)
term7=rotmat(i,3)*rotmat(j,1)*(rotmat(k,3)*rotmat(l,1)+rotmat(k,1)*rotmat(l,3))*twoelcl(16)

term8=rotmat(i,2)*rotmat(j,3)*(rotmat(k,2)*rotmat(l,3)+rotmat(k,3)*rotmat(l,2))*twoelcl(16)
term9=rotmat(i,3)*rotmat(j,2)*(rotmat(k,3)*rotmat(l,2)+rotmat(k,2)*rotmat(l,3))*twoelcl(16)

twoemol(irow,icol)=term1+term2+term3+term4+term5+term6+term7+term8+term9

end do
end do
end do
end do


!print*,'pp|pp',twoemol(5,5),twoemol(5,6),twoemol(5,7),twoemol(5,8),twoemol(5,9),twoemol(5,10)
!print*,'pp|pp',twoemol(6,5),twoemol(6,6),twoemol(6,7),twoemol(6,8),twoemol(6,9),twoemol(6,10)
!print*,'pp|pp',twoemol(7,5),twoemol(7,6),twoemol(7,7),twoemol(7,8),twoemol(7,9),twoemol(7,10)
!print*,'pp|pp',twoemol(8,5),twoemol(8,6),twoemol(8,7),twoemol(8,8),twoemol(8,9),twoemol(8,10)
!print*,'pp|pp',twoemol(9,5),twoemol(9,6),twoemol(9,7),twoemol(9,8),twoemol(9,9),twoemol(9,10)
!print*,'pp|pp',twoemol(10,5),twoemol(10,6),twoemol(10,7),twoemol(10,8),twoemol(10,9),twoemol(10,10)





!  compute overlaps in molecular coordinates
!(s|s)
overmol(1,1)=overlcl(1)
!(s|p)
overmol(1,2)=rotmat(1,3)*overlcl(2)
overmol(1,3)=rotmat(2,3)*overlcl(2)
overmol(1,4)=rotmat(3,3)*overlcl(2)

!(p|s)
overmol(2,1)=rotmat(1,3)*overlcl(3)
overmol(3,1)=rotmat(2,3)*overlcl(3)
overmol(4,1)=rotmat(3,3)*overlcl(3)

!(p|p) 
!xx
overmol(2,2)=(rotmat(1,1)*rotmat(1,1)+rotmat(1,2)*rotmat(1,2))*overlcl(5) + (rotmat(1,3)*rotmat(1,3)) * overlcl(4)
!xy
overmol(2,3)=(rotmat(1,1)*rotmat(2,1)+rotmat(1,2)*rotmat(2,2))*overlcl(5) + (rotmat(1,3)*rotmat(2,3) )* overlcl(4)
!xz
overmol(2,4)=(rotmat(1,1)*rotmat(3,1)+rotmat(1,2)*rotmat(3,2))*overlcl(5) + (rotmat(1,3)*rotmat(3,3) )* overlcl(4)

!yx
overmol(3,2)=(rotmat(2,1)*rotmat(1,1)+rotmat(2,2)*rotmat(1,2))*overlcl(5) + (rotmat(2,3)*rotmat(1,3) )* overlcl(4)
!yy
overmol(3,3)=(rotmat(2,1)*rotmat(2,1)+rotmat(2,2)*rotmat(2,2))*overlcl(5) + (rotmat(2,3)*rotmat(2,3) )* overlcl(4)
!yz
overmol(3,4)=(rotmat(2,1)*rotmat(3,1)+rotmat(2,2)*rotmat(3,2))*overlcl(5) + (rotmat(2,3)*rotmat(3,3) )* overlcl(4)

!zx
overmol(4,2)=(rotmat(3,1)*rotmat(1,1)+rotmat(3,2)*rotmat(1,2))*overlcl(5) + (rotmat(3,3)*rotmat(1,3) )* overlcl(4)
!zy
overmol(4,3)=(rotmat(3,1)*rotmat(2,1)+rotmat(3,2)*rotmat(2,2))*overlcl(5) + (rotmat(3,3)*rotmat(2,3) )* overlcl(4)
!zz
overmol(4,4)=(rotmat(3,1)*rotmat(3,1)+rotmat(3,2)*rotmat(3,2))*overlcl(5) + (rotmat(3,3)*rotmat(3,3) )* overlcl(4)









 end subroutine rotate
