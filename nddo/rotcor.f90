subroutine rotcor(coor,phi,theta,xout,yout,zout,matrix)
use constants
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(inout)::coor,matrix
double precision,intent(in)::phi,theta
double precision,dimension(:),intent(inout)::xout,yout,zout
double precision,dimension(3)::temp
double precision,dimension(3,3)::rotmat

! compute rotation matrix
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

matrix=rotmat
rotmat=transpose(rotmat)


do i=1,3
!call xgemm( 'N', 'N', 3, 1, 3, one,rotmat &
! , 3,coor(:,i) ,3,zero,temp,3 )
!xout(i)=temp(1)
!yout(i)=temp(2)
!zout(i)=temp(3)
xout(i)=rotmat(1,1)*coor(1,i) + rotmat(1,2)*coor(2,i)+rotmat(1,3)*coor(3,i)
yout(i)=rotmat(2,1)*coor(1,i) + rotmat(2,2)*coor(2,i)+rotmat(2,3)*coor(3,i)
zout(i)=rotmat(3,1)*coor(1,i) + rotmat(3,2)*coor(2,i)+rotmat(3,3)*coor(3,i)
end do

end subroutine rotcor



