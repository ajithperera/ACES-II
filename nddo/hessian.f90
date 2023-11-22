subroutine make_hess
use constants
use tables
use scratch_array
use control
implicit double precision (a-h,o-z)
interface
subroutine scfopt(x0,y0,z0,gout,natoms,ff)
double precision,intent(inout) ::ff
double precision, dimension(natoms),intent(inout) ::x0,y0,z0
 double precision, dimension(3*natoms),intent(inout) ::gout
end subroutine scfopt
end interface

double precision,allocatable,dimension(:,:)::hessian,gplus,gminus,hessold
double precision,allocatable,dimension(:)::vec2,amass
!double precision,allocatable,dimension(:)::x0,y0,z0
factor=5140.36636949d0 ! conversion factor from mass weighted hessian to cm(-1)
! to derive: (H*mol/bohr**2 grams) are units of eigenvalues of mass weighted
! hessian.  Conversion is then:
! 627.51 kcal/mol * 1000 cal/kcal * 4.184 J/cal * 1bohr**2/.529177(-10)**2 m 
! * 1000 g/kg
!take square root of all the above and divide by 2*pi*c where c is speed of light in cm/s
delta=1D-3
i3n=3*numat

print*,'-----------------------------------------------'
print*,'|                                             |'
print*,'|           VIBRATIONAL FREQUENCY             |'
print*,'-----------------------------------------------'
print*,''
allocate(hessian(i3n,i3n))
allocate(gplus(i3n,1))
allocate(gminus(i3n,1))
optimize=.false. ! set this to fool program into not printing xyz all the 
!time.  see subroutine scfopt.f90 
!allocate(x0(numat))
!allocate(y0(numat))
!allocate(z0(numat))
!x0=x
!y0=y
!x0=z
call cpusec(time1)
write(*,*)'NORMAL MODES WILL BE EVALUATED AT CURRENT GEOMETRY'
write(*,*)'IRRESPECTIVE OF MAGNITUDE OF FIRST DERIVATIVES...'
hessian=zero
icol=0
allocate(amass(i3n))
! evaluate initial gradient

SAVE_TREE=.true.
print*,''
write(*,*)'COMPUTING SECOND DERIVATIVE MATRIX...'
do i=1,numat
print*,i
icol=icol+1
! compute d2E/[dx(i) dq(j)] where q=x,y,z
xold=x(i)
x(i)=x(i)+delta
call scfopt(x,y,z,gplus,numat,ff)
x(i)=x(i)-2.0d0*delta
call scfopt(x,y,z,gminus,numat,ff)
gplus=(gplus-gminus)/two/delta
x(i)=xold
hessian(1:i3n,icol:icol)=gplus(1:i3n,1:1)
amass(icol)=atmass(zeff(i))

icol=icol+1
! compute d2E/[dy(i) dq(j)] where q=x,y,z
yold=y(i)
y(i)=y(i)+delta
call scfopt(x,y,z,gplus,numat,ff)
y(i)=y(i)-2.0d0*delta
call scfopt(x,y,z,gminus,numat,ff)
gplus=(gplus-gminus)/two/delta
y(i)=yold
hessian(1:i3n,icol:icol)=gplus(1:i3n,1:1)
amass(icol)=atmass(zeff(i))


icol=icol+1
! compute d2E/[dz(i) dq(j)] where q=x,y,z
zold=z(i)
z(i)=z(i)+delta
call scfopt(x,y,z,gplus,numat,ff)
z(i)=z(i)-2.0d0*delta
call scfopt(x,y,z,gminus,numat,ff)
gplus=(gplus-gminus)/two/delta
z(i)=zold
hessian(1:i3n,icol:icol)=gplus(1:i3n,1:1)
amass(icol)=atmass(zeff(i))
end do


hessian=hessian*.529177d0*.529177d0/627.51d0 ! convert to atomic units
! symmetrize hessian
hessian=hessian+transpose(hessian)
hessian=hessian/two
! mass weight hessian
amass=1.0d0/dsqrt(amass)
do i=1,i3n
do j=1,i3n
hessian(j,i)=hessian(j,i)*amass(j)*amass(i)
end do
end do
allocate(hessold(i3n,i3n))
hessold=hessian

!diagonalize hessian using tred2/tql2
if(allocated(scrmat))deallocate(scrmat)
if(allocated(scrvec))deallocate(scrvec)
allocate(scrmat(i3n,i3n))
allocate(scrvec(i3n))
allocate(vec2(i3n))

write(*,*)'DIAGONALIZING HESSIAN...'
!call tred2(i3n,i3n,hessian,scrvec,vec2,scrmat)
!call tql2(i3n,i3n,scrvec,vec2,scrmat,iout) ! eigenvals are in scrvec
call eig(hessian, scrmat, 0, i3n, 1)
do i = 1, i3n
call dcopy(1, hessian(i,i), 1, scrvec, 0) 
enddo
print*,'eigenvectors of mass weighted hessian'
!call eigsrt(scrvec,scrmat,i3n,i3n)
!call matprt(scrmat,i3n,i3n,i3n,i3n)

do i=1,i3n
print*,'vector',i
print*,'---------------'
do j=1,i3n
write(*,80)j,scrmat(j,i)
end do
end do
80 format(i2,f10.5)
vec2=zero
do k=1,i3n
total=zero
do j=1,i3n
total=total+scrmat(j,1)*hessold(j,k)
end do
vec2(k)=total
end do
total=zero
do i=1,i3n
total=total+vec2(i)*scrmat(i,1)
end do
print*,'eigenvalue',dsqrt(abs(total))*factor

do i=1,i3n
scrvec(i)=dsqrt(abs(scrvec(i)))*factor*dsign(one,scrvec(i))
end do
!scrvec=abs(scrvec)
!scrvec=dsqrt(scrvec)*factor

write(*,*)'MODE      ','FREQUENCY (cm-1)'
write(*,*)'----      ','----------------'
do i=1,i3n
write(*,90)i,'      ',scrvec(i)
end do
90 format(i4,f20.4)

call cpusec(time2)
write(*,*)'NORMAL MODE ANALYSIS REQUIRED',time2-time1,' seconds'



end subroutine make_hess
