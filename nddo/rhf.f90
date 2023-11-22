subroutine rhf
use constants
use indices
use scratch_array
use tables
use control
implicit double precision(a-h,o-z)

interface
subroutine mkvector(matrix,vector,idim1)
use indices
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(in)::matrix
double precision,dimension(:),intent(inout)::vector
integer,intent(in)::idim1
end subroutine mkvector

subroutine guess(density)
double precision,dimension(:),intent(inout)::density
end subroutine guess

subroutine f1(fock,density)
double precision,dimension(:),intent(inout)::fock
double precision,dimension(:),intent(in)::density
end subroutine f1

subroutine trace(S,N,out)
double precision,dimension(:,:),intent(in)::S
integer::N
double precision,intent(out)::out
end subroutine trace

subroutine f2(fock,density)
double precision,dimension(:),intent(inout)::fock
double precision,dimension(:),intent(in)::density
end subroutine f2

subroutine square(matrix,vector,idim1)
use indices
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(inout)::matrix
double precision,dimension(:),intent(in)::vector
integer,intent(in)::idim1
end subroutine square

              SUBROUTINE mmult(r,t,z,NBASIS,icol)
    DOUBLE PRECISION,dimension(:,:),intent(inout)::r,t,z
     double precision:: total
     integer::nbasis,icol
end subroutine mmult

end interface


!pdouble precision,dimension(:,:),allocatable::fock,density,vecs,test
!pdouble precision,dimension(:),allocatable::vec2,work
double precision,dimension(:,:),allocatable::vecs,test
double precision,dimension(:),allocatable::fock,density,vec2,work,scrvec3
integer,dimension(:),allocatable::iwork,ifail
!double precision,dimension(4,4)::plocal

!stdoutif(.not.save_tree)then
!stdoutprint*,'-----------------------------------------------'
!stdoutprint*,'|                                             |'
!stdoutprint*,'|    HARTREE FOCK SELF CONSISTENT FIELD       |'
!stdoutprint*,'-----------------------------------------------'
!stdoutprint*,''
!stdoutend if
!write(*,*)'@RHF-NDDO'
ihigh=num_basis-(nelectrons/2)+1

!pallocate(density(num_basis,num_basis))
allocate(density(ndim1))

!THIS WILL HAVE TO BE FIXED LATER ONCE WE DECIDE HOW TO HANDLE FINITE DIFFERENCE GRADIENTS
if(TRIAL)then
density=s3
else
call guess(density)
end if

! put this in the above if construct later




! build the diagonal fock matrix elements
if(.not.save_tree)then
!stdoutprint*,'Restricted Hartree Fock Field to be determined...'
end if
if(densityin)then
write(*,*)'Reading initial density from disk'
open(unit=1,file='DENSITY')
read(1,*)density
close(1)
end if



! do the scf iterations
allocate(scrvec(num_basis))
allocate(scrmat(num_basis,num_basis))
allocate(vecs(num_basis,num_basis))
allocate(vec2(num_basis))
allocate(scrvec3(ndim1))
allocate(fock(ndim1))
icol=nelectrons/2
if(lapack)then 
nroot=icol+1 !by default compute occupied subspace and 1 virtual unless user wants more (or less)! to be implemented -> note this can be done simply by letting user define nroot from the input deck.  
if(.not.save_tree)then
!stdoutwrite(*,*)'SCF will use incomplete diagonalizations'
!stdoutwrite(*,*)'Program will determine',nroot,'lowest energy orbitals of which',icol,'are occupied.'
end if
!lapack specific stuff
EVTOL = 2.0D0*DLAMCH('S')
lwork=8*num_basis
jwork=5*num_basis
else
if(.not.save_tree)then
!stdoutwrite(*,*)'SCF will determine ALL eigenvectors of occupied/virtual spaces.'
!stdoutwrite(*,*)'Note that the virtual space in not needed in SCF level calculations.'
!stdoutwrite(*,*)'If matrices are large it is useful to set Lapack keyword in input'
!stdoutwrite(*,*)'file which restricts diagonalization to the occupied subspace'
end if
end if


!allocate(test(num_basis,icol))

energy=zero
iter=0
if(.not.save_tree)then
!stdoutwrite(*,*)'Beginning SCF iterations. Convergence tolerance = ',scftol
!write(*,*)'CYCLE          ENERGY                    DELTA E                   DELTA P'
end if
! scf iteration
1 iter=iter+1
if(iter>40)then
print*,'NDDO SCF DID NOT CONVERGE AFTER',iter-1,' CYCLES'
!stop
goto 20
end if
call cpusec(time1)
call f1(fock,density)
call f2(fock,density)


old=energy
scrvec3=zero
!pscrmat=density*(fock+H)
scrvec3=density*(fock+H)

call square(scrmat,scrvec3,num_basis)
energy=half*sum(scrmat)
delta=abs(energy-old)
if(.not.save_tree)then
!write(*,40)iter,energy,delta,pdelta
end if
40 format(i3,f25.10,f25.10,f25.10)
if(pdelta.lt.scftol.and.delta.lt.scftol)goto 20

scrvec3=density


call cpusec(time1)
!do i=1,ndim1
!print*,i,fock(i),density(i)/two
!end do

!if(.not.Lapack)then ! if not subspace diagonalization or user hates lapack, do the full matrix with tred/tqli combo

!call square(scrmat,fock,num_basis)
!call cpusec(time1)

!call tred2(num_basis,num_basis,scrmat,scrvec,vec2,vecs)
!call tql2(num_basis,num_basis,scrvec,vec2,vecs,iout)
!call cpusec(time2)
!else

allocate(work(lwork))
allocate(iwork(jwork))
allocate(ifail(num_basis))
call cpusec(time1)
call square(scrmat,fock,num_basis)
!call dsyevx('V','I','U',num_basis,scrmat,num_basis,0.,0.,1,NROOT,EVTOL, &
!nevs,scrvec,vecs,num_basis,work,lwork,iwork,ifail,info)
print *, 'ERROR: no dsyevx'
call aces_exit(1)
call cpusec(time2)
deallocate(work)
deallocate(iwork)
deallocate(ifail)
do i=nroot+1,num_basis
scrvec(i)=1D20
end do
!end if

call eigsrt(scrvec,vecs,num_basis,num_basis)
call cpusec(time2)
44 continue
call xgemm( 'N', 'T', num_basis, num_basis, icol, two,vecs(1:num_basis,ihigh:num_basis) &
 , num_basis,vecs(1:num_basis,ihigh:num_basis) , num_basis,zero,scrmat,num_basis )


call mkvector(scrmat,density,num_basis)



pdelta=zero
do i=1,ndim1
pdelta=pdelta+(density(i)-scrvec3(i))**2
end do
!call matprt(density/2,num_basis,num_basis,num_basis,num_basis)


! compute energy
go to 1

20 etotal=energy+enuc
if(.not.save_tree)then
!stdout write(*,*)'***********************************'
!stdoutwrite(*,*)'            SCF RESULTS               '
!stdoutprint*,'Nuclear repulsion energy (eV) = ',enuc
!stdoutprint*,'Electronic energy (eV) = ',energy
!stdoutprint*,'SCF Total energy (eV) = ',etotal
!stdoutprint*,'homo energy (eV) = ',scrvec(ihigh)
!stdoutprint*,'lumo energy (eV) = ',scrvec(ihigh-1)
call square(scrmat,density,num_basis)
call trace(scrmat,num_basis,out1)
!stdoutprint*,'Total # electrons at convergence (Tr P)= ',nint(out1)
!stdoutwrite(*,*)'***********************************'
end if


if(keep)then
    if(allocated(s3))then
    deallocate(s3)
    end if
  allocate(s3(ndim1))
  s3=density
end if

if(densityout)then ! if true then write density to file

open(unit=1,file='DENSITY')
write(1,*)density
close(1)
write(*,*)'RHF density has been written to disk'
end if





deallocate(fock)
deallocate(density)
deallocate(vecs)
deallocate(vec2)
deallocate(scrvec)
deallocate(scrmat)
deallocate(scrvec3)


!10 format(i4,f15.8,f15.8,f15.8)

end subroutine rhf

     SUBROUTINE mmult(r,t,z,NBASIS,icol)
implicit double precision(a-h,o-z)
      DOUBLE PRECISION,dimension(:,:),intent(inout)::r,t,z
     double precision:: total
     integer::nbasis,icol
      do  i=1,NBASIS
         do  j=1,NBASIS
            total=0.0D0
            do  k=1,icol
               total=total+r(i,k)*t(k,j)
   end do
               z(i,j)=total
 end do
 end do
         
         end subroutine mmult



