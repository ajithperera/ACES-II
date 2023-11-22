subroutine uhf
!
!  this is a modification of the rhf subroutine
!
!
use constants
use indices
use scratch_array
use tables
use control
implicit double precision(a-h,o-z)

interface
subroutine guess(density)
double precision,dimension(:),intent(inout)::density
end subroutine guess

subroutine f1uhf(fock,density,abdens)
double precision,dimension(:),intent(inout)::fock
double precision,dimension(:),intent(in)::density,abdens
end subroutine f1uhf

subroutine trace(S,N,out)
double precision,dimension(:,:),intent(in)::S
integer::N
double precision,intent(out)::out
end subroutine trace

subroutine f2uhf(fock,density,abdens)
double precision,dimension(:),intent(inout)::fock
double precision,dimension(:),intent(in)::density,abdens
end subroutine f2uhf

              SUBROUTINE mmult(r,t,z,NBASIS,icol)
    DOUBLE PRECISION,dimension(:,:),intent(inout)::r,t,z
     double precision:: total
     integer::nbasis,icol
end subroutine mmult

subroutine square(matrix,vector,idim1)
use indices
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(inout)::matrix
double precision,dimension(:),intent(in)::vector
integer,intent(in)::idim1
end subroutine square

subroutine mkvector(matrix,vector,idim1)
use indices
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(in)::matrix
double precision,dimension(:),intent(inout)::vector
integer,intent(in)::idim1
end subroutine mkvector

end interface


double precision,dimension(:,:),allocatable::vecs,scrmat2
double precision,dimension(:),allocatable::vec2,work,focka,density,scrvec3
! extra stuff needed for uhf
double precision,dimension(:),allocatable::fockb
integer,dimension(:),allocatable::iwork,ifail

!stdoutprint*,'-----------------------------------------------'
!stdoutprint*,'|                                             |'
!stdoutprint*,'|    HARTREE FOCK SELF CONSISTENT FIELD       |'
!stdoutprint*,'-----------------------------------------------'
!stdoutprint*,'Unrestricted Hartree Fock Field to be determined...'
!ihigh=num_basis-(nelectrons/2)+1
write(*,*)'@UHF-NDDO'

allocate(density(ndim1))
if(TRIAL)then
density=s3
else
call guess(density)
end if



! get some memory
allocate(focka(ndim1))
allocate(fockb(ndim1))
if(allocated(adens))deallocate(adens)
if(allocated(bdens))deallocate(bdens)
allocate(adens(ndim1))
allocate(bdens(ndim1))
allocate(scrvec(num_basis))
allocate(scrmat(num_basis,num_basis))
allocate(scrmat2(num_basis,num_basis))
allocate(vecs(num_basis,num_basis))
allocate(vec2(num_basis))
allocate(scrvec3(ndim1))
!icol=nelectrons/2

! determine # of paired and unpaired electrons and max summation
! indices for alpha and beta density matrices
!print*,'Spin Multiplicity = ',int(mult)
iunpair=int(mult)-1
nclose=nelectrons-iunpair
numalpha=nclose/2+iunpair
numbeta=nclose/2
indexa=num_basis-numalpha+1
indexb=num_basis-numbeta+1
!print*,'# Alpha electrons =',numalpha
!print*,'# Beta  electrons =',numbeta

if(lapack)then 
nroota=numalpha+1 !by default compute occupied subspace and 1 virtual unless user wants more (or less)! to be implemented -> note this can be done simply by letting user define nroot from the input deck.  
nrootb=numbeta+1
!stdoutwrite(*,*)'SCF will use incomplete diagonalizations'
!stdoutwrite(*,*)'Program will determine',nroota,'lowest energy orbitals of alpha spin'
!stdoutwrite(*,*)'Program will determine',nrootb,'lowest energy orbitals of beta spin'
!lapack specific stuff
EVTOL = 2.0D0*DLAMCH('S')
lwork=8*num_basis
jwork=5*num_basis
else
!stdoutwrite(*,*)'SCF will determine ALL eigenvectors of occupied/virtual spaces.'
!stdoutwrite(*,*)'Note that the virtual space in not needed in SCF level calculations.'
!stdoutwrite(*,*)'If matrices are large it is useful to set Lapack keyword in input'
!stdoutwrite(*,*)'file which restricts diagonalization to the occupied subspace'
end if


adens=density/2
bdens=density/2
! put some arbitrary noise into the alpha and beta density matrices
! i got this from thiel's code.  i have no idea why it is done this
! way.
adens=adens*two*float(numalpha)/(numalpha+numbeta)
bdens=bdens*two*float(numbeta)/(numalpha+numbeta)
if( (numalpha==numbeta) .and. mult.ne.1)then
da=0.98d0
db=two-da
do i=1,num_basis
dc=da
da=db
db=dc
adens(i+offset1(i))=adens(i+offset1(i))*da
bdens(i+offset1(i))=bdens(i+offset1(i))*db
end do
end if

!call matprt(adens,num_basis,num_basis,num_basis,num_basis)
!call matprt(bdens,num_basis,num_basis,num_basis,num_basis)


!stdoutprint*,'Initial guess density is not N representable!'
! do the scf iterations

energy=zero
iter=0
!stdoutwrite(*,*)'Beginning SCF iterations. Convergence tolerance = ',scftol
!write(*,*)'CYCLE          ENERGY                    DELTA E                   DELTA P'
! scf iteration
1 iter=iter+1
if(iter>40)then
print*,'NDDO SCF DID NOT CONVERGE AFTER',iter-1,' CYCLES'
!stop
goto 20
end if
call f1uhf(focka,density,adens)
call f2uhf(focka,density,adens)
call f1uhf(fockb,density,bdens)
call f2uhf(fockb,density,bdens)
old=energy
scrvec3=density*h
scrvec3=scrvec3+adens*focka
scrvec3=scrvec3+bdens*fockb
call square(scrmat,scrvec3,num_basis)

energy=half*sum(scrmat)
delta=abs(energy-old)
!write(*,40)iter,energy,delta,pdelta
40 format(i3,f25.10,f25.10,f25.10)
if(pdelta.lt.scftol.and.delta.lt.scftol)goto 20

scrvec3=density
! diagonalize fock matrix with lapack or tqli/tred
!note that beta must be diagonalized before alpha b/c alpha homo energy is written last
!if(.not.Lapack)then 
!call square(scrmat,fockb,num_basis)
!call tred2(num_basis,num_basis,scrmat,scrvec,vec2,vecs)
!call tql2(num_basis,num_basis,scrvec,vec2,vecs,iout)
!call eigsrt(scrvec,vecs,num_basis,num_basis)
!call xgemm( 'N', 'T', num_basis, num_basis, numbeta, one,vecs(1:num_basis,indexb:num_basis) &
! , num_basis,vecs(1:num_basis,indexb:num_basis) , num_basis,zero,scrmat,num_basis )
!call mkvector(scrmat,bdens,num_basis)
!call square(scrmat,focka,num_basis)
!call tred2(num_basis,num_basis,scrmat,scrvec,vec2,vecs)
!call tql2(num_basis,num_basis,scrvec,vec2,vecs,iout)
!call eigsrt(scrvec,vecs,num_basis,num_basis)
!call xgemm( 'N', 'T', num_basis, num_basis, numalpha, one,vecs(1:num_basis,indexa:num_basis) &
! , num_basis,vecs(1:num_basis,indexa:num_basis) , num_basis,zero,scrmat,num_basis )
!call mkvector(scrmat,adens,num_basis)

!else



allocate(work(lwork))
allocate(iwork(jwork))
allocate(ifail(num_basis))
call square(scrmat,fockb,num_basis)
!call dsyevx('V','I','U',num_basis,scrmat,num_basis,0.,0.,1,nrootb,EVTOL, &
!nevs,scrvec,vecs,num_basis,work,lwork,iwork,ifail,info)
print *, 'ERROR: no dsyevx'
call aces_exit(1)
deallocate(work)
deallocate(iwork)
deallocate(ifail)
do i=nrootb+1,num_basis
scrvec(i)=1D20
end do
call eigsrt(scrvec,vecs,num_basis,num_basis)
call xgemm( 'N', 'T', num_basis, num_basis, numbeta, one,vecs(1:num_basis,indexb:num_basis) &
 , num_basis,vecs(1:num_basis,indexb:num_basis) , num_basis,zero,scrmat,num_basis )
call mkvector(scrmat,bdens,num_basis)
allocate(work(lwork))
allocate(iwork(jwork))
allocate(ifail(num_basis))
call square(scrmat,focka,num_basis)
!call dsyevx('V','I','U',num_basis,scrmat,num_basis,0.,0.,1,nroota,EVTOL, &
!nevs,scrvec,vecs,num_basis,work,lwork,iwork,ifail,info)
deallocate(work)
deallocate(iwork)
deallocate(ifail)
do i=nroota+1,num_basis
scrvec(i)=1D20
end do


call eigsrt(scrvec,vecs,num_basis,num_basis)
call xgemm( 'N', 'T', num_basis, num_basis, numalpha, one,vecs(1:num_basis,indexa:num_basis) &
 , num_basis,vecs(1:num_basis,indexa:num_basis) , num_basis,zero,scrmat,num_basis )
call mkvector(scrmat,adens,num_basis)


!end if





density=adens+bdens




pdelta=zero
do i=1,ndim1
pdelta=pdelta+(density(i)-scrvec3(i))**2
end do
! compute energy
go to 1

20 continue
!stdout write(*,*)'***********************************'
!stdoutwrite(*,*)'            SCF RESULTS               '
etotal=energy+enuc
!stdoutprint*,'Nuclear repulsion energy (eV) = ',enuc
!stdoutprint*,'Electronic energy (eV) = ',energy
!stdoutprint*,'SCF Total energy (eV) = ',etotal
!stdoutprint*,'homo energy (eV) = ',scrvec(indexa)
call square(scrmat,density,num_basis)
call trace(scrmat,num_basis,out)
call square(scrmat,adens,num_basis)
call trace(scrmat,num_basis,out1)
call square(scrmat,bdens,num_basis)
call trace(scrmat,num_basis,out2)
!stdoutprint*,'Total # electrons at convergence (Tr P)= ',nint(out)
!stdoutprint*,'Total # alpha electrons at convergence [Tr P(alpha)] = ',nint(out1)
!stdoutprint*,'Total # beta electrons at convergence [Tr P(beta)] = ',nint(out2)
!stdoutwrite(*,*)'***********************************'
!call matprt2(fock,num_basis,num_basis,num_basis,num_basis,scrmat)
if(keep)then
    if(allocated(s3))then
    deallocate(s3)
    end if
  allocate(s3(ndim1))
  s3=density
end if





deallocate(focka)
deallocate(fockb)
deallocate(density)
deallocate(vecs)
deallocate(vec2)
deallocate(scrvec)
deallocate(scrmat)
deallocate(scrvec3)


!10 format(i4,f15.8,f15.8,f15.8)

end subroutine uhf

 


