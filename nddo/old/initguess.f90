subroutine initguess
use tables
use constants
use indices
use scratch_array
implicit double precision(a-h,o-z)
double precision,dimension(:,:),allocatable::contr,expon,coorb, &
coef,expa,goverlap,tmatrix,temp,ginv,rectangle,ttdaginv,vecs,umatrix,d,ddag
integer,dimension(:),allocatable::lgto,ngau,ind
double precision,dimension(:),allocatable::ci,cj,ei,ej
double precision,dimension(3,2)::coor
CHARACTER(2)::LABEL(10)
DATA LABEL/'S','X','Y','Z','XX','YY','ZZ','XY','XZ','YZ'/
 


interface
subroutine normalize(contr,expon,lgto,ngau,ngto)
double precision,intent(inout),dimension(:,:)::contr,expon
integer,intent(in),dimension(:)::lgto,ngau
integer,intent(in)::ngto
end subroutine normalize

subroutine buildgto(contr,expon,lgto,ngau,ngto,goverlap,ci,cj,ei,ej,coef, &
coorb,expa)
double precision,intent(inout),dimension(:,:)::goverlap,coef,expa
double precision,intent(in),dimension(:,:)::contr,expon,coorb
integer,intent(in),dimension(:)::lgto,ngau
integer,intent(in)::ngto
double precision,dimension(:),intent(inout)::ci,cj,ei,ej
end subroutine buildgto

subroutine buildt(contr,expon,lgto,ngau,ngto,coef, &
coorb,expa,tmatrix)
double precision,intent(in),dimension(:,:)::contr,expon,coorb
integer,intent(in),dimension(:)::lgto,ngau
integer,intent(in)::ngto
double precision,intent(inout),dimension(:,:)::coef,expa,tmatrix
end subroutine buildt

end interface

print*,'-----------------------------------------------'
print*,'|                                             |'
print*,'|        PROJECTED MNDO INITIAL GUESS         |'
print*,'-----------------------------------------------'
print*,''
write(*,*)
open(unit=1,file='GTO')
read(1,*)ngto,nrows
write(*,*)'NUMBER OF GAUSSIAN ORBITALS FROM ACESX = ',NGTO

allocate(contr(nrows,ngto))
allocate(expon(nrows,ngto))
allocate(coorb(3,ngto))
allocate(lgto(ngto))
allocate(ngau(ngto))
read(1,*)lgto
read(1,*)ngau
read(1,*)contr
read(1,*)expon
read(1,*)coorb
close(1)


! debug print of basis functions in NDDO subroutine
! SHOULD MATCH THOSE FROM ACES
!      do i=1,ngto
!      write(*,*)'basis function number',i
!      write(*,*)'coordinates',coorb(1,i),coorb(2,i),coorb(3,i)
!      write(*,*)'type  ',label(lgto(i))
!      write(*,*)'------------------------'
!      do j=1,ngau(i)
!      write(*,*)expon(j,i),contr(j,i)
!      end do
!      end do
open(unit=1,file='coef')
      do i=1,ngto
      do j=1,6
      write(1,*)contr(j,i)
      end do
      end do
close(1)

open(unit=1,file='expa')
      do i=1,ngto
      do j=1,6
      write(1,*)expon(j,i)
      end do
      end do
close(1)
open(unit=1,file='ls')
      do i=1,ngto
      write(1,*)lgto(i)
      write(1,*)2
      end do
      
close(1)



allocate(ci(nrows))
allocate(cj(nrows))
allocate(ei(nrows))
allocate(ej(nrows))
allocate(coef(nrows,2))
allocate(expa(nrows,2))

! build the GTO overlap matrix
! probably could read from ACES but i'll just do 
! it myself


allocate(goverlap(ngto,ngto))
call normalize(contr,expon,lgto,ngau,ngto)
call buildgto(contr,expon,lgto,ngau,ngto,goverlap,ci,cj,ei,ej,coef, &
coorb,expa)
open(unit=1,file='goverlap')
write(1,*)goverlap
close(1)

allocate(tmatrix(ngto,num_basis))

call buildt(contr,expon,lgto,ngau,ngto,coef, &
coorb,expa,tmatrix)
allocate(scrmat(ngto,ngto))
call dgemm( 'N', 'T', ngto, ngto, num_basis, one,tmatrix &
 , ngto,tmatrix , ngto,zero,scrmat,ngto ) ! so scrmat hold TT^+
! invert TT+
allocate(umatrix(ngto,ngto))
call eig(scrmat,umatrix,ngto,ngto,1)

do i=1,ngto
if(abs(scrmat(i,i)).lt.1D-10)then
scrmat(i,i)=1D-4
end if
scrmat(i,i)=one/scrmat(i,i)
end do

allocate(rectangle(num_basis,ngto))
allocate(vecs(num_basis,num_basis))
open(unit=1,file='vecs')
read(1,*)vecs
close(1)
rectangle=zero
call dgemm( 'T', 'T', num_basis, ngto, num_basis, one,vecs &
 , num_basis,tmatrix ,ngto,zero,rectangle,num_basis )

deallocate(goverlap)
allocate(goverlap(num_basis,ngto))
call dgemm( 'N', 'N', num_basis, ngto, ngto, one,rectangle &
 , num_basis,umatrix , ngto,zero,goverlap,num_basis )

rectangle=zero
call dgemm( 'N', 'N', num_basis, ngto, ngto, one,goverlap &
 , num_basis,scrmat , ngto,zero,rectangle,num_basis )


goverlap=zero
call dgemm( 'N', 'T', num_basis, ngto, ngto, one,rectangle &
 , num_basis,umatrix , ngto,zero,goverlap,num_basis )

allocate(d(ngto,num_basis))
d=transpose(goverlap)


! check
deallocate(scrmat)
allocate(scrmat(num_basis,num_basis))
call dgemm( 'T', 'N', num_basis, num_basis, ngto, one,d &
 , ngto,tmatrix , ngto,zero,scrmat,num_basis )
deallocate(rectangle)
allocate(rectangle(num_basis,num_basis))
call dgemm( 'N', 'N', num_basis, num_basis, num_basis, one,scrmat &
 , num_basis,vecs , num_basis,zero,rectangle,num_basis )



print*,'please be unity'
call matprt(rectangle,num_basis,num_basis,num_basis,num_basis )








deallocate(scrmat)
allocate(scrmat(ngto,ngto))
open(unit=1,file='goverlap')
read(1,*)scrmat
close(1)
deallocate(rectangle)
allocate(rectangle(num_basis,ngto))

call dgemm( 'T', 'N', num_basis,ngto , ngto, one,d &
 , ngto,scrmat , ngto,zero,rectangle,num_basis )
vecs=zero
call dgemm( 'N', 'N', num_basis, num_basis, ngto, one,rectangle &
 , num_basis,d , ngto,zero,vecs,num_basis )

call matprt(vecs,num_basis,num_basis,num_basis,num_basis)

!do i=1,num_basis
!d(1:ngto,i:i)=d(1:ngto,i:i)/dsqrt(vecs(i,i))
!end do

deallocate(scrmat)
allocate(scrmat(ngto,ngto))
open(unit=1,file='goverlap')
read(1,*)scrmat
close(1)
deallocate(rectangle)
allocate(rectangle(num_basis,ngto))

call dgemm( 'T', 'N', num_basis,ngto , ngto, one,d &
 , ngto,scrmat , ngto,zero,rectangle,num_basis )
vecs=zero
call dgemm( 'N', 'N', num_basis, num_basis, ngto, one,rectangle &
 , num_basis,d , ngto,zero,vecs,num_basis )

call matprt(vecs,num_basis,num_basis,num_basis,num_basis)


! make density
deallocate(scrmat)
allocate(scrmat(ngto,ngto))
icol=nelectrons/2
ihigh=num_basis-(nelectrons/2)+1
print*,'indices',icol,ihigh
call dgemm( 'N', 'T', ngto, ngto, icol, two,d(1:ngto,ihigh:num_basis) &
 , ngto,d(1:ngto,ihigh:num_basis) , ngto,zero,scrmat,ngto )




deallocate(rectangle)
allocate(rectangle(ngto,ngto))
open(unit=1,file='goverlap')
read(1,*)rectangle
close(1)
deallocate(vecs)
allocate(vecs(ngto,ngto))

call dgemm( 'N', 'N', ngto,ngto , ngto, one,scrmat &
 , ngto,rectangle , ngto,zero,vecs,ngto )
total=zero
do i=1,ngto
total=total+vecs(i,i)
end do
print*,'projected count',total

open(unit=1,file='density')
write(1,*)scrmat
close(1)
















stop


end subroutine initguess


