subroutine initguess1
use control
use tables
use constants
use indices
use scratch_array
implicit double precision(a-h,o-z)
double precision,dimension(:,:),allocatable::contr,expon,coorb, &
coef,expa,goverlap,tmatrix,temp,ginv,rectangle,rectangle2
integer,dimension(:),allocatable::lgto,ngau,ind,istart,iend
double precision,dimension(:),allocatable::ci,cj,ei,ej,vec2
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


open(unit=1,file='GTO')
read(1,*)ngto,nrows
!write(*,*)'NUMBER OF GAUSSIAN ORBITALS FROM ACES= ',NGTO
!write(*,*)'@Project NDDO'
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
! allocate and build T matrix to hold overlaps between gto and sto
allocate(tmatrix(ngto,num_basis))
call buildt(contr,expon,lgto,ngau,ngto,coef, &
coorb,expa,tmatrix)
! compute S(-1)
allocate(ind(ngto))
allocate(ginv(ngto,ngto))
call migs(goverlap,ngto,ginv,ind)
allocate(temp(ngto,ngto))
temp=zero
! dont need this, just checking inverse
!call buildgto(contr,expon,lgto,ngau,ngto,goverlap,ci,cj,ei,ej,coef, &
!coorb,expa)
!call xgemm( 'N', 'N', ngto, ngto, ngto, one,goverlap &
! , ngto,ginv,ngto,zero,temp,ngto)
! temp should unity
!call matprt(temp,ngto,ngto,ngto,ngto)
! end of stuff not needed!!!!!!!!!!


! now do the rhf or uhf projections
jlength=ngto*(ngto+1)/2 ! this is the dimension of the upper triangle of the target matrix
CALL GETREC(-1,'JOBARC','NBASTOT ',1,MBASCOMP)! number of fns in comp basis
  if(restricted)then
       ! store rhf nddo density (held in s3) as a matrix
      deallocate(temp)
      allocate(temp(num_basis,num_basis))
      temp=zero
      call square(temp,s3,num_basis)
! do the multiplications: the equation is (S^-1)*T*P*(T^+)*(S^-1)=Q
! where Q is the projected nddo density to be written to jobarc
       allocate(rectangle(ngto,num_basis))
       call xgemm( 'N', 'N', ngto, num_basis, ngto, one,ginv &
       , ngto,tmatrix,ngto,zero,rectangle,ngto)
        allocate(rectangle2(ngto,num_basis))
        call xgemm( 'N', 'N', ngto, num_basis, num_basis, one,rectangle &
       ,    ngto,temp,num_basis,zero,rectangle2,ngto)

        deallocate(rectangle)
        allocate(rectangle(ngto,ngto))
        call xgemm( 'N', 'T', ngto, ngto, num_basis, one,rectangle2 &
       , ngto,tmatrix,ngto,zero,rectangle,ngto)


        deallocate(rectangle2)
        allocate(rectangle2(ngto,ngto))
        call xgemm( 'N', 'N', ngto, ngto, ngto, one,rectangle &
         , ngto,ginv,ngto,zero,rectangle2,ngto)

if(nirrep.gt.1)then ! if not C1 transform to SO basis
   deallocate(temp)
   deallocate(rectangle)
! MAKE NGTO NSO once you figure this stuff out
   allocate(temp(ngto,ngto))
   allocate(rectangle(mbascomp,ngto))  
deallocate(ginv)
allocate(ginv(mbascomp,ngto))
ginv=zero
! GET THE REQUIRED MATRICES
! TEMP HOLDS THE AO->SO TRANSFORMATION
! GINV HOLD THE ZMAT->CMP ORDER TRANSFORMATION
!
CALL GETREC(20,'JOBARC','ZMAT2CMP',MBASCOMP*NGTO*IRATIO,ginv)
  call xgemm( 'N', 'N', mbascomp, ngto, ngto, one,ginv &
 , mbascomp,rectangle2,ngto,zero,rectangle,mbascomp)
deallocate(rectangle2)
allocate(rectangle2(mbascomp,mbascomp))
   call xgemm( 'N', 'T', mbascomp, mbascomp, ngto, one,rectangle &
 , mbascomp,ginv,mbascomp,zero,rectangle2,mbascomp)






! now include logic to include lower half of each  non-zero block
	!first get offsets
	allocate(istart(nirrep))
	allocate(iend(nirrep))
	idum=0
	jdum=0
	istart=0
	iend=0
	do i=1,nirrep
	if(nbfirr(i).ne.0)then
	idum=jdum+1
	istart(i)=idum
	jdum=idum+nbfirr(i)-1
	iend(i)=jdum
	end if
	end do
!just to debug
!do i=1,nirrep
!write(*,*)'irrep =',i,'  # functions =',nbfirr(i),'with offsets ',istart(i),iend(i)
!end do

!compute length of vector to hold lower triangle of nonzero blocks
	jlength=0
	do i=1,nirrep
	jlength=jlength+nbfirr(i)*(nbfirr(i)+1)/2
	end do
allocate(scrvec(jlength))
! put the density in the vector
	icount=1
	do i=1,nirrep
	if(nbfirr(i).ne.0)then
	do j=istart(i),iend(i)
	do k=istart(i),j
	scrvec(icount)=rectangle2(k,j)
	icount=icount+1
	end do
	end do
	end if
	end do

CALL PUTREC(20,'JOBARC','NDDODENA',jlength*iratio,scrvec)
return
else

!store upper triangle and put the density on jobarc!
DEALLOCATE(GINV)
ALLOCATE(GINV(MBASCOMP,NGTO))
DEALLOCATE(RECTANGLE)
ALLOCATE(RECTANGLE(MBASCOMP,NGTO))
CALL GETREC(20,'JOBARC','ZMAT2CMP',MBASCOMP*NGTO*IRATIO,ginv)


  call xgemm( 'N', 'N', mbascomp, ngto, ngto, one,ginv &
 , mbascomp,rectangle2,ngto,zero,rectangle,mbascomp)

deallocate(rectangle2)
allocate(rectangle2(mbascomp,mbascomp))
   call xgemm( 'N', 'T', mbascomp, mbascomp, ngto, one,rectangle &
 , mbascomp,ginv,mbascomp,zero,rectangle2,mbascomp)




	jlength=mbascomp*(mbascomp+1)/2
   allocate(scrvec(jlength))
   call mkvector(rectangle2,scrvec,mbascomp)
   CALL PUTREC(20,'JOBARC','NDDODENA',jlength*iratio,scrvec)

end if



else ! do a uhf projection

   ! store uhf/nddo alpha density as a matrix called temp
   if(allocated(temp))deallocate(temp)
   allocate(temp(num_basis,num_basis))
   temp=zero
   call square(temp,adens,num_basis)

! do the multiplications: the equation is (S^-1)*T*P*(T^+)*(S^-1)=Q
! where Q is the projected nddo density to be written to jobarc
   if(allocated(rectangle))deallocate(rectangle)
   allocate(rectangle(ngto,num_basis))
   call xgemm( 'N', 'N', ngto, num_basis, ngto, one,ginv &
 , ngto,tmatrix,ngto,zero,rectangle,ngto)


   if(allocated(rectangle2))deallocate(rectangle2)
   allocate(rectangle2(ngto,num_basis))
   call xgemm( 'N', 'N', ngto, num_basis, num_basis, one,rectangle &
 , ngto,temp,num_basis,zero,rectangle2,ngto)

   deallocate(rectangle)
   allocate(rectangle(ngto,ngto))
   call xgemm( 'N', 'T', ngto, ngto, num_basis, one,rectangle2 &
 , ngto,tmatrix,ngto,zero,rectangle,ngto)


   deallocate(rectangle2)
   allocate(rectangle2(ngto,ngto))
   call xgemm( 'N', 'N', ngto, ngto, ngto, one,rectangle &
 , ngto,ginv,ngto,zero,rectangle2,ngto)



if(nirrep.gt.1)then ! if not C1 transform to SO basis
   deallocate(temp)
   deallocate(rectangle)
allocate(temp(mbascomp,ngto))
allocate(rectangle(mbascomp,ngto))

CALL GETREC(20,'JOBARC','ZMAT2CMP',MBASCOMP*NGTO*IRATIO,temp)

  call xgemm( 'N', 'N', mbascomp, ngto, ngto, one,temp &
 , mbascomp,rectangle2,ngto,zero,rectangle,mbascomp)
deallocate(rectangle2)
allocate(rectangle2(mbascomp,mbascomp))
   call xgemm( 'N', 'T', mbascomp, mbascomp, ngto, one,rectangle &
 , mbascomp,temp,mbascomp,zero,rectangle2,mbascomp)






! now include logic to include lower half of each  non-zero block
	!first get offsets
	allocate(istart(nirrep))
	allocate(iend(nirrep))
	idum=0
	jdum=0
	istart=0
	iend=0
	do i=1,nirrep
	if(nbfirr(i).ne.0)then
	idum=jdum+1
	istart(i)=idum
	jdum=idum+nbfirr(i)-1
	iend(i)=jdum
	end if
	end do
!just to debug
do i=1,nirrep
write(*,*)'irrep =',i,'  # functions =',nbfirr(i),'with offsets ',istart(i),iend(i)
end do

!compute length of vector to hold lower triangle of nonzero blocks
	jlength=0
	do i=1,nirrep
	jlength=jlength+nbfirr(i)*(nbfirr(i)+1)/2
	end do
allocate(scrvec(jlength))
! put the density in the vector
	icount=1
	do i=1,nirrep
	if(nbfirr(i).ne.0)then
	do j=istart(i),iend(i)
	do k=istart(i),j
	scrvec(icount)=rectangle2(k,j)
	icount=icount+1
	end do
	end do
	end if
	end do

CALL PUTREC(20,'JOBARC','NDDODENA',jlength*iratio,scrvec)

else
deallocate(temp)
allocate(temp(mbascomp,ngto))
CALL GETREC(20,'JOBARC','ZMAT2CMP',MBASCOMP*NGTO*IRATIO,temp)
deallocate(rectangle)
allocate(rectangle(mbascomp,ngto))

  call xgemm( 'N', 'N', mbascomp, ngto, ngto, one,temp &
 , mbascomp,rectangle2,ngto,zero,rectangle,mbascomp)
deallocate(rectangle2)
allocate(rectangle2(mbascomp,mbascomp))
   call xgemm( 'N', 'T', mbascomp, mbascomp, ngto, one,rectangle &
 , mbascomp,temp,mbascomp,zero,rectangle2,mbascomp)




        jlength=mbascomp*(mbascomp+1)/2
   allocate(scrvec(jlength))
  call mkvector(rectangle2,scrvec,mbascomp)
   CALL PUTREC(20,'JOBARC','NDDODENA',jlength*iratio,scrvec)


end if





 ! and now the beta projection



! store uhf/nddo beta density as a matrix called temp
   if(allocated(temp))deallocate(temp)
   allocate(temp(num_basis,num_basis))
   temp=zero
   call square(temp,bdens,num_basis)

! do the multiplications: the equation is (S^-1)*T*P*(T^+)*(S^-1)=Q
! where Q is the projected nddo density to be written to jobarc
   if(allocated(rectangle))deallocate(rectangle)
   allocate(rectangle(ngto,num_basis))
   call xgemm( 'N', 'N', ngto, num_basis, ngto, one,ginv &
 , ngto,tmatrix,ngto,zero,rectangle,ngto)


   if(allocated(rectangle2))deallocate(rectangle2)
   allocate(rectangle2(ngto,num_basis))
   call xgemm( 'N', 'N', ngto, num_basis, num_basis, one,rectangle &
 , ngto,temp,num_basis,zero,rectangle2,ngto)

   deallocate(rectangle)
   allocate(rectangle(ngto,ngto))
   call xgemm( 'N', 'T', ngto, ngto, num_basis, one,rectangle2 &
 , ngto,tmatrix,ngto,zero,rectangle,ngto)


   deallocate(rectangle2)
   allocate(rectangle2(ngto,ngto))
   call xgemm( 'N', 'N', ngto, ngto, ngto, one,rectangle &
 , ngto,ginv,ngto,zero,rectangle2,ngto)

   !put the beta density on jobarc and we are done


if(nirrep.gt.1)then ! if not C1 transform to SO basis
   deallocate(temp)
   deallocate(rectangle)
! MAKE NGTO NSO once you figure this stuff out
   allocate(temp(mbascomp,ngto))
   allocate(rectangle(mbascomp,ngto))  
CALL GETREC(20,'JOBARC','ZMAT2CMP',mbascomp*NGTO*IRATIO,temp)

  call xgemm( 'N', 'N', mbascomp, ngto, ngto, one,temp &
 , mbascomp,rectangle2,ngto,zero,rectangle,mbascomp)
deallocate(rectangle2)
allocate(rectangle2(mbascomp,mbascomp))
   call xgemm( 'N', 'T', mbascomp, mbascomp, ngto, one,rectangle &
 , mbascomp,temp,mbascomp,zero,rectangle2,mbascomp)






! now include logic to include lower half of each  non-zero block
	!first get offsets

	icount=1
	do i=1,nirrep
	if(nbfirr(i).ne.0)then
	do j=istart(i),iend(i)
	do k=istart(i),j
	scrvec(icount)=rectangle2(k,j)
	icount=icount+1
	end do
	end do
	end if
	end do

CALL PUTREC(20,'JOBARC','NDDODENB',jlength*iratio,scrvec)
else

!store upper triangle and put the density on jobarc!
deallocate(temp)
allocate(temp(mbascomp,ngto))
CALL GETREC(20,'JOBARC','ZMAT2CMP',MBASCOMP*NGTO*IRATIO,temp)
deallocate(rectangle)
allocate(rectangle(mbascomp,ngto))

  call xgemm( 'N', 'N', mbascomp, ngto, ngto, one,temp &
 , mbascomp,rectangle2,ngto,zero,rectangle,mbascomp)
deallocate(rectangle2)
allocate(rectangle2(mbascomp,mbascomp))
   call xgemm( 'N', 'T', mbascomp, mbascomp, ngto, one,rectangle &
 , mbascomp,temp,mbascomp,zero,rectangle2,mbascomp)




	jlength=mbascomp*(mbascomp+1)/2
!   allocate(scrvec(jlength))
  call mkvector(rectangle2,scrvec,mbascomp)
   CALL PUTREC(20,'JOBARC','NDDODENB',jlength*iratio,scrvec)
 end if




end if

return



end subroutine initguess1
