subroutine buildgto(contr,expon,lgto,ngau,ngto,goverlap,ci,cj,ei,ej,coef, &
coorb,expa)
use constants
implicit double precision(a-h,o-z)
double precision,intent(inout),dimension(:,:)::goverlap,coef,expa
double precision,intent(in),dimension(:,:)::contr,expon,coorb
integer,intent(in),dimension(:)::lgto,ngau
integer,intent(in)::ngto
double precision,dimension(:),intent(inout)::ci,cj,ei,ej
double precision,dimension(3,2)::xyz
integer,dimension(1,2)::imax

interface
subroutine ssinitial(coef,expa,xyz,ssout,n1,n2)
double precision,dimension(:,:),intent(in)::coef,expa,xyz
integer,intent(in)::n1,n2
double precision,intent(inout)::ssout
end subroutine ssinitial

subroutine psinitial(coef,expa,xyz,psout,lone,ltwo,imax)
double precision,dimension(:,:),intent(in)::coef,expa,xyz
integer,intent(in)::lone,ltwo
integer,intent(in),dimension(:,:)::imax
double precision,intent(inout)::psout
end subroutine psinitial

subroutine dsinitial(coef,expa,xyz,dsout,lone,ltwo,imax)
use constants
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(in)::coef,expa,xyz
integer,intent(in)::lone,ltwo
integer,intent(in),dimension(:,:)::imax
double precision,intent(inout)::dsout
end subroutine dsinitial

subroutine ppinitial(coef,expa,xyz,ppout,lone,ltwo,imax)
use constants
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(in)::coef,expa,xyz
integer,intent(in)::lone,ltwo
integer,intent(in),dimension(:,:)::imax
double precision,intent(inout)::ppout
end subroutine ppinitial

subroutine dpinitial(coef,expa,xyz,dpout,lone,ltwo,imax)
use constants
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(in)::coef,expa,xyz
integer,intent(in)::lone,ltwo
integer,intent(in),dimension(:,:)::imax
double precision,intent(inout)::dpout
end subroutine dpinitial

subroutine ddinitial(coef,expa,xyz,ddout,lone,ltwo,imax)
use constants
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(in)::coef,expa,xyz
integer,intent(in)::lone,ltwo
integer,intent(in),dimension(:,:)::imax
double precision,intent(inout)::ddout
end subroutine ddinitial

end interface

do i=1,ngto
itype=lgto(i)
coef(1:ngau(i),1)=contr(1:ngau(i),i)
expa(1:ngau(i),1)=expon(1:ngau(i),i)
xyz(1,1)=coorb(1,i)
xyz(2,1)=coorb(2,i)
xyz(3,1)=coorb(3,i)
imax(1,1)=ngau(i)


do j=1,i


jtype=lgto(j)
coef(1:ngau(j),2)=contr(1:ngau(j),j)
expa(1:ngau(j),2)=expon(1:ngau(j),j)
xyz(1,2)=coorb(1,j)
xyz(2,2)=coorb(2,j)
xyz(3,2)=coorb(3,j)
imax(1,2)=ngau(j)

if(itype.eq.1)then
   if(jtype.eq.1)then
   call ssinitial(coef,expa,xyz,goverlap(j,i),ngau(i),ngau(j))
!   write(*,*)j,i,goverlap(j,i),'ss'
   elseif(jtype.ge.2 .and. jtype.le.4)then
   call psinitial(coef,expa,xyz,goverlap(j,i),itype,jtype,imax)
!   write(*,*)j,i,goverlap(j,i),'sp'
   else
   call dsinitial(coef,expa,xyz,goverlap(j,i),itype,jtype,imax)
!   write(*,*)j,i,goverlap(j,i),'sd'
   end if


elseif(itype.ge.2 .and. itype.le.4)then
   if(jtype.eq.1)then
      call psinitial(coef,expa,xyz,goverlap(j,i),itype,jtype,imax)
!   write(*,*)j,i,goverlap(j,i),'ps'
   elseif(jtype.ge.2 .and. jtype.le.4)then
      call ppinitial(coef,expa,xyz,goverlap(j,i),itype,jtype,imax)
!      write(*,*)j,i,goverlap(j,i),'pp'
   else
   call dpinitial(coef,expa,xyz,goverlap(j,i),itype,jtype,imax)
!   write(*,*)j,i,goverlap(j,i),'pd'
      end if

else
   if(jtype.eq.1)then
   call dsinitial(coef,expa,xyz,goverlap(j,i),itype,jtype,imax)
!   write(*,*)j,i,goverlap(j,i),'ds'
   elseif(jtype.ge.2 .and. jtype.le.4)then
   call dpinitial(coef,expa,xyz,goverlap(j,i),itype,jtype,imax)
!   write(*,*)j,i,goverlap(j,i),'pd'
   else
   call ddinitial(coef,expa,xyz,goverlap(j,i),itype,jtype,imax)
!   write(*,*)j,i,goverlap(j,i),'dd'
      end if




end if

end do
end do

do i=1,ngto
do j=1,ngto
goverlap(j,i)=goverlap(i,j)
end do
end do



end subroutine buildgto
