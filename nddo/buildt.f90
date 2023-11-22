subroutine buildt(contr,expon,lgto,ngau,ngto,coef, &
coorb,expa,tmatrix)
use tables
use constants
use spherical
use indices
implicit double precision(a-h,o-z)
double precision,intent(in),dimension(:,:)::contr,expon,coorb
integer,intent(in),dimension(:)::lgto,ngau
integer,intent(in)::ngto
double precision,intent(inout),dimension(:,:)::coef,expa,tmatrix
double precision,dimension(6,10)::coefsto,expasto
double precision,dimension(6)::ci,ei
integer:: lsto(10)=(/1,2,3,4,5,6,7,8,9,10/)
double precision,dimension(3,2)::xyz
integer,dimension(1,2)::imax
double precision,allocatable,dimension(:,:)::section,t1

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

!double precision,dimension(:,:),allocatable::t1
!double precision,dimension(:,:),allocatable::t2
ncoef=6

! will have to loop over atoms and do each section
do jatom=1,numat

jcol=nbas(species(jatom))

if(jcol>4)jcol=10

allocate(section(ngto,jcol))
section=zero


call stofit(1,nqs(zeff(jatom)),zs(species(jatom)),ncoef,1,ci,ei)
coefsto(1:ncoef,1)=ci
expasto(1:ncoef,1)=ei

if(jcol>1)then
call stofit(2,nqp(zeff(jatom)),zp(species(jatom)),ncoef,2,ci,ei)
coefsto(1:ncoef,2)=ci
expasto(1:ncoef,2)=ei
coefsto(1:ncoef,3)=ci
expasto(1:ncoef,3)=ei
coefsto(1:ncoef,4)=ci
expasto(1:ncoef,4)=ei


if(jcol>4)then
call stofit(3,nqd(zeff(jatom)),zetad(species(jatom)),ncoef,5,ci,ei)
coefsto(1:ncoef,5)=ci
expasto(1:ncoef,5)=ei
coefsto(1:ncoef,6)=ci
expasto(1:ncoef,6)=ei
coefsto(1:ncoef,7)=ci
expasto(1:ncoef,7)=ei
call stofit(3,nqd(zeff(jatom)),zetad(species(jatom)),ncoef,8,ci,ei)
coefsto(1:ncoef,8)=ci
expasto(1:ncoef,8)=ei
coefsto(1:ncoef,9)=ci
expasto(1:ncoef,9)=ei
coefsto(1:ncoef,10)=ci
expasto(1:ncoef,10)=ei
end if
end if


!do i=1,jcol
!write(32,*)lsto(i)
!write(32,*)2
!do j=1,ncoef
!write(30,*)coefsto(j,i)
!write(31,*)expasto(j,i)
!end do
!end do



! so now the sto info is in coefsto. what remains is to loop over 
! sto's and gto's and do the overlaps

xyz(1,1)=x(jatom)/autoang
xyz(2,1)=y(jatom)/autoang
xyz(3,1)=z(jatom)/autoang

do islater=1,jcol ! loop over slaters
itype=lsto(islater)
coef(1:ncoef,1)=coefsto(1:ncoef,islater)
expa(1:ncoef,1)=expasto(1:ncoef,islater)

imax(1,1)=ncoef


do jgauss=1,ngto


jtype=lgto(jgauss)
coef(1:ngau(jgauss),2)=contr(1:ngau(jgauss),jgauss)
expa(1:ngau(jgauss),2)=expon(1:ngau(jgauss),jgauss)
xyz(1,2)=coorb(1,jgauss)
xyz(2,2)=coorb(2,jgauss)
xyz(3,2)=coorb(3,jgauss)
imax(1,2)=ngau(jgauss)

if(itype.eq.1)then
   if(jtype.eq.1)then
   call ssinitial(coef,expa,xyz,section(jgauss,islater),ncoef,ngau(jgauss))
   elseif(jtype.ge.2 .and. jtype.le.4)then
   call psinitial(coef,expa,xyz,section(jgauss,islater),itype,jtype,imax)
   else
   call dsinitial(coef,expa,xyz,section(jgauss,islater),itype,jtype,imax)
   end if


elseif(itype.ge.2 .and. itype.le.4)then
   if(jtype.eq.1)then
      call psinitial(coef,expa,xyz,section(jgauss,islater),itype,jtype,imax)
   elseif(jtype.ge.2 .and. jtype.le.4)then
      call ppinitial(coef,expa,xyz,section(jgauss,islater),itype,jtype,imax)
   else
   call dpinitial(coef,expa,xyz,section(jgauss,islater),itype,jtype,imax)
      end if

else
   if(jtype.eq.1)then
   call dsinitial(coef,expa,xyz,section(jgauss,islater),itype,jtype,imax)
   elseif(jtype.ge.2 .and. jtype.le.4)then
   call dpinitial(coef,expa,xyz,section(jgauss,islater),itype,jtype,imax)
   else
   call ddinitial(coef,expa,xyz,section(jgauss,islater),itype,jtype,imax)
      end if




end if

end do 
end do ! end loops over gaussians and sto basis function


! transform from cartesian d to sperical harmonic d used in mndo

if(jcol==10)then
allocate(t1(ngto,9))
call matprt(section,ngto,10,ngto,10)
call xgemm( 'N', 'N', ngto, 9, 10, one,section, &
 ngto,transpose(trans_columns) , 10,zero,t1,ngto )
deallocate(section)
allocate(section(ngto,9))
section=t1
call matprt(section,ngto,9,ngto,9)
deallocate(t1)
end if

tmatrix(1:ngto,ifirst(jatom):ilast(jatom))=section
deallocate(section)



end do ! end loop over atoms

return
end subroutine buildt
