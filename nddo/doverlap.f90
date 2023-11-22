 subroutine doverlap(iatom,jatom,rau,dlocal)
use tables
use constants
use spherical
implicit double precision(a-h,o-z)
double precision,dimension(6,4)::coefi,coefj,expai,expaj
double precision,dimension(6)::ci,ei,cj,ej
double precision,intent(in)::rau
double precision,intent(inout),dimension(:,:)::dlocal
integer,intent(in)::iatom,jatom
integer,parameter::ncoef=6
integer:: label(10)=(/1,2,3,4,5,6,7,8,9,10/)
integer:: map(10)=(/1,2,2,2,3,3,3,3,3,3/)
double precision,dimension(3,2)::coor
double precision,dimension(6,2)::coef,expa
double precision,allocatable,dimension(:,:)::cartesian

double precision,dimension(:,:),allocatable::t1
double precision,dimension(:,:),allocatable::t2



interface
subroutine insert(ci,ei,cj,ej,NGAUSS,coef,expa,x1,y1,z1,x2,y2,z2,coor)
double precision,dimension(:),intent(in)::ci,ei,cj,ej
double precision,dimension(:,:),intent(inout)::coef,expa,coor
double precision,intent(in)::x1,y1,z1,x2,y2,z2
integer,intent(in)::ngauss
end subroutine insert

subroutine dsover(coor,lone,ltwo,dsout,NGAUSS,coef,expa,rau)
double precision,intent(in):: rau
double precision,dimension(:,:),intent(in)::coef,expa,coor
integer,intent(in)::lone,ltwo,ngauss
double precision,intent(inout)::dsout
end subroutine dsover

      subroutine dpover(coor,lone,ltwo,dpout,NGAUSS,coef,expa,rau)
double precision,intent(in):: rau
double precision,dimension(:,:),intent(in)::coef,expa,coor
integer,intent(in)::lone,ltwo,ngauss
double precision,intent(inout)::dpout
end subroutine dpover

      subroutine ddover(coor,lone,ltwo,ddout,NGAUSS,coef,expa,rau)
double precision,intent(in):: rau
double precision,dimension(:,:),intent(in)::coef,expa,coor
integer,intent(in)::lone,ltwo,ngauss
double precision,intent(inout)::ddout
end subroutine ddover
end interface
irow=nbas(species(iatom))
if(irow>4)irow=10
jcol=nbas(species(jatom))
if(jcol>4)jcol=10
allocate(cartesian(irow,jcol))
cartesian=zero
call stofit(1,nqs(zeff(iatom)),zs(species(iatom)),ncoef,1,ci,ei)
coefi(1:ncoef,1)=ci
expai(1:ncoef,1)=ei
if(irow>1)then
call stofit(2,nqp(zeff(iatom)),zp(species(iatom)),ncoef,2,ci,ei)
coefi(1:ncoef,2)=ci
expai(1:ncoef,2)=ei
if(nbas(species(iatom))>4)then
call stofit(3,nqd(zeff(iatom)),zetad(species(iatom)),ncoef,5,ci,ei)
coefi(1:ncoef,3)=ci
expai(1:ncoef,3)=ei
call stofit(3,nqd(zeff(iatom)),zetad(species(iatom)),ncoef,8,ci,ei)
coefi(1:ncoef,4)=ci
expai(1:ncoef,4)=ei
end if
end if


call stofit(1,nqs(zeff(jatom)),zs(species(jatom)),ncoef,1,ci,ei)
coefj(1:ncoef,1)=ci
expaj(1:ncoef,1)=ei
if(jcol>1)then
call stofit(2,nqp(zeff(jatom)),zp(species(jatom)),ncoef,2,ci,ei)
coefj(1:ncoef,2)=ci
expaj(1:ncoef,2)=ei
if(nbas(species(jatom))>4)then
call stofit(3,nqd(zeff(jatom)),zetad(species(jatom)),ncoef,5,ci,ei)
coefj(1:ncoef,3)=ci
expaj(1:ncoef,3)=ei
call stofit(3,nqd(zeff(jatom)),zetad(species(jatom)),ncoef,8,ci,ei)
coefj(1:ncoef,4)=ci
expaj(1:ncoef,4)=ei
end if
end if





! do the s(i) d(j) block
if(nbas(species(jatom)).gt.4)then
ci=coefi(1:ncoef,1)
ei=expai(1:ncoef,1)
cj=coefj(1:ncoef,3)
ej=expaj(1:ncoef,3)
call insert(ci,ei,cj,ej,ncoef,coef,expa,x(iatom),y(iatom),z(iatom),x(jatom),y(jatom),z(jatom),coor)
coor=coor/autoang
icount=1
do j=5,7
call dsover(coor,1,label(j),dsout,ncoef,coef,expa,rau)
cartesian(icount,j)=dsout
end do

cj=coefj(1:ncoef,4)
ej=expaj(1:ncoef,4)
call insert(ci,ei,cj,ej,ncoef,coef,expa,x(iatom),y(iatom),z(iatom),x(jatom),y(jatom),z(jatom),coor)
coor=coor/autoang
do j=8,10
call dsover(coor,1,label(j),dsout,ncoef,coef,expa,rau)
cartesian(icount,j)=dsout
end do
end if 


if(nbas(species(iatom)).ge.4)then

if(nbas(species(jatom)).gt.4)then
! do the p d block
ci=coefi(1:ncoef,2)
ei=expai(1:ncoef,2)
cj=coefj(1:ncoef,3)
ej=expaj(1:ncoef,3)
call insert(ci,ei,cj,ej,ncoef,coef,expa,x(iatom),y(iatom),z(iatom),x(jatom),y(jatom),z(jatom),coor)
coor=coor/autoang
do k=5,7
do j=2,4
call dpover(coor,label(j),label(k),dpout,ncoef,coef,expa,rau)
cartesian(j,k)=dpout
end do
end do

cj=coefj(1:ncoef,4)
ej=expaj(1:ncoef,4)
call insert(ci,ei,cj,ej,ncoef,coef,expa,x(iatom),y(iatom),z(iatom),x(jatom),y(jatom),z(jatom),coor)
coor=coor/autoang
do k=8,10
do j=2,4
call dpover(coor,label(j),label(k),dpout,ncoef,coef,expa,rau)
cartesian(j,k)=dpout
end do
end do

end if
end if

if(nbas(species(iatom)).gt.4)then
! do the d s block
icount=1
ci=coefi(1:ncoef,3)
ei=expai(1:ncoef,3)
cj=coefj(1:ncoef,1)
ej=expaj(1:ncoef,1)
call insert(ci,ei,cj,ej,ncoef,coef,expa,x(iatom),y(iatom),z(iatom),x(jatom),y(jatom),z(jatom),coor)
coor=coor/autoang
do j=5,7
call dsover(coor,label(j),1,dsout,ncoef,coef,expa,rau)
cartesian(j,1)=dsout
end do
ci=coefi(1:ncoef,4)
ei=expai(1:ncoef,4)
call insert(ci,ei,cj,ej,ncoef,coef,expa,x(iatom),y(iatom),z(iatom),x(jatom),y(jatom),z(jatom),coor)
coor=coor/autoang
do j=8,10
call dsover(coor,label(j),1,dsout,ncoef,coef,expa,rau)
cartesian(j,1)=dsout
end do


if(nbas(species(jatom)).ge.4)then
! do the d p block
ci=coefi(1:ncoef,3)
ei=expai(1:ncoef,3)
cj=coefj(1:ncoef,2)
ej=expaj(1:ncoef,2)
call insert(ci,ei,cj,ej,ncoef,coef,expa,x(iatom),y(iatom),z(iatom),x(jatom),y(jatom),z(jatom),coor)
coor=coor/autoang
do j=2,4
do k=5,7
call dpover(coor,label(k),label(j),dpout,ncoef,coef,expa,rau)
cartesian(k,j)=dpout
end do
end do
ci=coefi(1:ncoef,4)
ei=expai(1:ncoef,4)
call insert(ci,ei,cj,ej,ncoef,coef,expa,x(iatom),y(iatom),z(iatom),x(jatom),y(jatom),z(jatom),coor)
coor=coor/autoang
do j=2,4
do k=8,10
call dpover(coor,label(k),label(j),dpout,ncoef,coef,expa,rau)
cartesian(k,j)=dpout
end do
end do

end if


if(nbas(species(jatom)).gt.4)then
! do the dd block
ci=coefi(1:ncoef,3)
ei=expai(1:ncoef,3)
cj=coefj(1:ncoef,3)
ej=expaj(1:ncoef,3)
call insert(ci,ei,cj,ej,ncoef,coef,expa,x(iatom),y(iatom),z(iatom),x(jatom),y(jatom),z(jatom),coor)
coor=coor/autoang
do j=5,7
do k=5,7
call ddover(coor,label(k),label(j),ddout,ncoef,coef,expa,rau)
cartesian(k,j)=ddout
end do
end do

cj=coefj(1:ncoef,4)
ej=expaj(1:ncoef,4)
call insert(ci,ei,cj,ej,ncoef,coef,expa,x(iatom),y(iatom),z(iatom),x(jatom),y(jatom),z(jatom),coor)
coor=coor/autoang
do j=8,10
do k=5,7
call ddover(coor,label(k),label(j),ddout,ncoef,coef,expa,rau)
cartesian(k,j)=ddout
end do
end do

ci=coefi(1:ncoef,4)
ei=expai(1:ncoef,4)
cj=coefj(1:ncoef,3)
ej=expaj(1:ncoef,3)
call insert(ci,ei,cj,ej,ncoef,coef,expa,x(iatom),y(iatom),z(iatom),x(jatom),y(jatom),z(jatom),coor)
coor=coor/autoang
do j=5,7
do k=8,10
call ddover(coor,label(k),label(j),ddout,ncoef,coef,expa,rau)
cartesian(k,j)=ddout
end do
end do

ci=coefi(1:ncoef,4)
ei=expai(1:ncoef,4)
cj=coefj(1:ncoef,4)
ej=expaj(1:ncoef,4)
call insert(ci,ei,cj,ej,ncoef,coef,expa,x(iatom),y(iatom),z(iatom),x(jatom),y(jatom),z(jatom),coor)
coor=coor/autoang
do j=8,10
do k=8,10
call ddover(coor,label(k),label(j),ddout,ncoef,coef,expa,rau)
cartesian(k,j)=ddout
end do
end do
end if
end if


!call matprt(cartesian,irow,jcol,irow,jcol)


if(irow==10.and.jcol==4)then
allocate(t1(9,4))
call xgemm( 'N', 'N', 9, 4, 10, one,trans_columns, &
 9,cartesian , 10,zero,t1,9 )
dlocal=t1
deallocate(t1)

elseif(irow==4.and.jcol==10)then
allocate(t1(4,9))
call xgemm( 'N', 'N', 4, 9, 10, one,cartesian, &
 4,transpose(trans_columns) , 10,zero,t1,4 )
dlocal=t1
deallocate(t1)

elseif(irow==10.and.jcol==10)then
allocate(t1(9,10))
allocate(t2(9,9))

call xgemm( 'N', 'N', 9, 10, 10, one,trans_columns, &
 9,cartesian , 10,zero,t1,9 )
call xgemm( 'N', 'N', 9, 9, 10, one,t1, &
 9,transpose(trans_columns), 10,zero,t2,9 )
dlocal=t2
deallocate(t2)

elseif(irow==1.and.jcol==10)then
allocate(t1(1,9))
call xgemm( 'N', 'N', 1, 9, 10, one,cartesian, &
 1,transpose(trans_columns) , 10,zero,t1,1 )
dlocal=t1
deallocate(t1)

elseif(irow==10.and.jcol==1)then
allocate(t1(9,4))
call xgemm( 'N', 'N', 9, 1, 10, one,trans_columns, &
 9,cartesian , 10,zero,t1,9 )
dlocal=t1
deallocate(t1)


end if




deallocate(cartesian)
return
end subroutine doverlap




















