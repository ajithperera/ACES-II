subroutine overlap(iatom,jatom,iorbs,jorbs,rau,overlcl)
use tables
use constants
implicit double precision (a-h,o-z)
integer,intent(in)::iatom,jatom,iorbs,jorbs
double precision,dimension(:),intent(inout)::overlcl
double precision,intent(in)::rau
double precision,dimension(6)::ci,ei,cj,ej
double precision,dimension(6,2)::coef,expa
double precision,dimension(3,2)::coor
integer,parameter::ncoef=6

!
! storage in overlcl is as follows
!  1=(s|s) 2=(s|o) 3=(o|s) 4=(o|o) 5=(p|p)
!
!
!






interface
subroutine ssover(dist,ssout,NGAUSS,ci,ei,cj,ej)
double precision,dimension(6),intent(in):: ci,ei,cj,ej
integer,intent(in)::ngauss
double precision,intent(in)::dist
double precision,intent(inout)::ssout
end subroutine ssover

subroutine insert(ci,ei,cj,ej,NGAUSS,coef,expa,x1,y1,z1,x2,y2,z2,coor)
double precision,dimension(:),intent(in)::ci,ei,cj,ej
double precision,dimension(:,:),intent(inout)::coef,expa,coor
double precision,intent(in)::x1,y1,z1,x2,y2,z2
integer,intent(in)::ngauss
end subroutine insert


subroutine psover(coor,lone,ltwo,psout,NGAUSS,coef,expa,rau)
double precision,intent(in):: rau
double precision,dimension(:,:),intent(in)::coef,expa,coor
integer,intent(in)::lone,ltwo,ngauss
double precision,intent(inout)::psout
end subroutine psover

subroutine ppover(coor,lone,ltwo,ppout,NGAUSS,coef,expa,rau)
double precision,intent(in):: rau
double precision,dimension(:,:),intent(in)::coef,expa,coor
integer,intent(in)::lone,ltwo,ngauss
double precision,intent(inout)::ppout
end subroutine ppover

end interface

!print*,'overlaps'
! hydrogen hydrogen only has one integral
if(iorbs==1.and.jorbs==1)then
call stofit(1,nqs(zeff(iatom)),zs(species(iatom)),ncoef,1,ci,ei)
call stofit(1,nqs(zeff(jatom)),zs(species(jatom)),ncoef,1,cj,ej)
call ssover(rau,ss,ncoef,ci,ei,cj,ej)
overlcl(1)=ss
return
end if


! hydrogen - heavy atom
if(iorbs==1.and.jorbs==4)then
call stofit(1,nqs(zeff(iatom)),zs(species(iatom)),ncoef,1,ci,ei)
call stofit(1,nqs(zeff(jatom)),zs(species(jatom)),ncoef,1,cj,ej)
call ssover(rau,ss,ncoef,ci,ei,cj,ej)
overlcl(1)=ss

call stofit(2,nqp(zeff(jatom)),zp(species(jatom)),ncoef,4,cj,ej)
call insert(ci,ei,cj,ej,ncoef,coef,expa,zero,zero,zero,zero,zero,rau,coor)
call psover(coor,1,4,psout,ncoef,coef,expa,rau)
overlcl(2)=psout
return
end if


! heavy atom - hydrogen

if(iorbs==4.and.jorbs==1)then
call stofit(1,nqs(zeff(iatom)),zs(species(iatom)),ncoef,1,ci,ei)
call stofit(1,nqs(zeff(jatom)),zs(species(jatom)),ncoef,1,cj,ej)
call ssover(rau,ss,ncoef,ci,ei,cj,ej)
overlcl(1)=ss

call stofit(2,nqp(zeff(iatom)),zp(species(iatom)),ncoef,4,ci,ei)
call insert(ci,ei,cj,ej,ncoef,coef,expa,zero,zero,zero,zero,zero,rau,coor)
call psover(coor,4,1,psout,ncoef,coef,expa,rau)
overlcl(3)=psout
return
end if


! heavy - heavy
if(iorbs==4.and.jorbs==4)then
! (s|s)
call stofit(1,nqs(zeff(iatom)),zs(species(iatom)),ncoef,1,ci,ei)
call stofit(1,nqs(zeff(jatom)),zs(species(jatom)),ncoef,1,cj,ej)
call ssover(rau,ss,ncoef,ci,ei,cj,ej)
overlcl(1)=ss
!print*,'ss',ss

! (s|o)
call stofit(2,nqp(zeff(jatom)),zp(species(jatom)),ncoef,4,cj,ej)
call insert(ci,ei,cj,ej,ncoef,coef,expa,zero,zero,zero,zero,zero,rau,coor)
call psover(coor,1,4,psout,ncoef,coef,expa,rau)
overlcl(2)=psout
!print*,'so',psout
! (o|s)
call stofit(2,nqp(zeff(iatom)),zp(species(iatom)),ncoef,4,ci,ei)
call stofit(1,nqs(zeff(jatom)),zs(species(jatom)),ncoef,1,cj,ej)
call insert(ci,ei,cj,ej,ncoef,coef,expa,zero,zero,zero,zero,zero,rau,coor)
call psover(coor,4,1,psout,ncoef,coef,expa,rau)
overlcl(3)=psout
!print*,'os',psout
! (o|o)
call stofit(2,nqp(zeff(jatom)),zp(species(jatom)),ncoef,4,cj,ej)
call insert(ci,ei,cj,ej,ncoef,coef,expa,zero,zero,zero,zero,zero,rau,coor)
call ppover(coor,4,4,ppout,ncoef,coef,expa,rau)
overlcl(4)=ppout
!print*,'oo',ppout
! (p|p)
call ppover(coor,2,2,ppout,ncoef,coef,expa,rau)
overlcl(5)=ppout
!print*,'pp',ppout
end if














end subroutine overlap
