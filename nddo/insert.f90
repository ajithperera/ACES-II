subroutine insert(ci,ei,cj,ej,NGAUSS,coef,expa,x1,y1,z1,x2,y2,z2,coor)
implicit double precision (a-h,o-z)
double precision,dimension(:),intent(in)::ci,ei,cj,ej
double precision,dimension(:,:),intent(inout)::coef,expa,coor
double precision,intent(in)::x1,y1,z1,x2,y2,z2
integer,intent(in)::ngauss
do i=1,ngauss
coef(i,1)=ci(i)
coef(i,2)=cj(i)
expa(i,1)=ei(i)
expa(i,2)=ej(i)
end do
coor(1,1)=x1
coor(2,1)=y1
coor(3,1)=z1
coor(1,2)=x2
coor(2,2)=y2
coor(3,2)=z2
end subroutine insert
