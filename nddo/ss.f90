subroutine ssover(dist,ssout,NGAUSS,ci,ei,cj,ej)
use constants
implicit double precision (a-h,o-z)
double precision,dimension(6),intent(in):: ci,ei,cj,ej
integer,intent(in)::ngauss
double precision,intent(in)::dist
double precision,intent(inout)::ssout
w=zero



distsq=dist*dist
do  it=1,NGAUSS
!do  jt=1,NGAUSS
r=ei(it)
t=ej(1)
rpt=r+t
t1=pi/rpt
t1=t1*t1*t1
t1=sqrt(t1)
t2=-r*t/rpt
t2=t2*distsq
t2=exp(t2)
w  =  w  +   t1*t2*ci(it)*cj(1)

r=ei(it)
t=ej(2)
rpt=r+t
t1=pi/rpt
t1=t1*t1*t1
t1=sqrt(t1)
t2=-r*t/rpt
t2=t2*distsq
t2=exp(t2)
w  =  w  +   t1*t2*ci(it)*cj(2)

r=ei(it)
t=ej(3)
rpt=r+t
t1=pi/rpt
t1=t1*t1*t1
t1=sqrt(t1)
t2=-r*t/rpt
t2=t2*distsq
t2=exp(t2)
w  =  w  +   t1*t2*ci(it)*cj(3)

r=ei(it)
t=ej(4)
rpt=r+t
t1=pi/rpt
t1=t1*t1*t1
t1=sqrt(t1)
t2=-r*t/rpt
t2=t2*distsq
t2=exp(t2)
w  =  w  +   t1*t2*ci(it)*cj(4)

r=ei(it)
t=ej(5)
rpt=r+t
t1=pi/rpt
t1=t1*t1*t1
t1=sqrt(t1)
t2=-r*t/rpt
t2=t2*distsq
t2=exp(t2)
w  =  w  +   t1*t2*ci(it)*cj(5)

r=ei(it)
t=ej(6)
rpt=r+t
t1=pi/rpt
t1=t1*t1*t1
t1=sqrt(t1)
t2=-r*t/rpt
t2=t2*distsq
t2=exp(t2)
w  =  w  +   t1*t2*ci(it)*cj(6)


!end do
end do
ssout=w
end subroutine ssover

