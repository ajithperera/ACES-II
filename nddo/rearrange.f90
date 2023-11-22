subroutine rearrange(WW,ipair,jpair)
implicit double precision(a-h,o-z)
double precision,intent(inout),dimension(:,:)::WW
integer,intent(in)::ipair,jpair
double precision,allocatable,dimension(:,:)::r,s
! rearrange rows first
if(ipair>1)then
allocate(r(1,jpair))
r(1:1,1:jpair)=ww(3:3,1:jpair)
ww(3:3,1:jpair)=ww(4:4,1:jpair)
ww(4:4,1:jpair)=ww(7:7,1:jpair)
ww(7:7,1:jpair)=ww(8:8,1:jpair)
ww(8:8,1:jpair)=ww(6:6,1:jpair)
ww(6:6,1:jpair)=ww(5:5,1:jpair)
ww(5:5,1:jpair)=r
end if

if(jpair>1)then
allocate(s(ipair,1))
s(1:ipair,1:1)=ww(1:ipair,3:3)
ww(1:ipair,3:3)=ww(1:ipair,4:4)
ww(1:ipair,4:4)=ww(1:ipair,7:7)
ww(1:ipair,7:7)=ww(1:ipair,8:8)
ww(1:ipair,8:8)=ww(1:ipair,6:6)
ww(1:ipair,6:6)=ww(1:ipair,5:5)
ww(1:ipair,5:5)=s
end if

return
end


