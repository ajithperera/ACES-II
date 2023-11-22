subroutine loadparam
use tables
character(6)::a
double precision::value
integer::zz
write(*,*)
write(*,*)
write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++'
write(*,*)'Parameters to be read from external source...'
write(*,*)
write(*,*)'User has requested that the following parameters overwrite the defaults:'
itype=size(ntype,1)
open(unit=10,file='PARAMETERS')
! read parameter file and load values until there are no more lines to be read
do
read(10,*,iostat=io)a,zz,value
if(io<0)exit
write(*,20)a,zz,value

! locate row on which parameter is stored
do i=1,itype
if(ntype(i)==zz)exit
end do
!now assign parameter to the proper row
if(a=='USS')then
uss(i)=value
else if(a=='UPP')then
upp(i)=value
else if(a=='UDD')then
udd(i)=value
else if(a=='ZS')then
zs(i)=value
else if(a=='ZP')then
zp(i)=value
else if(a=='ZD')then
zetad(i)=value
else if(a=='BETAS')then
betas(i)=value
else if(a=='BETAP')then
betap(i)=value
else if(a=='BETAD')then
betad(i)=value
else if(a=='GSS')then
gss(i)=value
else if(a=='GPP')then
gpp(i)=value
else if(a=='GSP')then
gsp(i)=value
else if(a=='GP2')then
gp2(i)=value
else if(a=='HSP')then
hsp(i)=value
else if(a=='ZSN')then
zsone(i)=value
else if(a=='ZPN')then
zpone(i)=value
else if(a=='ZDN')then
zdone(i)=value
else if(a=='ALPHA')then
alpha(i)=value

else if(a=='A1')then
g1(1,i)=value
else if(a=='B1')then
g2(1,i)=value
else if(a=='C1')then
g3(1,i)=value
else if(a=='A2')then
g1(2,i)=value
else if(a=='B2')then
g2(2,i)=value
else if(a=='C2')then
g3(2,i)=value
else if(a=='A3')then
g1(3,i)=value
else if(a=='B3')then
g2(3,i)=value
else if(a=='C3')then
g3(3,i)=value
else if(a=='ATHREE')then
alpha3(i)=value
else if(a=='A1C')then
ac(1,i)=value
else if(a=='B1C')then
bc(1,i)=value
else if(a=='A2C')then
ac(2,i)=value
else if(a=='B2C')then
bc(2,i)=value
else if(a=='A4')then
g1(4,i)=value
else if(a=='B4')then
g2(4,i)=value
else if(a=='C4')then
g3(4,i)=value





end if

end do





20 format(A10,I2,F10.5)
end subroutine loadparam

