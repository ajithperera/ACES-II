 subroutine mndod_param(k,ussa,uppa,udda,betasa,betapa,betada,zsa,zpa,zda,gssa,gspa,gppa,    &
 gp2a,hspa,nbasis,alp,gdda,zsn,zpn,zdn,a11,a21,a31,a12,a22,a32,a13,a23,a33)
 implicit double precision (a-h,o-z)
 double precision,intent(inout):: ussa,uppa,udda,betasa,betapa,betada,zsa,zpa,zda  &
 ,gssa,gspa,gppa,gp2a,hspa,alp,gdda,zsn,zpn,zdn,a11,a21,a31,a12,a22,a32,a13,a23,a33
 integer,intent(inout)::nbasis

! hydrogen am1
 if(k==1)then
!      USSA=-13.07332100D0
!       UPPA=0.0d0
!       UDDA= 0.0d0
!       BETASA=-5.62651200D0
!       BETAPA=0.0d0
!       BETADA=0.0d0
!       ZSA =0.96780700d0
!       ZPA =0.0d0
!       ZDA =0.0d0
!       GSSA=14.79420800D0
!       GSPA=0.0
!       GPPA=0.0
!       GP2A=0.0
!       HSPA=0.0
!      GDDA=0.0d0
!       nbasis=1
!       alp=3.35638600

USSA=-13.07332100D0
UPPA=-5.86d0
UDDA= 0.0d0
BETASA=-5.62651200D0
BETAPA=-1.720d0
BETADA=0.0d0
ZSA =0.96780700d0
ZPA =.4d0
ZDA =0.0d0
GSSA=14.79420800D0
GSPA=6.76000000d0
GPPA=4.90000000d0
GP2A=3.90000000d0
HSPA=0.20000000d0
GDDA=0.0d0
nbasis=4
alp=3.35638600d0
a11=     1.12875000D0
a21=     5.09628200D0
a31=             1.53746500D0
a12=           -1.06032900D0
a22=         6.00378800D0
a32=           1.57018900D0
a13=       0.0
a23=         0.0
a33=         0.0

elseif(k==9)then
USSA=-110.4353030D0
UPPA=-105.6850470D0
BETASA=-48.4059390D0
BETAPA=-27.7446600D0
ZSA=4.7085550D0
ZPA=2.4911780D0
GSSA=10.4966670D0
GSPA=16.0736890D0
GPPA=14.8172560D0
GP2A=14.4183930D0
HSPA=0.7277630D0
nbasis=4
alp=3.3589210D0
a11=   0.2420790D0
a21=        4.8000000D0
a31=          0.9300000D0
a12=          0.0036070D0
a22=        4.6000000D0
a32=          1.6600000D0
a13=0.0d0
a23=0.0d0
a33=0.0d0
 



elseif(k==14)then
      USSA=-36.05153000
        UPPA =-27.53569100
        UDDA= -14.67743900
        BETASA=-8.21073420
        BETAPA=-4.88462030
        BETADA=-2.60801150
        ZSA =1.91565460
        ZPA =1.68161130
        ZDA =0.96677166
        GSSA=10.74164704
        GSPA=8.07954333
        GPPA=7.43649667
        GP2A=6.56775126
       HSPA=1.44930367
        nbasis=9
        alp=1.66006930
        zsn=1.52929180
        zpn=0.97628075
        zdn=0.93816441
        a11=  -0.39060000D0
        a21=  6.00005400D0
        a31=  0.63226200D0
        a12= 0.05725900D0
        a22=    6.00718300D0
        a32=     2.01998700D0
        a13=   0.0d0
        a23=   0.0d0
        a33=   0.0d0
 elseif(k==16)then
      USSA=-56.88912800
        UPPA =-47.27474500
        UDDA= -25.09511800
        BETASA=-10.99954500
        BETAPA=-12.21543700
        BETADA=-1.88066950
        ZSA =2.22585050
        ZPA =2.09970560
        ZDA =1.23147250 
          nbasis=9
        alp=2.02305950
        zsn=1.73639140
        zpn=1.12118170
        zdn=1.05084670


elseif(k==8)then   
!    USSA=-86.9930020D0
!       UPPA= -71.8795800D0
!       UDDA= 0.0d0
!       BETASA=-45.2026510D0
!       BETAPA=-24.7525150D0
!       BETADA=0.0d0
!       ZSA =3.7965440D0
!       ZPA =2.3894020D0
!       ZDA =0.0d0
!       GSSA=15.7557600D0
!       GSPA=10.6211600D0
!       GPPA=13.6540160D0
!       GP2A=12.4060950D0
!       HSPA=0.5938830D0
!       GDDA=0.0d0
!       nbasis=4
!       alp=3.21710200

        USSA=-86.9930020D0
        UPPA= -71.8795800D0
        UDDA= -5.0d0
        BETASA=-45.2026510D0
        BETAPA=-24.7525150D0
        BETADA=-1.00000000d0
       ZSA =3.7965440D0
       ZPA =2.3894020D0
       ZDA =1.8d0
       GSSA=15.75576000D0
       GSPA=10.62116000D0
       GPPA=13.65401600D0
       GP2A=12.40609500D0
       HSPA=0.59388300D0
       GDDA=0.0d0
       nbasis=9
       alp=3.21710200
       zsn=.3
       zpn=.1
       zdn=.6
       a11=-1.13112800D0
       a21=6.00247700D0
       a31=1.60731100D0
       a12=1.13789100D0
       a22=5.95051200D0
       a32=1.59839500D0
       a13=0.0d0
       a23=0.0d0
       a33=0.0d0


       
 end if

 
end subroutine mndod_param





 
