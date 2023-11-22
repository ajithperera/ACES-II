 module tables
! define iratio which equals iintfp from aces
integer::iratio,nirrep
integer,dimension(8)::nbfirr
integer,dimension(:),allocatable::zeff,ntype,species,zstore
integer,dimension(:,:),allocatable::opt,ref
double precision,dimension(:,:),allocatable::g1,g2,g3,q,ac,bc
double precision,dimension(:),allocatable::alpha,g,alpha3

double precision ,dimension(:), allocatable ::x,y,z,uss,upp,betas,betap,zs,zp,gss,gsp,gpp,    &
 gp2,hsp,d1,d2,p0,p1,p2,d3,d4,d5,d0,p3,p4,p5,p6,p7,zetad,udd,betad,gdd,hpp,zsone,zpone,zdone,alps,alpp,alpd &
,bsom1,bpom1,bdom1
! dipole moment factors
double precision ,dimension(:),allocatable::hyfsp,hyfpd
double precision,dimension(:),allocatable ::twoe
!vectorize make this a vector and not matrix
double precision ,dimension(:), allocatable ::S,H

 integer ,dimension(:), allocatable ::nbas,pairs

! d orbital stuff
double precision ,dimension(:), allocatable :: F0DD,F2DD,F4DD,F0SD,G2SD, &
              F0PD,F2PD,G1PD,G3PD, &
    IF0SD,IG2SD
double precision,dimension(:,:),allocatable ::REPD
double precision,dimension(30,30)::B
double precision,dimension(30)::F


! list the effective core charges
double precision:: eff_core(86)=(/ &
1.0,2.0, &
1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0, &
1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0, &
1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0, &
1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0, &
1.,2.,3.,3.,3.,3.,3.,3.,3.,&
                 3.,3.,3.,3.,3.,3.,3.,&
                 3.,4.,5.,6.,7.,8.,9.,10.,11.,2.,&
                 3.,4.,5.,6.,7.,1./)

! list the quantum numbers for spd. Note that Hydrogen and Helium use 1s 2p for polarization functions
!integer:: nqs(18)=(/1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3/)
!integer:: nqp(18)=(/2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3/)
!integer:: nqd(18)=(/3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3/)
integer,dimension(86)::nqs,nqp,nqd
data nqs /2*1,8*2,8*3,18*4,18*5,32*6/
data nqp /2*2,8*2,8*3,18*4,18*5,32*6/
data nqd /30*3,18*4,32*5,5*6,3/
character(2),dimension(103)::periodic=(/'H ','HE',& 
                'LI','BE','B ','C ','N ','O ','F ','NE',& 
                 'NA','MG','AL','SI','P ','S ','CL','AR',& 
                 'K ','CA',& 
                 'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',& 
                           'GA','GE','AS','SE','BR','KR',& 
                 'RB','SR',& 
                 'Y ','ZR','NB','MO','TC','RU','RH','PD','AG','CD',& 
                           'IN','SN','SB','TE','I ','XE',& 
                 'CS','BA',& 
                 'LA','CE','PR','ND','PM','SM','EU','GD',& 
                 'TB','DY','HO','ER','TM','YB','LU' ,&
                      'HF','TA','W ','RE','OS','IR','PT','AU','HG',& 
                           'TL','PB','BI','PO','AT','RN',& 
                 'FR','RA',& 
                 'AC','TH','PA','U ','NP','PU','AM','CM',& 
                 'BK','CF','ES','FM','MD','NO','LR'                /)

double precision,dimension(103)::atmass=(/ 1.007825D+00 , 4.00260D+00 , 7.01600D+00 , &
          9.01218D+00  , 11.00931D+00 , 12.00000D+00 , &
          14.00307D+00 , 15.99491D+00 , 18.99840D+00 , &
          19.99244D+00 , 22.98980D+00 , 23.98504D+00 , &
          26.98153D+00 , 27.97693D+00 , 30.97376D+00 , &
          31.97207D+00 , 34.96885D+00 , 39.96238D+00  , &
          38.96371D+00 , 39.96259D+00 , 44.95591D+00  , &
          47.94795D+00 , 50.94396D+00 , 51.94051D+00  , &
          54.93805D+00 , 55.93494D+00 , 58.93320D+00  , &
          57.93535D+00 , 62.93960D+00 , 63.92915D+00  , &
          68.92558D+00 , 73.92118D+00 , 74.92159D+00  , &
          79.91650D+00 , 78.91830D+00 , 83.91150D+00  , &
          84.91180D+00 , 87.90560D+00 , 88.90580D+00  , &
          89.90470D+00 , 92.90640D+00 , 97.90540D+00  , &
          00.00000D+00 , 101.9043D+00 , 102.9055D+00  , &
          105.9032D+00 , 106.9050D+00 , 113.9036D+00  , &
          114.9041D+00 , 117.9034D+00 , 120.9038D+00  , &
          129.9067D+00 , 126.9044D+00 , 131.9042D+00  , &
          132.9054D+00 , 137.9052D+00 , 138.9063D+00  , &
          139.9054D+00 , 140.9076D+00 , 141.9077D+00  , &
          144.9127D+00 , 151.9197D+00 , 152.9212D+00  , &
          157.9241D+00 , 158.9253D+00 , 163.9292D+00  , &
          164.9303D+00 , 165.9320D+00 , 168.9342D+00  , &
          173.9389D+00 , 174.9408D+00 , 179.9465D+00  , &
          180.9480D+00 , 183.9509D+00 , 186.9557D+00  , &
          191.9615D+00 , 192.9629D+00 , 194.9648D+00  , &
          196.9665D+00 , 201.9706D+00 , 204.9744D+00  , &
          207.9766D+00 , 208.9804D+00 , 208.9824D+00  , &
          209.9875D+00 , 222.0157D+00 , 223.0197D+00  , &
          226.0254D+00 , 227.0277D+00 , 232.0381D+00  , &
          231.0359D+00 , 238.0508D+00 , 237.0482D+00  , &
          244.0642D+00 , 243.0614D+00 , 247.0703D+00  , &
          247.0703D+00 , 251.0796D+00 , 252.0829D+00  , &
          257.0751D+00 , 258.0986D+00 , 259.1009D+00  , &
          260.1053D+00 /)
! sparkle atoms
integer,dimension(4)::spcharge(4)=(/2.0,1.0,-2.0,-1.0/)
double precision,dimension(4)::p0sparkle(4)=(/1.,1.,1.,1./)
double precision ,dimension(:), allocatable ::xsparkle,ysparkle,zsparkle
integer,dimension(:),allocatable::zeffsp
double precision,dimension(4)::alphasp(4)=(/1.5,1.5,1.5,1.5/)




 end module tables
