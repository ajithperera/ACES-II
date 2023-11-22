module constants
double precision,parameter::autoev=27.21D0,autoang=.529167d0,zero=0.0d0
double precision,parameter ::one=1.0d0,  half=0.5d0,  two=2.0d0 ,three=3.0d0 &
,four=4.0d0  ,  pt25=.25d0,  pt125=.125d0,  pt06=.0625d0
double precision , parameter ::   pi=3.141592653589793d0
! globally accessed integers (# atoms, # basis fns, # one center pairs etc
integer::  numat,num_basis,num_pairs,nelectrons,itype,ninput
double precision::enuc,etotal
! sparkle atom stuff
integer::nsparkle

end module constants
