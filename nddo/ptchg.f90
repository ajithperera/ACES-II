subroutine ptchg
use constants
use tables
use control
implicit double precision (a-h,o-z)


interface
subroutine chgsep(L,out,condon,D)
integer,intent(in)::L
double precision,intent(inout)::out
double precision,intent(in)::condon,D
end subroutine chgsep

subroutine factorial(i,total)
implicit double precision(a-h,o-z)
integer,intent(inout)::total
integer,intent(in)::i
end subroutine factorial
end interface



d1=zero
d2=zero
p0=zero
p1=zero
p2=zero
!p2=zero
if(method=='MNDOD')then
d0=zero
d3=zero
d4=zero
d5=zero
p3=zero
p4=zero
p5=zero
p6=zero
p7=zero
end if

!stdoutprint*,''
!stdoutwrite(*,*)'Table of computed charge separations and additive terms used to compute'
!stdoutwrite(*,*)'two electron repulsion integrals.  See (TCA (1977) 46,89-104) for details'
!stdoutwrite(*,*)'Ensure that the computed values are equal to the reference values.  If this'
!stdoutwrite(*,*)'is found to be in error, a finer grid will be necessary in routine limit.f90'
!stdoutwrite(*,*)'------------------------------------------------------------------------------'
do i=1,itype
if(nbas(i)>1)then


! comment out the old code b/c need a more general section to handle 1s-2p basis set cases
! new code follows commented block
!xx = float(nq(ntype(i)))
!rterm=two * xx + 1
!term=rterm * (four * zs(i) * zp(i) ) ** (xx + half)
!term = term / (zs(i) + zp(i)) ** (two*xx + two)
!term = term * .577350269d0
!d1(i)=term

!term= rterm * (two*xx + two) * .05d0
!term=dsqrt(term)
!term=term/zp(i)
!d2(i)=term
!end if
xx = float(nqs(ntype(i)))
yy = float(nqp(ntype(i)))
term1=(  two*zs(i)   )  ** ( xx + half )
term2=(  two*zp(i)   )  ** ( yy + half )
term3=( zs(i)  +  zp(i)  ) ** (xx + yy + two)
term3=one/term3
call factorial(2*nqs(ntype(i)),iterm4)
call factorial(2*nqp(ntype(i)),iterm5)
aa  =  float(iterm4)
bb=float(iterm5)
term6=one/dsqrt(aa*bb)
call factorial(nqs(ntype(i))+nqp(ntype(i))+1,iterm7)
term8 = term1 * term2 * term3 * term6 * float(iterm7)
d1(i)=term8*.577350269d0
hyfsp(i)=d1(i)*5.0832d0

xx = float(nqp(ntype(i)))
yy = float(nqp(ntype(i)))
term1=(  two*zp(i)   )  ** ( xx + half )
term2=(  two*zp(i)   )  ** ( yy + half )
term3=( zp(i)  +  zp(i)  ) ** (xx + yy + three)
term3=one/term3
call factorial(2*nqp(ntype(i)),iterm4)
call factorial(2*nqp(ntype(i)),iterm5)
aa  =  float(iterm4)
bb=float(iterm5)
term6=one/dsqrt(aa*bb)
call factorial(nqp(ntype(i))+nqp(ntype(i))+2,iterm7)
term8 = term1 * term2 * term3 * term6 * float(iterm7)
d2(i)=dsqrt(term8)*0.447213595d0
end if




! 0 = ss  1 = sp  2 = pp 3 = ds  4 = dp 5 =dd
v1=zero
v2=zero
c1=zero
c2=zero
    if(nbas(i)==1)then
    
        condon=gss(i)/autoev
        call chgsep(0,p0(i),condon,zero)
        value=autoev/(two*p0(i))
        c0=condon*autoev
    else if(nbas(i)>1)then
        condon=gss(i)/autoev
        call chgsep(0,p0(i),condon,zero)
        value=autoev/(two*p0(i))
        c0=condon*autoev

        c1=hsp(i)/autoev
        call chgsep(1,p1(i),c1,d1(i))
        call func(p1(i),v1,1,zero,d1(i))
        v1=v1*autoev
        c1=c1*autoev

        c2=half * ( gpp(i) - gp2(i) )
        c2=c2/autoev
        call chgsep(2,p2(i),c2,d2(i))
        call func(p2(i),v2,2,zero,d2(i))
        v2=v2*autoev
        c2=c2*autoev
        if(method=='MNDOD')p6(i)=p0(i)
    end if
! do the d orbital additive and charge separations if necessary
  if(nbas(i)>4)then
! d3 holds the sd charge separation
xx= float(nqs(ntype(i)))
yy = float(nqd(ntype(i)))
term1=(  two*zs(i)   )  ** ( xx + half )
term2=(  two*zetad(i)   )  ** ( yy + half )
term3=( zs(i)  +  zetad(i)  ) ** (xx + yy + three)
term3=one/term3
call factorial(2*nqs(ntype(i)),iterm4)
call factorial(2*nqd(ntype(i)),iterm5)
aa  =  float(iterm4)
bb=float(iterm5)
term6=one/dsqrt(aa*bb)
call factorial(nqs(ntype(i))+nqd(ntype(i))+2,iterm7)
term8 = term1 * term2 * term3 * term6 * float(iterm7)
d3(i)=dsqrt(term8)*.508132748d0

condon= REPD(19,i)
condon=condon/autoev
call chgsep(2,p3(i),condon,d3(i))
call func(p3(i),v3,2,zero,d3(i))
v3=v3*autoev
c3=condon*autoev








xx= float(nqp(ntype(i)))
! pd charge separation
term1=(  two*zp(i)   )  ** ( xx + half )
term2=(  two*zetad(i)   )  ** ( yy + half )
term3=( zp(i)  +  zetad(i)  ) ** (xx + yy + two)
term3=one/term3
call factorial(2*nqp(ntype(i)),iterm4)
call factorial(2*nqd(ntype(i)),iterm5)
aa  =  float(iterm4)
bb=float(iterm5)
term6=one/dsqrt(aa*bb)
call factorial(nqp(ntype(i))+nqd(ntype(i))+1,iterm7)
term8 = term1 * term2 * term3 * term6 * float(iterm7)
d4(i)=term8*.447213595d0
condon=REPD(23,I)-1.8D0*REPD(35,I) 
condon=condon/autoev
call chgsep(1,p4(i),condon,d4(i))
call func(p4(i),v4,1,zero,d4(i))
v4=v4*autoev
c4=condon*autoev
hyfpd(i)=d4(i)*5.0832d0







!dd charge separation
term1=(  two*zetad(i)   )  ** ( yy + half )
term2=(  two*zetad(i)   )  ** ( yy + half )
term3=( zetad(i)  +  zetad(i)  ) ** (yy + yy + three)
term3=one/term3
call factorial(2*nqd(ntype(i)),iterm4)
call factorial(2*nqd(ntype(i)),iterm5)
aa  =  float(iterm4)
bb=float(iterm5)
term6=one/dsqrt(aa*bb)
call factorial(nqd(ntype(i))+nqd(ntype(i))+2,iterm7)
term8 = term1 * term2 * term3 * term6 * float(iterm7)
d5(i)=dsqrt(term8)*.377964473d0
condon=REPD(44,I)-(20.0D0/35.0D0)*REPD(52,I)
condon=condon/autoev
call chgsep(2,p5(i),condon,d5(i))
call func(p5(i),v5,2,zero,d5(i))
v5=v5*autoev
c5=condon*autoev
! do the pp additive term
p6(i)=p0(i)
! do the dd additive term ! see TCA paper 
condon=0.2D0*(REPD(29,I)+TWO*REPD(30,I)+TWO*REPD(31,I))
condon=condon/autoev
call chgsep(0,p7(i),condon,zero)
v7=autoev*half/p7(i)
c7=condon*autoev









end if


     
     



!stdoutwrite(*,*)''
!stdoutwrite(*,*)'Charge separations and additive terms for ',periodic(ntype(i))
!stdoutwrite(*,*)'--------------------------------------------------------------'
!stdoutwrite(*,10)'L','D','P','Calc','Ref'
!stdoutwrite(*,11)0,zero,p0(i),value,c0
!stdoutwrite(*,11)1,d1(i),p1(i),v1,c1
!stdoutwrite(*,11)2,d2(i),p2(i),v2,c2
if(nbas(i)>4)then
!stdoutwrite(*,11)2,d3(i),p3(i),v3,c3
!stdoutwrite(*,11)1,d4(i),p4(i),v4,c4
!stdoutwrite(*,11)2,d5(i),p5(i),v5,c5
!stdoutwrite(*,11)0,zero,p6(i),value,value
!stdoutwrite(*,11)0,zero,p7(i),v7,c7
end if
end do
10 format(A1,A7,A10,A15,A10)
11 format(i1,f10.4,f10.4,f15.10,f15.10)
end subroutine ptchg



