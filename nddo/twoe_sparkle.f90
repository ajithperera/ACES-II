 subroutine twoe_sparkle(i,j,wout)
  use constants
  use indices
  use control
  use scratch_array
  use tables
 implicit double precision (a-h,o-z)
! double precision, parameter :: autoev=27.21D0,autoang=.529167d0,zero=0.0d0
double precision,dimension(:,:),intent(inout)::wout
integer,intent(in)::i,j
 double precision,dimension(22) :: twoelcl(22)
 double precision,dimension(10,10)::twoemol
 double precision,dimension(5)::overlcl
 double precision,dimension(4,4)::overmol
 double precision,dimension(3,2)::coor
 double precision,dimension(15,45)::YY
 double precision,allocatable,dimension(:,:)::WW,W2,dlocal,smatrix
 logical,allocatable,dimension(:,:)::r2cent
 double precision,dimension(8,2)::po
 double precision,dimension(6,2)::dd
 
 
!********************************************
!         this section lists the interfaces for all the subroutines

interface 

  subroutine local(twoelcl,ni,nj,d1a,d2a,d1b,d2b,p0a,p0b,p1a,p1b,p2a,p2b,rsq,r) 
  double precision,intent(in) ::d1a,d2a,d1b,d2b,p0a,p0b,p1a,p1b,p2a,p2b,rsq,r
  integer,intent(in) ::ni,nj
  double precision , dimension(22),intent(inout) ::twoelcl
  end subroutine local

  subroutine quadrant(x1,y1,z1,phi,sinth)
  double precision, intent(inout)::phi,x1,y1,z1,sinth
  end subroutine quadrant

  subroutine rotate(twoemol,theta,phi,twoelcl,overlcl,overmol)
  double precision,dimension(:,:),intent(inout) ::twoemol,overmol
  double precision,intent(in)::theta,phi
  double precision,dimension(:),intent(in)::twoelcl,overlcl
  end subroutine rotate

  subroutine overlap(iatom,jatom,iorbs,jorbs,rau,overlcl)
  integer,intent(in)::iatom,jatom,iorbs,jorbs
  double precision,dimension(:),intent(inout)::overlcl
  double precision,intent(in)::rau
  end subroutine overlap

 subroutine repel(nuclear,r,iatom,jatom,ssss)
double precision,intent(inout)::nuclear,r
double precision,intent(in)::ssss
integer,intent(in)::iatom,jatom
 end subroutine repel

subroutine rearrange(WW,ipair,jpair)
double precision,intent(inout),dimension(:,:)::WW
integer,intent(in)::ipair,jpair
end subroutine rearrange

      SUBROUTINE WSTORE (ww,NI,MODE,ip)
            integer,intent(in)::mode,ip,NI
          double precision,dimension(:,:),intent(inout)::WW
      end subroutine wstore

subroutine doverlap(iatom,jatom,rau,dlocal)
double precision,intent(in)::rau
integer,intent(in)::iatom,jatom
double precision,intent(inout),dimension(:,:)::dlocal
end subroutine doverlap

end interface
!  end  section 
!********************************************8   




!***********************************************
!            loop over atoms and compute 2eri's

!initialize S


if(METHOD=='MNDOD')then
allocate(r2cent(45,45))
call bdata4(r2cent)
end if



rsq=(x(i)-xsparkle(j))**2+(y(i)-ysparkle(j))**2+(z(i)-zsparkle(j))**2
   r=dsqrt(rsq)
   if(r>1025.0d0)then
   go to 90
   end if
   r=r/autoang
   rau=r
   rsq=r*r
   nj=1
ni=nbas(species(i))

! compute the 22 nonzero integrals for an sp basis set in local coordinates where 
! atom i is always at the origin and atom j is up the positive z axiz which fixes
! the phase in the local coordinates
twoelcl=zero
! do the sp portion
if(ni>4)then
isp=4
else
isp=ni
end if
if(nj>4)then
jsp=4
else
jsp=nj
end if
 call local(twoelcl,isp,jsp,d1(species(i)),d2(species(i)),zero,zero &
,p0(species(i)),p0sparkle(zeffsp(j)),p1(species(i)),zero,p2(species(i)),zero,rsq,r)
!  the 22 local integrals are stored in twoelcl in the order given in local.f90
!  now we need to compute the rotation matrix from local to molecular coordinates
! compute   r,theta,phi, and orientation (quadrant) of two atoms
 r=r*autoang
 x1=xsparkle(j)-x(i)
 y1=ysparkle(j)-y(i)
 z1=zsparkle(j)-z(i)
 costh=z1/r
if(abs(costh).gt.one)costh=dsign(one,costh)
 theta=dacos(costh)
 sinth=dsin(theta)

! special case if two atoms both lie on the x,y,or z axes. Because of my not so cleverly
! chosen convention for rotation, i have to handle these cases separately.
!print*,'atom ji',j,i
!if on + or - z axis,then

 if(abs(sinth)<1D-10)sinth=zero
 if(abs(sinth)==zero)then  ! there is a problem here with machine precision. 1d-20 effectively turns
! this check off
     phi=zero
     if(z1.gt.zero)then
     theta=zero
     else
     theta=pi
     end if
! if not on z axis then...
else
         cosphi=x1/(r*sinth)
         if(abs(cosphi).gt.one)cosphi=dsign(one,cosphi)
         phi=dacos(cosphi)
         cosphi=dcos(phi)
! if on y axis or in plane of y axis then
         if(abs(cosphi)<1D-10)cosphi=zero
         if(cosphi==zero)then
              if(y1.gt.zero)then
              phi=pi/two
              else
              phi=three*pi/two
              end if
         elseif(y1==zero.and.z1==zero)then
               if(x1<zero)then
                  phi=pi
                  else
                  phi=zero
                  end if
         else
         call quadrant(x1,y1,z1,phi,sinth)
          end if


end if

! now we have r,theta, and phi in local coordinates.  use
! these angles to rotate the integrals counter clockwise about z axis thru phi
! followed by a counter clockwise rotation thru theta about y axis
 
! twoemol will return a 10 x 10 matrix of two electron repulsion integrals between
! atoms i and j in molecular coordinates.  these will then be placed into the two electron
! integral vector
! do the overlaps for this pair as well
call rotate(twoemol,theta,phi,twoelcl,overlcl,overmol)


! do the d orbital integrals if necessary)
if( ni>4  .or. nj>4)then
! pass the coordinates in a 2 column array
coor(1,1)=x(i)
coor(2,1)=y(i)
coor(3,1)=z(i)
coor(1,2)=xsparkle(j)
coor(2,2)=ysparkle(j)
coor(3,2)=zsparkle(j)
dd=zero
po=zero

po(1,1)=p0(species(i))
po(2,1)=p1(species(i))
po(3,1)=p2(species(i))
po(4,1)=p3(species(i))
po(5,1)=p4(species(i))
po(6,1)=p5(species(i))
po(7,1)=p6(species(i))
po(8,1)=p7(species(i))
dd(1,1)=zero
dd(2,1)=d1(species(i))
dd(3,1)=d2(species(i))
dd(4,1)=d3(species(i))
dd(5,1)=d4(species(i))
dd(6,1)=d5(species(i))

po(1,2)=p0sparkle(zeffsp(j))
 call rotmat(1,2,ni,nj,2,coor,rau,yy)
 allocate(ww(pairs(species(i)),1))
ww=zero
 CALL REPPD (2,1,rau,WW,1,pairs(species(i)),1,dd,po)
 CALL ROTD  (WW,YY,1,pairs(species(i)),r2cent)
call rearrange(WW,pairs(species(i)),1)

end if









! compute the core-core repulsion for this pair which requires two electron integral (ss|ss)
call repelsparkle(t1,r,i,j,twoemol(1,1))
enuc=enuc+t1
! put the overlap chunk into the overlap vector



itop=isp*(isp+1)/2
jtop=jsp*(jsp+1)/2
allocate(w2(pairs(species(i)),1))
W2(1:itop,1:jtop)=twoemol(1:itop,1:jtop)
  if(method=='MNDOD')then
  if(ni>4)then
  w2(itop+1:pairs(species(i)),1:jtop)=ww(itop+1:pairs(species(i)),1:jtop)
  end if
  if(allocated(ww))deallocate(ww)
end if

wout=w2
if(allocated(w2))deallocate(w2)




90 continue

if(method=='MNDOD')deallocate(r2cent)
  end subroutine twoe_sparkle

