 subroutine two_electron
  use constants
  use indices
  use control
  use scratch_array
  use tables
 implicit double precision (a-h,o-z)
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
   subroutine three_core(angular)
double precision,intent(inout)::angular
end subroutine three_core

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

ir=0
jc=0
ir2=0
jc2=0
enuc=zero
twoe=zero
h=zero
allocate(ifirst(numat))
allocate(ilast(numat))
allocate(ifirst2(numat))
allocate(ilast2(numat))

if(METHOD=='MNDOD')then
allocate(r2cent(45,45))
call bdata4(r2cent)
end if



  do i=1,numat-1

           ni=nbas(species(i))
           istart=ir+1           
           iend=istart+nbas(species(i))-1
           ifirst(i)=istart
           ilast(i)=iend
           ir=nbas(species(i))+ir  
           jc=nbas(species(i))+jc
           jend=jc

           
           istart2=ir2+1
           iend2=istart2+pairs(species(i))-1
           ir2=pairs(species(i))+ir2
           jc2=pairs(species(i))+jc2
           jend2=jc2
           ifirst2(i)=istart2
           ilast2(i)=iend2
i1=istart2+1
i2=istart2+2
i3=istart2+3
i4=istart2+4
i5=istart2+5
i6=istart2+6
i7=istart2+7
i8=istart2+8
i9=istart2+9
! put in diagonal two electron pieces
    nhere=pairs(species(i))
        if(nhere>10)then
       allocate(ww(pairs(species(i)),pairs(species(i))))
       ww=zero
       call wstore(ww,species(i),0,1)
       call rearrange(WW,pairs(species(i)),pairs(species(i)))

     ii=0
     do m=istart2,istart2+9 ! labels sp block
     ii=ii+1
     jj=10
     do n=istart2+10,iend2 ! labels d block
     jj=jj+1
     twoe(m+offset2(n))=ww(ii,jj)
     end do
     end do

     ii=10
     do m=istart2+10,iend2 ! labels d block rows
     ii=ii+1
     jj=ii-1
     do n=m,iend2 ! labels d block
     jj=jj+1
     twoe(m+offset2(n))=ww(ii,jj)
     end do
     end do



    



       deallocate(ww)
       end if
    if(nhere==1)then
       twoe(istart2+offset2(istart2))=gss(species(i))

    elseif(nhere.ge.10)then
  
! ssss
twoe(istart2+offset2(istart2))=gss(species(i))
! ss|pp
    twoe(offset2(i4)+istart2)=gsp(species(i))
    twoe(offset2(i7)+istart2)=gsp(species(i))
    twoe(offset2(i9)+istart2)=gsp(species(i))
! psps

    twoe(i1+offset2(i1))=hsp(species(i))
    twoe(i2+offset2(i2))=hsp(species(i))
    twoe(i3+offset2(i3))=hsp(species(i))
! pppp
    twoe(i4+offset2(i4))=gpp(species(i))
    twoe(i7+offset2(i7))=gpp(species(i))
    twoe(i9+offset2(i9))=gpp(species(i))
! ppp'p'

   twoe(i4+offset2(i7))=gp2(species(i))
    twoe(i4+offset2(i9))=gp2(species(i))
    twoe(i7+offset2(i9))=gp2(species(i))



! pp'pp'
    gppgp2=half*(gpp(species(i))-gp2(species(i)))
    twoe(i5+offset2(i5))=gppgp2
    twoe(i6+offset2(i6))=gppgp2
    twoe(i8+offset2(i8))=gppgp2
    end if


  do j=i+1,numat

rsq=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
   r=dsqrt(rsq)
   if(r>1025.0d0)then
   go to 90
   end if
   r=r/autoang
   rau=r
   rsq=r*r
   nj=nbas(species(j))


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
 call local(twoelcl,isp,jsp,d1(species(i)),d2(species(i)),d1(species(j)),d2(species(j)) &
,p0(species(i)),p0(species(j)),p1(species(i)),p1(species(j)),p2(species(i)),p2(species(j)),rsq,r)
!  the 22 local integrals are stored in twoelcl in the order given in local.f90
!  now we need to compute the rotation matrix from local to molecular coordinates
! compute   r,theta,phi, and orientation (quadrant) of two atoms
 r=r*autoang
 x1=x(j)-x(i)
 y1=y(j)-y(i)
 z1=z(j)-z(i)
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
!         print*,'so phi is',phi*180.0d0/pi

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
call overlap(i,j,isp,jsp,rau,overlcl)

call rotate(twoemol,theta,phi,twoelcl,overlcl,overmol)

! do the d orbital integrals if necessary)
if( ni>4  .or. nj>4)then
! pass the coordinates in a 2 column array
coor(1,1)=x(i)
coor(2,1)=y(i)
coor(3,1)=z(i)
coor(1,2)=x(j)
coor(2,2)=y(j)
coor(3,2)=z(j)
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

po(1,2)=p0(species(j))
po(2,2)=p1(species(j))
po(3,2)=p2(species(j))
po(4,2)=p3(species(j))
po(5,2)=p4(species(j))
po(6,2)=p5(species(j))
po(7,2)=p6(species(j))
po(8,2)=p7(species(j))
dd(1,2)=zero
dd(2,2)=d1(species(j))
dd(3,2)=d2(species(j))
dd(4,2)=d3(species(j))
dd(5,2)=d4(species(j))
dd(6,2)=d5(species(j))
 call rotmat(1,2,ni,nj,2,coor,rau,yy)
 allocate(ww(pairs(species(i)),pairs(species(j))))
ww=zero
 CALL REPPD (2,1,rau,WW,pairs(species(j)),pairs(species(i)),1,dd,po)
 CALL ROTD  (WW,YY,pairs(species(j)),pairs(species(i)),r2cent)
call rearrange(WW,pairs(species(i)),pairs(species(j)))


if(allocated(dlocal))deallocate(dlocal)
allocate(dlocal(ni,nj))
call doverlap(i,j,rau,dlocal)

end if









! compute the core-core repulsion for this pair which requires two electron integral (ss|ss)
call repel(t1,r,i,j,twoemol(1,1))
enuc=enuc+t1
! put the overlap chunk into the overlap vector
allocate(smatrix(ni,nj))
smatrix(1:isp,1:jsp)=overmol(1:isp,1:jsp)




itop=isp*(isp+1)/2
jtop=jsp*(jsp+1)/2
allocate(w2(pairs(species(i)),pairs(species(j))))
W2(1:itop,1:jtop)=twoemol(1:itop,1:jtop)
if(method=='MNDOD')then
  if(nj>4)then
   w2(1:itop,jtop+1:pairs(species(j)))=ww(1:itop,jtop+1:pairs(species(j)))
   smatrix(1:isp,jsp+1:nj)=dlocal(1:isp,jsp+1:nj)
   end if
  if(ni>4)then
  w2(itop+1:pairs(species(i)),1:jtop)=ww(itop+1:pairs(species(i)),1:jtop)
  smatrix(isp+1:ni,1:jsp)=dlocal(isp+1:ni,1:jsp)
  end if
  if(ni>4 .and. nj >4)then
  w2(itop+1:pairs(species(i)),jtop+1:pairs(species(j)))=ww(itop+1:pairs(species(i)),jtop+1:pairs(species(j)))
smatrix(isp+1:ni,jsp+1:nj)=dlocal(isp+1:ni,jsp+1:nj)
  end if
if(allocated(ww))deallocate(ww)
end if

      jstart=jend+1     
            
      jend=jstart+nbas(species(j))-1
      ii=0              
      do m=istart,iend
       ii=ii+1
       jj=0
     do n=jstart,jend
      jj=jj+1
       S(m+offset1(n))=smatrix(ii,jj)
       if(dabs(S(m+offset1(n)))<1D-7)S(m+offset1(n))=zero
end do 
end do
deallocate(smatrix)


! put the two electron chunk into the matrix
    jstart2=jend2+1
    jend2=jstart2+pairs(species(j))-1
    ii=0

    do m=istart2,iend2
    ii=ii+1
    jj=0
    do n=jstart2,jend2
    jj=jj+1
    twoe(m+offset2(n))=w2(ii,jj)
end do
end do

deallocate(w2)





! compute two center hcore elements while we are here
nj=nbas(species(j))
if(ni==1)then
       if(nj==1)then
          H(offset1(jstart)+istart)=half*S(offset1(jstart)+istart)*(betas(species(i))+betas(species(j)))
      else
          H(offset1(jstart)+istart)=half*S(offset1(jstart)+istart)*(betas(species(i))+betas(species(j)))
          H(offset1(jstart+1)+istart)=half*S(offset1(jstart+1)+istart)*(betas(species(i))+betap(species(j)))
          H(offset1(jstart+2)+istart)=half*S(offset1(jstart+2)+istart)*(betas(species(i))+betap(species(j)))
          H(offset1(jstart+3)+istart)=half*S(offset1(jstart+3)+istart)*(betas(species(i))+betap(species(j)))


           ! s(i) with d(j)
               if(nj>4)then
                  do kk=1,5
     H(offset1(jstart+3+kk)+istart)=half*S(offset1(jstart+3+kk)+istart)*(betas(species(i))+betad(species(j)))
                  end do
               end if
               
        end if

else
   if(nj==1)then
           H(offset1(jstart)+istart)=half*S(offset1(jstart)+istart)*(betas(species(i))+betas(species(j)))
           H(offset1(jstart)+istart+1)=half*S(offset1(jstart)+istart+1)*(betap(species(i))+betas(species(j)))
           H(offset1(jstart)+istart+2)=half*S(offset1(jstart)+istart+2)*(betap(species(i))+betas(species(j)))
           H(offset1(jstart)+istart+3)=half*S(offset1(jstart)+istart+3)*(betap(species(i))+betas(species(j)))
                 if(ni>4)then
                  do kk=1,5

    H(offset1(jstart)+istart+3+kk)=half*S(offset1(jstart)+istart+3+kk)*(betas(species(j))+betad(species(i)))
                 end do
                 end if
             
        


    else

       H(offset1(jstart)+istart)=half*S(offset1(jstart)+istart)*(betas(species(i))+betas(species(j)))
                if(nj>4)then  
                  do kk=1,5

  H(offset1(jstart+3+kk)+istart)=half*S(offset1(jstart+3+kk)+istart)*(betas(species(i))+betad(species(j)))
                  end do
               end if
                if(ni>4)then  
                  do kk=1,5

       H(offset1(jstart)+istart+3+kk)=half*S(offset1(jstart)+istart+3+kk)*(betas(species(j))+betad(species(i)))
                  end do
               end if


             do l=1,3
 H(offset1(jstart)+istart+l)=half*S(offset1(jstart)+istart+l)*(betas(species(j))+betap(species(i)))
 H(offset1(jstart+l)+istart)=half*S(offset1(jstart+l)+istart)*(betap(species(j))+betas(species(i)))

                 do ll=1,3

H(offset1(jstart+l)+istart+ll)=half*S(offset1(jstart+l)+istart+ll)*(betap(species(i))+betap(species(j)))
             end do
                if(ni>4)then
                do ll=1,5

H(offset1(jstart+l)+istart+3+ll)=half*S(offset1(jstart+l)+istart+3+ll)*(betap(species(j))+betad(species(i)))
                end do
                end if

                if(nj>4)then
                do ll=1,5

H(offset1(jstart+3+ll)+istart+l)=half*S(offset1(jstart+3+ll)+istart+l)*(betad(species(j))+betap(species(i)))
                end do
                end if

             
             end do
            ! do the d(i) d(j) block here
     
             if(ni>4.and.nj>4)then
                do kk=4,8
                do jj=4,8

H(offset1(jstart+kk)+istart+jj)=half*S(offset1(jstart+kk)+istart+jj)*(betad(species(j))+betad(species(i)))
              end do
              end do
              end if




     end if
end if










90 continue

 end do
 end do




! do the one that was excluded from the loop above
           

           istart=ir+1           
           iend=istart+nbas(species(numat))-1
           ifirst(numat)=istart
           ilast(numat)=iend
           istart2=ir2+1
           iend2=istart2+pairs(species(numat))-1
           ifirst2(numat)=istart2
           ilast2(numat)=iend2
i1=istart2+1
i2=istart2+2
i3=istart2+3
i4=istart2+4
i5=istart2+5
i6=istart2+6
i7=istart2+7
i8=istart2+8
i9=istart2+9
! put in diagonal two electron pieces
    nhere=pairs(species(numat))
        if(nhere>10)then
       allocate(ww(pairs(species(numat)),pairs(species(numat))))
       ww=zero
       call wstore(ww,species(numat),0,1)
       call rearrange(WW,pairs(species(numat)),pairs(species(numat)))

            ii=0
     do m=istart2,istart2+9 ! labels sp block
     ii=ii+1
     jj=10
     do n=istart2+10,iend2 ! labels d block
     jj=jj+1
     twoe(m+offset2(n))=ww(ii,jj)
     end do
     end do

     ii=10
     do m=istart2+10,iend2 ! labels d block rows
     ii=ii+1
     jj=ii-1
     do n=m,iend2 ! labels d block
     jj=jj+1
     twoe(m+offset2(n))=ww(ii,jj)
     end do
     end do

       deallocate(ww)
       end if
    if(nhere==1)then

twoe(istart2+offset2(istart2))=gss(species(numat))

    elseif(nhere.ge.10)then
! ssss

twoe(istart2+offset2(istart2))=gss(species(numat))
! ss|pp
    twoe(offset2(i4)+istart2)=gsp(species(numat))
    twoe(offset2(i7)+istart2)=gsp(species(numat))
    twoe(offset2(i9)+istart2)=gsp(species(numat))
    

! psps

    twoe(i1+offset2(i1))=hsp(species(numat))
    twoe(i2+offset2(i2))=hsp(species(numat))
    twoe(i3+offset2(i3))=hsp(species(numat))
! pppp

    twoe(i4+offset2(i4))=gpp(species(numat))
    twoe(i7+offset2(i7))=gpp(species(numat))
    twoe(i9+offset2(i9))=gpp(species(numat))
! ppp'p'

    twoe(offset2(i7)+i4)=gp2(species(numat))
    twoe(offset2(i9)+i4)=gp2(species(numat))
    twoe(offset2(i9)+i7)=gp2(species(numat))





! pp'pp'
    gppgp2=half*(gpp(species(numat))-gp2(species(numat)))

    twoe(i5+offset2(i5))=gppgp2
    twoe(i6+offset2(i6))=gppgp2
    twoe(i8+offset2(i8))=gppgp2
    end if



if(method=='MNDOD')deallocate(r2cent)
  end subroutine two_electron

