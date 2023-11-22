      subroutine fuzzy2(x,y,z,r,x1,y1,z1,wn,natoms,rmat)
      implicit double precision (a-h,o-z)
      dimension x(natoms),y(natoms),z(natoms),wn(natoms)
     $,rmat(natoms,natoms)
      parameter(half=0.5d0,one=1.0d0,two=2.0d0)
c iatom is reference so translate to origin
c   so make u(ij) where iatom=i is at origin
        do i=1,natoms
         wn(i)=1.0d0
         ri=(x1-x(i))**2+(y1-y(i))**2+(z1-z(i))**2
      ri=dsqrt(ri)
      do j=1,natoms
         if(j.eq.i)then
         goto 10
         end if
      rj=(x1-x(j))**2+(y1-y(j))**2+(z1-z(j))**2
      rj=dsqrt(rj)
      rij=rmat(j,i)
      uij=(ri-rj)/rij
      f1=1.5d0*uij-half*uij*uij*uij
      f2=1.5d0*f1-half*f1*f1*f1
      f3=1.5d0*f2-half*f2*f2*f2
      s3=half-half*f3
      wn(i)=wn(i)*s3
 10   continue
      end do
      end do


c     now, with all wtemp,compute normalized weight function
      total=0.0d0
      do i=1,natoms
      total=total+wn(i)
      end do
      do i=1,natoms
      wn(i)=wn(i)/total
      end do
 


      return
      end
