      subroutine zjkapxia(Amat,Z,Xia,nocc,vrt,aa)
      implicit none
      integer nocc,vrt
      double precision Amat(nocc,nocc,nocc,vrt),Z(nocc,nocc)
      double precision Xia(nocc,vrt),aa(nocc,nocc,nocc,vrt)
      integer i,j,k,a
      
      do i=1,nocc
       do a=1,vrt
        do j=1,nocc
          do k=1,nocc
c           if(j .gt. k) then
             xia(i,a)=xia(i,a)+Z(j,k)*(
c     &                 Amat(j,i,k,a)
c     &                +Amat(k,i,j,a)+
     &                 2.0d0*aa(j,i,k,a)) 
       write(*,*)j,i,k,a,aa(j,i,k,a),z(j,k)
c          end if
            end do
          end do
       write(*,*) i,a,xia(a,i)
        end do
      end do 
                  
c      call kkk(nocc*nocc,z)
      return
      end  
