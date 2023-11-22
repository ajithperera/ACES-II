      subroutine freezelem(dvva,dvvb)
      implicit none
C     Create state-averaged density matrix for ROHF case.
C     Freeze the virtual orbitals in the beta side which correspond
C     to occupied alpha orbitals, and then average within the remaining
C     subblock
C     Common blocks
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer nstart,nirrep
      common/syminf/nstart,nirrep
C     Input/Output Variables
      double precision dvva(*),dvvb(*)
C     Local Variables
      integer irrep,ii,virt(2),idvva,idvvb
C---------------------------------------------------------------
      idvva=1
      idvvb=1
      do irrep=1,nirrep
         virt(1)=vrt(irrep,1)
         virt(2)=vrt(irrep,2)
         if (virt(1).ne.virt(2)) then
            call freezelem_hlp(dvva(idvva),dvvb(idvvb),virt(1),virt(2))
         else
            call daxpy(virt(1)**2,0.5d0,dvva(idvva),1,dvvb(idvvb),1)
            call dcopy(virt(1)**2,dvvb(idvvb),1,dvva(idvva),1)
         endif
         idvva=idvva+virt(1)**2
         idvvb=idvvb+virt(2)**2
      end do
            
      end

C     SUBROUTINE FREEZELEM_HLP

      subroutine freezelem_hlp(Daa,Dbb,dima,dimb)
c
      implicit none
C     Input Variables
      integer dima,dimb
      double precision Daa(dima,dima),Dbb(dimb,dimb)
C     Local Variables
      integer i,j,offset,k
c
      offset=dimb-dima
c
      do 10 k=1,offset
c - subst. rows and columns of Dbb
         do j=1,dimb
            if (j.eq.k) then
               Dbb(k,j)=1.0d0
            else 
               Dbb(k,j)=0.0d0
               Dbb(j,k)=0.0d0
            end if
         end do
 10   continue
c
c - make : Daa <-  1/2 (Daa + Dbb)	  
c        : Dbb <-  1/2 (Daa + Dbb)	  
c
      do i=1,dima
         do j=1,dima
            Daa(i,j)=0.5d0*(Daa(i,j)+Dbb(i+offset,j+offset))
            Dbb(i+offset,j+offset)=Daa(i,j)
         end do
      end do
      
      return
      end
