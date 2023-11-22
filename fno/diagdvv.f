      subroutine diagdvv(dvv,u,occ,scr,maxcor,uhf)
c
c diagonalize the virtual virtual relaxed density matrix
c the beta density (if uhf) begins at dvv(nd2(1)+1) 
c
c===============variable declarations and common blocks======================
      implicit none
C     Common blocks
      integer nocco(2),nvrto(2)
      common/info/nocco,nvrto
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
C     Input variables
      integer uhf,maxcor
      double precision dvv(*)
C     Output variables
      double precision u(*),occ(*)
C     Pre-allocated Local Variables
      double precision scr(maxcor)
C     Local Variables
      integer virt,i000,i010,irrep,ii,iocc,idvv,spin,iu,i020
      character*5 cspin(2)
      data cspin /'Alpha',' Beta'/
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c idvv is the offset in the symmetry packed density matrix
c iocc is the offset in the eigenvalues of the density matrix
c iu is the offset in the eigenvectors of the density matrix
      idvv=1
      iocc=1
      iu=1
      write(6,1000)
      write(6,1001)
      write(6,1000)
      do 10 spin=1,uhf+1
         write(6,1002) cspin(spin)
         write(6,1000)
         do 20 irrep=1,nirrep
            virt=vrt(irrep,spin)
            if(virt.gt.0) then
               i000=1
               i010=i000+virt**2
               i020=i010+virt**2
               if(i020.gt.maxcor)call insmem('diagdvv',i020,maxcor)
               call scopy(virt**2,dvv(idvv),1,scr(i000),1)
               call eig(scr(i000),u(iu),1,virt,-1)
               call scopy(virt,scr(i000),virt+1,occ(iocc),1)
            endif
            write(6,1003) irrep
            write(6,2000)(occ(iocc+ii-1),ii=1,virt)
            idvv=idvv+virt**2
            iu=iu+virt**2
            iocc=iocc+virt
 20      continue         
         write(6,1000)
 10   continue
      return
 1000 format(T3,70('-'))
 1001 format(T6,'Natural orbital occupation numbers for virtual ',
     &   'virtual block')
 1002 format(T29,A5,' spin')
 1003 format(T3,'Irrep ',I2)
 2000 format((6(2X,F9.7)))
      end
