      subroutine dropvrt(dvv,drop,modrop,redrop,reord,evals,nmo,uhf)
      implicit none
C     Common blocks
      integer nocco(2),nvrto(2)
      common/info/nocco,nvrto
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
C     Input Variables
      integer uhf,drop,nmo

C     Input/Output Variables
      double precision dvv(*)

C     Pre-allocated Local Variables
      integer modrop(drop),redrop(drop),reord(nmo)
      double precision evals(nmo)

C     Local Variables
      integer dropocc,ii,jj,idvv,orboff,irrep,spin,virt
      logical trunc
      character*8 scfeval(2)
      data scfeval /'SCFEVLA0','SCFEVLB0'/
C----------------------------------------------------------------------     
      call getrec(20,'JOBARC','MODROPA ',drop,modrop)
      dropocc=0
      do ii=1,drop
         if (modrop(ii).le.nocco(1)) then
            dropocc=dropocc+1
         endif
         if (modrop(ii).gt.nmo) then
            drop=drop-1
         endif
      end do
      if (dropocc.lt.drop) then
         do spin=1,uhf+1         
            call getrec(20,'JOBARC',scfeval(spin),nmo*iintfp,evals)
            do ii=1,nmo
               reord(ii)=ii
            end do
            call piksr2(nmo,evals,reord)
            do ii=dropocc+1,drop
               redrop(ii-dropocc)=reord(modrop(ii))
            end do
            idvv=1+nd2(1)*(spin-1)
            orboff=nocco(spin)
            do irrep=1,nirrep
               virt=vrt(irrep,spin)
               do ii=1,virt
                  trunc=.false.
                  do jj=1,drop-dropocc
                     if (redrop(jj).eq.(ii+orboff)) then
                        trunc=.true.
                     endif
                  end do                  
                  if (trunc) then
                     do jj=1,virt
                        dvv(idvv+(ii-1)*virt+jj-1)=0.0d0
                        dvv(idvv+(jj-1)*virt+ii-1)=0.0d0
                     end do
                  endif
               end do
               idvv=idvv+virt**2
               orboff=orboff+virt
            end do
         end do
      endif
      call putrec(20,'JOBARC','FNOFREEZ',1,dropocc)
      drop=drop-dropocc
      end
