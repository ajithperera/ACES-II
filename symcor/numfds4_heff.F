      Subroutine numfds4_heff(icore, icrsiz)
C   
c loop over a variety of effective Hamiltonians and t-moments
c At this point we need to reserve space on jobarc for various quantities.
c Also we need to keep everything in core. Run through this twice.
C
      implicit none
C
#include <mrcc_heff.com>
#include <machsp.com>
#include <jodaflags.com>

      character*8 irrnmf(14), irrnmc(14) 
      integer icore(icrsiz)
      integer ioff, icase, nbasx, nstart, nirrep, irreps, itop, 
     &        dirprd, nrow, ndim, icrsiz, ione, ndim2, npoints
      double precision eref
C
      ioff = 1
      ione = 1
c        
      call getrec(20,'JOBARC','PARENERG',iintfp, icore(ioff))         
      ioff = ioff + iintfp
c        
c diabatic orbitals
c
      call getrec(20, 'JOBARC', 'NUMBASIS', 1, nbasx)
      ndim = 0
      call getrec(20, 'JOBARC', 'NDIABOCC', 1, ndim)
      if (ndim .ne. 0) then
      call getrec(20, 'JOBARC', 'DIABCPI0', nbasx*ndim*iintfp,
     $     icore(ioff))
      ioff = ioff + nbasx*ndim*iintfp
      endif
      ndim = 0
      call getrec(20, 'JOBARC', 'NDIABVRT', 1, ndim)
      if (ndim .ne. 0) then
      call getrec(20, 'JOBARC', 'DIABCPA0', nbasx*ndim*iintfp,
     $     icore(ioff))
      ioff = ioff + nbasx*ndim*iintfp
      endif
      
      do icase = 1, ncases
         ndim = 0
         call getrec(-1, 'JOBARC', ndimheff(icase), ione, ndim)
         call getrec(20, 'JOBARC', 'NUMBASIS', 1, nbasx)
         if (icase .eq. 1 .or. icase .eq. 2) then 
            nrow = nbasx
         else
            nrow = nbasx*nbasx
         endif
C
         if (ndim .ne. 0) then
            write(6,*) '@numfds4_heff ', icase, ndimheff(icase), ndim
c           
c process heff part
c     
            ndim2 = ndim*ndim
            call getrec(20,'JOBARC',nameheff(icase),
     $                  ndim2*iintfp, icore(ioff))
            ioff = ioff + ndim2*iintfp
c              
            call getrec(20,'JOBARC',namehdiab(icase),
     $                  nrow*ndim*iintfp, icore(ioff))
#ifdef _DEBUG_LVL0
         write(6,*) ' namehdiab in heff, icase, ndim ',
     $           icase, ndim
         call output(icore(ioff), 1, ndim, 1, ndim, ndim,
     $               ndim, 1)
#endif
            ioff = ioff + nrow*ndim*iintfp
c
            call getrec(20,'JOBARC',nametran(icase),
     $                ndim2*iintfp, icore(ioff))
            ioff = ioff + ndim2*iintfp
c
c now process transition moments
c
            call getrec(20,'JOBARC',nametmom(icase),
     $                  ndim*3*iintfp, icore(ioff))
            ioff = ioff + 3*ndim*iintfp
         endif
      enddo
c
c do same for extra set of diabatization vectors
c
      do icase = 1, ncases
         ndim = 0
         call getrec(-1, 'JOBARC', namenvec(icase), ione, ndim)
         call getrec(20, 'JOBARC', 'NUMBASIS', 1, nbasx)
         if (icase .eq. 1 .or. icase .eq. 2) then 
            nrow = nbasx
         else
            nrow = nbasx*nbasx
         endif
C
         if (ndim .ne. 0) then
            write(6,*) '@numfds4_heff ', icase, namenvec(icase), ndim
            write(6,*) ' extra vec ', namehvec(icase), ioff
            call getrec(20,'JOBARC',namehvec(icase),
     $                  nrow*ndim*iintfp, icore(ioff))
            ioff = ioff + nrow*ndim*iintfp
         endif
      enddo
c
c zero out JOBARC records: This is not necessary as far as I can
c tell. Revisit during Debuggiing. Ajith Perera
c
       itop = ioff
       call aces_ja_truncate('JOBARCST', 1)
c
c now set everything up for the calculation
c
      CALL GETREC(20,'JOBARC','NUMPOINT',IONE,npoints)
#ifdef _DEBUG_LVL0
      write(6,*) ' npoints when initializing JOBARC', npoints
#endif
      itop = ioff
      ioff = 1
C
      call putrec(20,'JOBARC','REFPAR_E', iintfp, icore(ioff))
      call getrec(20,'JOBARC','REFPAR_E', iintfp, eref)
      write(6,*) ' @numfds4_heff, REFPAR_E to JOBARC ', eref
c
      ioff = ioff+iintfp
      call zero(icore(itop), npoints)
      call putrec(20,'JOBARC','PNTENERG', npoints*iintfp, icore(itop))
c
c        
c diabatic orbitals
c
      call getrec(20, 'JOBARC', 'NUMBASIS', 1, nbasx)
      ndim = 0
      call getrec(20, 'JOBARC', 'NDIABOCC', 1, ndim)
      if (ndim .ne. 0) then
      call putrec(20, 'JOBARC', 'DIABCPI0', nbasx*ndim*iintfp,
     $     icore(ioff))
      ioff = ioff + nbasx*ndim*iintfp
      endif
      ndim = 0
      call getrec(20, 'JOBARC', 'NDIABVRT', 1, ndim)
      if (ndim .ne. 0) then
      call putrec(20, 'JOBARC', 'DIABCPA0', nbasx*ndim*iintfp,
     $     icore(ioff))
      ioff = ioff + nbasx*ndim*iintfp
      endif

       do icase = 1, ncases
          ndim = 0
          call getrec(-1, 'JOBARC', ndimheff(icase), ione, ndim)
          call getrec(20, 'JOBARC', 'NUMBASIS', 1, nbasx)
          if (icase .eq. 1 .or. icase .eq. 2) then
              nrow = nbasx
          else
              nrow = nbasx*nbasx
          endif
          if (ndim .ne. 0) then
c
c process heff part
c
             ndim2 = ndim*ndim
             call putrec(20,'JOBARC',refheff(icase),
     $                   ndim2*iintfp, icore(ioff))
             ioff = ioff + ndim2*iintfp
c
             call putrec(20,'JOBARC',namehdiab(icase),
     $                   nrow*ndim*iintfp, icore(ioff))
C
#ifdef _DEBUG_LVL0
            write(6,*) ' namehdiab in heff putrec, nrow*ndim',
     $            nrow*ndim
#endif
C
             ioff = ioff + nrow*ndim*iintfp
c
             call putrec(20,'JOBARC',nametran(icase),
     $                   ndim2*iintfp, icore(ioff))
             ioff = ioff + ndim2*iintfp
c
c also create the pnt records
c
             call zero(icore(itop), ndim2*npoints)
             call putrec(20,'JOBARC',pntheff(icase),
     $                   ndim2*npoints*iintfp, icore(itop))
c
c now process transition moments
c
             call putrec(20,'JOBARC',reftmom(icase),
     $                    ndim*3*iintfp, icore(ioff))
c            call putrec(20,'JOBARC',nametmom(icase),
c    $                    ndim*3*iintfp, icore(ioff))
             ioff = ioff + ndim*3*iintfp
             call zero(icore(itop), ndim*3*npoints)
             call putrec(20,'JOBARC',pnttmom(icase),
     $                    ndim*3*npoints*iintfp, icore(itop))
c
          endif
      enddo
c
c same for extra vectors
c
      do icase = 1, ncases
         ndim = 0
         call getrec(-1, 'JOBARC', namenvec(icase), ione, ndim)
         call getrec(20, 'JOBARC', 'NUMBASIS', 1, nbasx)
         if (icase .eq. 1 .or. icase .eq. 2) then 
            nrow = nbasx
         else
            nrow = nbasx*nbasx
         endif
C
         if (ndim .ne. 0) then
            write(6,*) ' extra vec 2 ', namehvec(icase), ioff
            call putrec(20,'JOBARC',namehvec(icase),
     $                  nrow*ndim*iintfp, icore(ioff))
            ioff = ioff + nrow*ndim*iintfp
         endif
      enddo         
c
      return
      end  

