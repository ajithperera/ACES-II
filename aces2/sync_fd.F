
c This routine gathers all the data needed for parallel finite differences.
c The goto routine for what should be done is symcor/upd_fd.F

#include "flags.h"

      subroutine sync_fd
      implicit none

#include "mpif.h"

c INTERNAL PARAMETERS
      integer iMemInc
      parameter (iMemInc=1048576)

c INTERNAL VARIABLES
      integer iCore(2), i0, iCrSiz, iMemMin
      double precision dDipXYZ(3), dPolXYZ(3,3)
      integer nDbls, nAtom, nSize, nPoint, iOff, iLast, iTmp
      integer lSend, lRecv, lFree, lGeom, lGrdXYZ, lGrdInt, lIMAP
      character*4 szType
      integer iCount, iRoot, iError
      logical bGeomOpt, bAnGrad

c COMMON BLOCKS
#include "machsp.com"
#include "parallel_aces.com"
#include "flags.com"
#include "flags2.com"

c ----------------------------------------------------------------------

c   o start the JOBARC subsystem
      call aces_ja_init

c   o load system size and number of displacements
      call getrec(1,'JOBARC','NUMPOINT',1,nPoint)
      call getrec(1,'JOBARC','NATOMS',1,nAtom)
      nSize = 3*nAtom

c   o set action flags
      bGeomOpt = iflags2(h_IFLAGS2_geom_opt).ne.0
      bAnGrad  = iflags2(h_IFLAGS2_grad_calc).eq.0

c   o estimate minimum memory requirements
      if (bAnGrad) then
c      - getting the derivatives
         iTmp = iintfp*nSize*3
         iTmp = iTmp + nAtom + iAnd(nAtom,1)
         iTmp = iTmp + iintfp*nSize*(1+nSize)
         iMemMin = iTmp
c      - updating the derivatives
         iTmp = iintfp*(nSize+(nSize+12)*nPoint)
         iMemMin = max(iMemMin,iTmp)
      else
         iMemMin = iintfp*nPoint*2
      end if

c   o allocate memory
      iCrSiz = iflags(h_IFLAGS_mem)
      iCore(1) = 0
      iCore(2) = 0
      do while ((iCore(1).eq.0.and.iCore(2).eq.0).and.
     &          (iCrSiz.gt.iMemMin))
         call aces_malloc(iCrSiz,iCore,i0)
         if (iCore(1).eq.0.and.iCore(2).eq.0) iCrSiz = iCrSiz - iMemInc
      end do
      if (iCore(1).eq.0.and.iCore(2).eq.0) then
         print *, '@SYNC_FD: unable to allocate at least ',
     &            iMemMin,' integers of memory'
         call aces_exit(1)
      end if
c      print *, '@SYNC_FD: Allocated ',iCrSiz,' integers @ iCore(',i0,')'

c   o initialize the iCore free space pointer
      lFree = i0

c   o find last displacement which was calculated
      call getrec(1,'JOBARC','FDCALCTP',nPoint,iCore(lFree))
      iOff = nPoint-1
      do while ((iCore(lFree+iOff).ge.0).and.(iOff.ge.0))
         iOff = iOff - 1
      end do
      iLast = 1+iOff
      if (iLast.eq.0.and.(bGeomOpt.or.irank.ne.0)) then
         print *, '@SYNC_FD: Assertion failed.'
         print *, '          There was no previous calculation.'
         call aces_exit(1)
      end if

c CONCEPT
c    There are 4 records that contain information needed to form numerical
c 1st or 2nd derivatives: ENGPOINT, GRDPOINT, DIPPOINT, and POLPOINT.
c    Numerical gradients require ENGPOINT. Technically, fully numerical
c Hessians also require the REFENERG record, but this point is always done by
c root process, which forms the Hessian, so we don't have to consider it here.
c Analytical gradients, which imply numerical Hessians, require GRD-, DIP-,
c and POL-POINT records.
c    sync_fd will marshall all the records into a continuous array and reduce
c it onto the master process.

      if (bAnGrad) then

c      o determine the scope of symmetry
         call getrec(-1,'JOBARC','DANGERUS',1,iTmp)
         if (iflags(h_IFLAGS_fd_usegroup).eq.0.and.iTmp.eq.0) then
            szType='FULL'
         else
            szType='COMP'
         end if

c      o assign indices
         lGrdInt = lFree
         lFree = lFree + iintfp*nSize
         lGrdXYZ = lFree
         lFree = lFree + iintfp*nSize
         lGeom = lFree
         lFree = lFree + iintfp*nSize
         lIMAP = lFree
         lFree = lFree + nAtom + iand(nAtom,1)
         nDbls = (iCrSiz+i0-lFree)/iintfp

c      o load the original coordinates for the last displacement
         call getrec(1,'JOBARC','FDCOORDS',iintfp*nSize*iLast,
     &               iCore(lGrdInt))
         call c_memmove(iCore(lGeom),
     &                  iCore(lGrdInt+iintfp*nSize*iOff),
     &                  ifltln*nSize)

c      o load and rotate the derivatives
         call getgrd(nAtom,szType,iCore(lGeom),
     &               iCore(lGrdXYZ),iCore(lGrdInt),
     &               dDipXYZ,dPolXYZ,iCore(lIMAP),
     &               iCore(lFree),nDbls,.false.)

c      o reclaim scratch after GRDINT
         lSend = lGrdXYZ
         lRecv = lSend
         iCount = 0

c      o energy gradient
         nDbls = nSize*nPoint
         call getrec(1,'JOBARC','GRDPOINT',iintfp*nDbls,iCore(lRecv))
         call dcopy(nSize,iCore(lGrdInt),1,
     &                    iCore(lRecv+iintfp*nSize*iOff),1)
         lRecv = lRecv + iintfp*nDbls
         iCount = iCount + nDbls

c      o dipole
         nDbls = 3*nPoint
         call getrec(1,'JOBARC','DIPPOINT',iintfp*nDbls,iCore(lRecv))
         call dcopy(3,dDipXYZ,1,iCore(lRecv+iintfp*3*iOff),1)
         lRecv = lRecv + iintfp*nDbls
         iCount = iCount + nDbls

c      o polarizability
         nDbls = 9*nPoint
         call getrec(1,'JOBARC','POLPOINT',iintfp*nDbls,iCore(lRecv))
         call dcopy(9,dPolXYZ,1,iCore(lRecv+iintfp*9*iOff),1)
         lRecv = lRecv + iintfp*nDbls
         iCount = iCount + nDbls

         iRoot = 0
         call MPW_Reduce(iCore(lSend),iCore(lRecv),
     &                   iCount,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                   iRoot,MPI_COMM_WORLD,iError)
         if (iError.ne.0) then
            print *, '@SYNC_FD: MPI reduction failed.'
            print *, '          return code = ',iError
            call aces_exit(1)
         end if

         if (irank.eq.0) then
            call putrec(1,'JOBARC',
     &                  'GRDPOINT',iintfp*nSize*nPoint,iCore(lRecv))
            lRecv = lRecv + iintfp*nSize*nPoint
            call putrec(1,'JOBARC',
     &                  'DIPPOINT',iintfp*3*nPoint,iCore(lRecv))
            lRecv = lRecv + iintfp*3*nPoint
            call putrec(1,'JOBARC',
     &                  'POLPOINT',iintfp*9*nPoint,iCore(lRecv))
         end if

c     else if (.not.bAnGrad) then
      else

         call getrec(1,'JOBARC','ENGPOINT',iintfp*nPoint,iCore(lFree))
         if (iLast.ne.0) then
            if (iflags(h_IFLAGS_excite).eq.0) then
               call getrec(1,'JOBARC','TOTENERG',iintfp,
     &                     iCore(lFree+iintfp*iOff))
            else
               call getrec(1,'JOBARC','TOTENER2',iintfp,
     &                     iCore(lFree+iintfp*iOff))
            end if
         end if
         iCount = nPoint
         iRoot  = 0
         if (bGeomOpt) then
            call MPW_Allreduce(iCore(lFree),iCore(lFree+iintfp*nPoint),
     &                         iCount,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                         MPI_COMM_WORLD,iError)
         else
            call MPW_Reduce(iCore(lFree),iCore(lFree+iintfp*nPoint),
     &                      iCount,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                      iRoot,MPI_COMM_WORLD,iError)
         end if
         if (iError.ne.0) then
            print *, '@SYNC_FD: MPI reduction failed.'
            print *, '          return code = ',iError
            call aces_exit(1)
         end if
         if (bGeomOpt.or.irank.eq.0) then
            call putrec(1,'JOBARC',
     &                  'ENGPOINT',iintfp*nPoint,
     &                  iCore(lFree+iintfp*nPoint))
         end if

c     end if (bAnGrad)
      end if

c   o free memory
      call c_free(iCore)

c   o finalize the JOBARC subsystem and free memory
      call aces_ja_fin

      return
c     end subroutine sync_fd
      end

