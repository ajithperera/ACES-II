
c This routine drives the AO-based particle-particle ladder contraction(s)
c using a single-pass algorithm.

#include "parallel_aces.h"

c#define _DEBUG_AOLAD2

      subroutine aolad2(iCore,iCoreDim,iUHF,bTau,irp_x,
     &                  list_ao,list_ao_inc)
      implicit none

      integer iCore(*), iCoreDim, iUHF, irp_x, list_ao, list_ao_inc
      logical bTau

#include "mpif.h"

      integer iUnit, iLnBuf
      integer nSizeT, iSizT2Mix
      integer iSpin, irrep, nMO, nAO, i
      integer irp_ij, max_ij
      integer irp_ab, max_ab
      integer iOff, ndx, iTmp, iTmp1, iTmp2, iTmp3, iTmp4
      integer iRoots(8), iBchOff(8), iBchLen(8)
      integer displs(_PA_MAX_PROCS), recvcounts(_PA_MAX_PROCS), iErr

      integer pFree,
     &        pdBf, piBf,
     &        pIK0, pIL0, pJK0, pJL0, pKI0, pKJ0,
     &        pVal, pSym, pTyp, pWhr, pAOS, pMap,
     &        p_T2, p_dZ
      character*80 szFileName
      integer       iFileName
      logical bInCore

      INTEGER        IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     &               IRPDPDAO(8),IRPDPDAOS(8),ISTART(8,8),ISTARTMO(8,3)
      COMMON /AOSYM/ IAOPOP,   IOFFAO,   IOFFV,     IOFFO,
     &               IRPDPDAO,   IRPDPDAOS,   ISTART,     ISTARTMO
#include "machsp.com" /* for iintfp */
#include "syminf.com" /* for nirrep and dirprd(,) */
#include "sympop.com" /* for irpdpd(,) and isytyp(,) */
#include "sym.com" /* for pop(,) and nt() */
#include "info.com" /* for nocco() and nvrto() */
#include "parallel_aces.com"

c ----------------------------------------------------------------------

c   o assign a unit number to the AO integral file(s)
      iUnit = 10

c   o retrieve the numbers of MOs and AOs
      nMO = nocco(1) + nvrto(1)
      call getrec(1,'JOBARC','NBASTOT',1,nAO)

c   o fill in the AOSYM common block
      call getaoinf(iUHF,irp_x)

c   o set the buffer length
      iLnBuf = 600

c   o calculate the total number of back-transformed T2 amplitudes
      if (iUHF.eq.0) then
         nSizeT = 0
         do irp_ij = 1, nirrep
            irp_ab = dirprd(irp_ij,irp_x)
            nSizeT =   nSizeT
     &               + IRPDPDAO(irp_ab)*irpdpd(irp_ij,isytyp(2,46))
         end do
      else
         nSizeT = 0
         do irp_ij = 1, nirrep
            irp_ab = dirprd(irp_ij,irp_x)
            iTmp =   irpdpd(irp_ij,isytyp(2,44))
     &             + irpdpd(irp_ij,isytyp(2,45))
     &             + irpdpd(irp_ij,isytyp(2,46))
            nSizeT = nSizeT + IRPDPDAO(irp_ab)*iTmp
         end do
      end if

c ----------------------------------------------------------------------

c   o create the memory map
      pFree = 1

c   o allocate the integral and index buffers
      pdBf  = pFree
      piBf  = pdBf + max(NT(1),NT(2),iLnBuf,nAO*nMO)*iintfp
      pFree = piBf + max(NT(1),NT(2),iLnBuf,nAO*nMO*iintfp)

c   o allocate temp arrays for AO pair-indices
      pIK0  = pFree
      pIL0  = pIK0 + iLnBuf*8
      pJK0  = pIL0 + iLnBuf*8
      pJL0  = pJK0 + iLnBuf*8
      pKI0  = pJL0 + iLnBuf*8
      pKJ0  = pKI0 + iLnBuf*8
      pFree = pKJ0 + iLnBuf*8

c   o ???
      pVal  = pFree
      pSym  = pVal + iLnBuf*8*iintfp
      pTyp  = pSym + iLnBuf*8
      pWhr  = pTyp + iLnBuf*8
      pAOS  = pWhr + iLnBuf*8
c      pMap  = pAOS + max(iLnBuf*8,nAO)
      pMap  = pAOS + nAO
      pFree = pMap + nAO*nAO

c   o shift pFree to a double boundary (assuming iCore(1) is aligned)
      pFree = pFree + 1 - iand(pFree,1)

c   o restrict the implementation to incore (holding all T2 and Z arrays)
      iTmp = (iCoreDim+1-pFree)/iintfp
      if (iTmp.lt.nSizeT*2) then
         print *, '@AOLAD2: insufficient memory'
         print *, '         need ',nSizeT*2,' dbls, have ',iTmp
         print *, '         retry using AO_LADDERS=MULTIPASS'
         call aces_exit(1)
      end if
      bInCore = .true.
      p_T2  = pFree
      p_dZ  = p_T2 + nSizeT*iintfp
      pFree = p_dZ + nSizeT*iintfp

c ----------------------------------------------------------------------

c   o initialize iAOSym
      ndx = pAOS
      do irrep = 1, nirrep
         do i = 1, IAOPOP(irrep)
            iCore(ndx) = irrep
            ndx = ndx + 1
         end do
      end do

c   o initialize iMap
      call aosymvec(iCore(pMap),nAO)

c   o load the T2(XX,I<J) amplitudes and transpose them to T2(I<J,XX)
      iOff = 0
      ndx  = p_T2
      do irp_ij = 1, nirrep
         irp_ab = dirprd(irp_ij,irp_x)
         max_ij = irpdpd(irp_ij,isytyp(2,46))
         max_ab = IRPDPDAO(irp_ab)
         call getlst(iCore(p_dZ),1,max_ij,1,irp_ij,list_ao+3)
         call transp(iCore(p_dZ),iCore(ndx),max_ij,max_ab)
         ndx = ndx + iintfp*max_ab*max_ij
         iRoots(irp_ab) = iOff
         iOff = iOff + mod(max_ab,nprocs)
         iOff = mod(iOff,nprocs)
         call paces_batch_stat(irank,nprocs,iRoots(irp_ab),max_ab,
     &                         iBchOff(irp_ab),iBchLen(irp_ab))
      end do
      if (iUHF.ne.0) then
         do iSpin = 2, 1, -1
         do irp_ij = 1, nirrep
            irp_ab = dirprd(irp_ij,irp_x)
            max_ij = irpdpd(irp_ij,isytyp(2,43+iSpin))
            max_ab = IRPDPDAO(irp_ab)
            call getlst(iCore(p_dZ),1,max_ij,1,irp_ij,list_ao+iSpin)
            call transp(iCore(p_dZ),iCore(ndx),max_ij,max_ab)
            ndx = ndx + iintfp*max_ab*max_ij
         end do
         end do
      end if
#ifdef _DEBUG_AOLAD2
      print *, 'iRank = ',irank,'; nProcs = ',nprocs
      print *, 'irpdpdao(*) = ',(irpdpdao(irp_ab),irp_ab=1,nirrep)
      print *, 'iRoots(*)   = ',(iRoots(irp_ab),irp_ab=1,nirrep)
      print *, 'iBchOff(*)  = ',(iBchOff(irp_ab),irp_ab=1,nirrep)
      print *, 'iBchLen(*)  = ',(iBchLen(irp_ab),irp_ab=1,nirrep)
#endif

c   o initialize the T2 increments
      call zero(iCore(p_dZ),nSizeT)

c   o open the integral files and process each

      call gfname('IIII',szFileName,iFileName)
      open(unit=iUnit,file=szFileName(1:iFileName),
     &     form='UNFORMATTED',status='OLD')
      call rdaoijkl2(iCore(p_T2),iCore(p_dZ),iCore(pdBf),iCore(piBf),
     &               iCore(pIK0),iCore(pIL0),iCore(pJK0),iCore(pJL0),
     &               iCore(pKI0),iCore(pKJ0),iCore(pVal),iCore(pSym),
     &               iCore(pTyp),iCore(pWhr),iCore(pAOS),iCore(pMap),
     &               iLnBuf,iUnit,iUHF,nAO,iSizT2Mix,irp_x,
     &               iBchOff,iBchLen)
      close(unit=iUnit,status='KEEP')

      if (nirrep.gt.1) then

      call gfname('IJIJ',szFileName,iFileName)
      open(unit=iUnit,file=szFileName(1:iFileName),
     &     form='UNFORMATTED',status='OLD')
      call rdaoijkl2(iCore(p_T2),iCore(p_dZ),iCore(pdBf),iCore(piBf),
     &               iCore(pIK0),iCore(pIL0),iCore(pJK0),iCore(pJL0),
     &               iCore(pKI0),iCore(pKJ0),iCore(pVal),iCore(pSym),
     &               iCore(pTyp),iCore(pWhr),iCore(pAOS),iCore(pMap),
     &               iLnBuf,iUnit,iUHF,nAO,iSizT2Mix,irp_x,
     &               iBchOff,iBchLen)
      close(unit=iUnit,status='KEEP')

      call gfname('IIJJ',szFileName,iFileName)
      open(unit=iUnit,file=szFileName(1:iFileName),
     &     form='UNFORMATTED',status='OLD')
      call rdaoijkl2(iCore(p_T2),iCore(p_dZ),iCore(pdBf),iCore(piBf),
     &               iCore(pIK0),iCore(pIL0),iCore(pJK0),iCore(pJL0),
     &               iCore(pKI0),iCore(pKJ0),iCore(pVal),iCore(pSym),
     &               iCore(pTyp),iCore(pWhr),iCore(pAOS),iCore(pMap),
     &               iLnBuf,iUnit,iUHF,nAO,iSizT2Mix,irp_x,
     &               iBchOff,iBchLen)
      close(unit=iUnit,status='KEEP')

      if (nirrep.gt.2) then

      call gfname('IJKL',szFileName,iFileName)
      open(unit=iUnit,file=szFileName(1:iFileName),
     &     form='UNFORMATTED',status='OLD')
      call rdaoijkl2(iCore(p_T2),iCore(p_dZ),iCore(pdBf),iCore(piBf),
     &               iCore(pIK0),iCore(pIL0),iCore(pJK0),iCore(pJL0),
     &               iCore(pKI0),iCore(pKJ0),iCore(pVal),iCore(pSym),
     &               iCore(pTyp),iCore(pWhr),iCore(pAOS),iCore(pMap),
     &               iLnBuf,iUnit,iUHF,nAO,iSizT2Mix,irp_x,
     &               iBchOff,iBchLen)
      close(unit=iUnit,status='KEEP')

c     end if (nirrep.gt.2)
      end if

c     end if (nirrep.gt.1)
      end if

      if (bInCore) then
         if (iUHF.eq.0) then
            ndx = p_dZ
            do irp_ij = 1, nirrep
               irp_ab = dirprd(irp_ij,irp_x)
               max_ij = irpdpd(irp_ij,isytyp(2,46))
               max_ab = IRPDPDAO(irp_ab)
               iTmp   = iintfp * max(max_ab,max_ij)
               iTmp1  = pdBf
               iTmp2  = iTmp1 + iTmp
               iTmp3  = iTmp2 + iTmp
               iTmp4  = iTmp3 + iTmp
               if (iTmp4.gt.p_T2) then
                  print *, '@AOLAD2: insufficient memory'
                  call aces_exit(1)
               end if
      do i = 1, nprocs
         call paces_batch_stat(i-1,nprocs,iRoots(irp_ab),max_ab,
     &                         displs(i),recvcounts(i))
         recvcounts(i) = max_ij*recvcounts(i)
         displs(i)     = max_ij*displs(i)
      end do
      call MPI_Allgatherv(
     &         iCore(ndx+iintfp*max_ij*iBchOff(irp_ab)),
     &                                           max_ij*iBchLen(irp_ab),
     &              MPI_DOUBLE_PRECISION,
     &         iCore(ndx),recvcounts,displs,
     &              MPI_DOUBLE_PRECISION,
     &         MPI_COMM_WORLD,iErr)
      if (iErr.ne.0) then
         print *, '@AOLAD2: unable to gather Z(I<J,XX)'
         call aces_exit(iErr)
      end if
               call symrhf3(irp_ab,irp_ij,pop(1,1),IAOPOP,max_ij,
     &                      iCore(ndx),
     &                      iCore(iTmp1),iCore(iTmp2),iCore(iTmp3))
               call transp(iCore(ndx),iCore(p_T2),max_ab,max_ij)
               call putlst(iCore(p_T2),1,max_ij,1,irp_ij,list_ao_inc+3)
               ndx = ndx + iintfp*max_ab*max_ij
            end do
         else
            ndx = p_dZ
            do iSpin = 3, 1, -1
            do irp_ij = 1, nirrep
               irp_ab = dirprd(irp_ij,irp_x)
               max_ij = irpdpd(irp_ij,isytyp(2,43+iSpin))
               max_ab = IRPDPDAO(irp_ab)
      do i = 1, nprocs
         call paces_batch_stat(i-1,nprocs,iRoots(irp_ab),max_ab,
     &                         displs(i),recvcounts(i))
         recvcounts(i) = max_ij*recvcounts(i)
         displs(i)     = max_ij*displs(i)
      end do
      call MPI_Allgatherv(
     &         iCore(ndx+iintfp*max_ij*iBchOff(irp_ab)),
     &                                           max_ij*iBchLen(irp_ab),
     &              MPI_DOUBLE_PRECISION,
     &         iCore(ndx),recvcounts,displs,
     &              MPI_DOUBLE_PRECISION,
     &         MPI_COMM_WORLD,iErr)
      if (iErr.ne.0) then
         print *, '@AOLAD2: unable to gather Z(I<J,XX)'
         call aces_exit(iErr)
      end if
               call transp(iCore(ndx),iCore(p_T2),max_ab,max_ij)
               call putlst(iCore(p_T2),1,max_ij,1,irp_ij,
     &                     list_ao_inc+iSpin)
               ndx = ndx + iintfp*max_ab*max_ij
            end do
            end do
         end if
c     end if (bInCore)
      end if

      return
c     end subroutine aolad2
      end

