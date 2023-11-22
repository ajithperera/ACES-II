      subroutine int_gen_init(Nshells, Ivshloff, ivAngMom, Nfps,
     &                        Npfps, Ncfps, Ixshells, Ixalpha,
     &                        Ixpcoef, End_nfps, Temp, Ispherical)

      implicit none

      integer i, j, n, iCenter, iShell, iOffset, iTmp,
     &        nFps, ivShlOff, temp, End_nfps, Ispherical,
     &        nContAOs, Ncfps, Npfps, Nshells

      integer nSphrs, nCarts, ivAngMom, Ixshells, Ixalpha, 
     &        Ixpcoef, previous_end

      Dimension nFpS(*), ivShlOff(*), Temp(*), End_nfps(*), 
     &          Npfps(*), Ncfps(*), ivAngMom(*), Ixshells(*),
     &          Ixalpha(*), Ixpcoef(*)

      if (iSpherical.ne.0) then
         iOffset = 0
         do iShell = 1, nShells
            ivShlOff(iShell) = iOffset
            nSphrs = 2*ivAngMom(iShell)+1 
            nFpS(iShell)     = nSphrs*nCFpS(iShell)
            iOffset          = iOffset + nFpS(iShell)
         end do
         nContAOs = iOffset
      else
         iOffset = 0
         do iShell = 1, nShells
            ivShlOff(iShell) = iOffset
            nCarts = (ivAngMom(iShell)+1)*(ivAngMom(iShell)+2)/2 
            nFpS(iShell)     = nCarts*nCFpS(iShell)
            iOffset          = iOffset + nFpS(iShell)
         end do
         nContAOs = iOffset
      end if

      do i = 1, nshells
         ixshells(i) = i
      enddo

      do i = 1, nshells
         temp(i) = nfps(i)
      enddo

      do i = 1, nshells
         nfps(i) = temp(ixshells(i))
      enddo

      do i = 1, nshells
         temp(i) = ixalpha(i)
      enddo

      do i = 1, nshells
         ixalpha(i) = temp(ixshells(i))
      enddo

      do i = 1, nshells
         temp(i) = ixpcoef(i)
      enddo

      do i = 1, nshells
         ixpcoef(i) = temp(ixshells(i))
      enddo

      do i = 1, nshells
         temp(i) = ncfps(i)
      enddo

      do i = 1, nshells
         Ncfps(i) = temp(ixshells(i))
      enddo

      do i = 1, nshells
         temp(i) = npfps(i)
      enddo

      do i = 1, nshells
         npfps(i) = temp(ixshells(i))
      enddo

      previous_end = 0
      do iShell = 1, nShells
         end_nfps(iShell) = previous_end + nfps(iShell)
         previous_end     = end_nfps(iShell)
      enddo

c---------------------------------------------------------------------------
c   Shell offsets are different for Flocke package.
c---------------------------------------------------------------------------

      ivShlOff(1) = 0
      do iShell = 2, nShells
         ivShlOff(iShell) = end_nfps(iShell-1)
      enddo 

      return
c     end subroutine int_gen_init
      end
