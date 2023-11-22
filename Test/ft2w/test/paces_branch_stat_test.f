
      program main
      implicit none

      integer iNode, nNodes, iRoot, nTotal, iSum
      integer iOff(8,8), nEls(8,8)

      do nTotal = 0, 3200, 640
      do nNodes = 1, 8
      print *, 'nNodes =',nNodes,'; nTotal =',nTotal

      do iRoot = 0, nNodes-1
      do iNode = 0, nNodes-1
      call paces_batch_stat(iNode,nNodes,iRoot,nTotal,
     &                      iOff(1+iNode,1+iRoot),
     &                      nEls(1+iNode,1+iRoot))
c      print *, iNode,'/(',iRoot,')',nNodes,': ',1+iOff,'-',iOff+nEls,
c     &         ' of ',nTotal
      end do
      end do

      if (.false.) then
      print *, 'iOff:'
      do iNode = 1, nNodes
         print '(i1,a,8(x,i5))',
     &         iNode,':',(iOff(iNode,iRoot),iRoot=1,nNodes)
      end do
      end if
      if (.true.) then
      print *, 'nEls:'
      do iNode = 1, nNodes
         print '(i1,a,8(x,i5))',
     &         iNode,':',(nEls(iNode,iRoot),iRoot=1,nNodes)
      end do
      end if
      print '(/)'

      do iRoot = 1, nNodes
         iSum = 0
         do iNode = 1, nNodes
            iSum = iSum + nEls(iNode,iRoot)
         end do
         if (iSum.ne.nTotal) then
            print *, 'ERROR: col ',iRoot,iSum,nTotal
         end if
      end do
      do iNode = 1, nNodes
         iSum = 0
         do iRoot = 1, nNodes
            iSum = iSum + nEls(iNode,iRoot)
         end do
         if (iSum.ne.nTotal) then
            print *, 'ERROR: row ',iNode,iSum,nTotal
         end if
      end do

      end do
      end do

      end

