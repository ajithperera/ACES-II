
      SUBROUTINE reorder(inmat,outmat,reord,n1,n2,trans)
      implicit none
C     Input Variables
      integer reord(1),n1,n2
      character*1 trans
      double precision inmat(n1,n2)
C     Output Variables
      double precision outmat(n1,n2)
C     Local Variables
      integer ii,jj
C------------------------------------------------------------
      
      if ((trans.eq.'n').or.(trans.eq.'N')) then
         do ii=1,n1
            do jj=1,n2
               outmat(ii,jj)=inmat(ii,reord(jj))
            end do
         end do
      else if ((trans.eq.'t').or.(trans.eq.'T')) then
         do ii=1,n1
            do jj=1,n2
               outmat(reord(ii),jj)=inmat(ii,jj)
            end do
         end do
      else
         call errex
      endif

      end
      
