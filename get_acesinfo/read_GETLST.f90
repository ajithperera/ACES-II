subroutine read_GETLST(string,fileout)
   implicit none
   character(len=*) string
   integer fileout, fileread
   integer i,j,k,Num0
   integer NMO,Nocc,Nvrt

   integer, allocatable :: Num1(:)
   double precision, allocatable :: R1MAT(:)
   integer iMem
   integer, allocatable :: iCore(:)

! COMMON BLOCKS
!#include "flags.com" /* for iflags */
!#include "flags2.com" /* for iflags2 */
!#include "machsp.com" /* for iintfp */
!#include "aces_time.com" /* for timing stats */

   fileread = 21


   if (string.eq.'R1MAT') then
!     iMem = 100*1024*1024
!     allocate(iCore(iMem))
!     call aces_io_init(iCore,1024,iMem,.true.)

      allocate (Num1(2))
      call getrec(fileread,'JOBARC','NOCCORB',2,Num1)  !NOCC
      Nocc = Num1(1) 
      call getrec(fileread,'JOBARC','NVRTORB',2,Num1)  !NVRT
      Nvrt = Num1(1) 
!     NMO = Nocc(1)+Nvrt(1)
      allocate (R1MAT(Nocc*Nvrt))
      call zero(R1MAT,Nocc*Nvrt)
      call getlst(R1MAT,1,1,1,1,490)
      write (fileout,'(/,a)') 'R1 matrix'
      call print_rmat(R1MAT,Nvrt*Nocc,Nvrt,Nocc,1)
      deallocate (R1MAT)
      deallocate (Num1)

!     call aces_io_fin
!     deallocate(iCore)
   end if

end subroutine read_GETLST
