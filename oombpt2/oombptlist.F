
      subroutine oombptlist(vector,length,new_full_diff,spin,cycle,
     &                      read_or_write)
      implicit none
C Common Blocks
      integer ooUnit,cyclim
      common/oombpt/ooUnit,cyclim
C Input variables
      integer length,spin,cycle
      character*1 new_full_diff,read_or_write
C Output if read_or_write='w',Input if read_or_write='r'
      double precision vector(length)
C Local variables
      integer record,stat
C - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - -- - - - -
      if ((new_full_diff.eq.'N').or.(new_full_diff.eq.'n')) then
        record=2+(spin-1)*(2*cyclim+1)
      else if ((new_full_diff.eq.'F').or.(new_full_diff.eq.'f')) then
        record=(spin-1)*(2*cyclim+1)+2*cycle+1
      else if ((new_full_diff.eq.'D').or.(new_full_diff.eq.'d')) then
        record=(spin-1)*(2*cyclim+1)+2*cycle+2
      else
        write(6,*) 'Vector can only be new, full or diff'
        call errex
      endif
      if ((read_or_write.eq.'R').or.(read_or_write.eq.'r')) then
         read(unit=ooUnit,rec=record,err=666,iostat=stat) vector
      else if ((read_or_write.eq.'W').or.(read_or_write.eq.'w')) then
         write(unit=ooUnit,rec=record,err=666,iostat=stat) vector
      else
        write(6,*) 'Can only read or write'
        call errex
      endif
      if (stat.ne.0) write(6,*) 'STAT NOT EQUAL TO ZERO'
      return
666   call aces_io_error('oombptlist',ooUnit,Stat)
      end
