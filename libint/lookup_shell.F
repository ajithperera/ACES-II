      subroutine lookup_shell(table, nshells, a, m)
c---------------------------------------------------------------------------
c   Converts basis function coordinates into shell block coordinates,
c   using a lookup table that relates the shell index to the number of
c   functions per shell.
c
c   Arguments:
c   	table		Array containing the ending functions per shell.
c       nshells         Number of shells.
c	a		Function coordinate to be converted.
c	m		Return value for the shell block coordinate 
c			corresponding to a.
c---------------------------------------------------------------------------

      implicit none
      integer nshells
      integer table(nshells)
      integer i, a, m, sum

      do I = 1, nshells
 
         If (table(i) .ge. a) then
            m = i
            return
         endif
      enddo

      m = nshells
      return
      end 
     
