      subroutine rdradgrd(radgrid)

c This will determine which grid to use for each radial point.  It will
c be called twice.  The first time (determined by numgrid=0 at the start
c of the routine) it reads constants and determine numgrid.  The second
c time (numgrid>0) it actually fills in the grids.

      implicit none
#include <aces.h>
c TODO: needs proper gfname lookup
      integer radgio
      character*(*) radgfil
      parameter (radgio=40)
      parameter (radgfil='RADGRD')
#include <intgrtflags.com>
#include <grid.com>
      integer radgrid(int_numradpts)

      integer i,radpt,type,subtype,err,j,flag
      logical exists
      character*80 line

      call callstack_push('RDRADGRD')

      flag=numgrid

c If a file containing specific radial grids (RADGRD) exists, use those.

      inquire(file=radgfil,exist=exists)
      if (exists) then
        open(unit=radgio,file=radgfil,form='formatted')
        rewind(radgio)
        err=0

c       Read the angular grid associated with each radial point.  The first
c       time this is called (numgrid.eq.0), calculate numgrid.

        read (radgio,'(a)') line
        do 30 i=1,int_numradpts
          read (radgio,9000) radpt,type,subtype
 9000     format(i4,i2,i2)

          if (flag.eq.0) then
            numgrid=0
            do 10 j=1,numgrid
              if (type.eq.gridlist(j,1) .and.
     &            subtype.eq.gridlist(j,2)) then
                gridlist(j,3)=gridlist(j,3)+1
                goto 30
              endif
   10       continue
            numgrid=numgrid+1
            gridlist(numgrid,1)=type
            gridlist(numgrid,2)=subtype
            gridlist(numgrid,3)=1
          else
            do 20 j=1,numgrid
              if (type.eq.gridlist(j,1) .and.
     &            subtype.eq.gridlist(j,2)) then
                radgrid(i)=j
                goto 30
              end if
   20       continue
          end if
   30   continue
        close(unit=radgio)

c If the file does not exist, use a default way of assigning grids.

      else
        if (flag.eq.0) then
          numgrid=1

          if (int_ks) then
            gridlist(1,1)=int_defpottype
            gridlist(1,2)=int_defpotgrid
          else
            gridlist(1,1)=int_defenetype
            gridlist(1,2)=int_defenegrid
         end if

         gridlist(1,3)=int_numradpts
        else
          do 40 i=1,int_numradpts
            radgrid(i)=1
   40     continue
        endif
      endif
      call callstack_pop
      return

  900 write (*,9010)
 9010 format(t3,'@RDRADGRD - error in RADGRD file')
      call errex
      return
      end
