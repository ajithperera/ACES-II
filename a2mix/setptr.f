      function setptr(last,max,type,num)

      implicit none

      integer iintln,ifltln,iintfp,ialone,ibitwd
      common /machsp/ iintln,ifltln,iintfp,ialone,ibitwd

c This function returns a pointer to some location in an array (which
c has max elements).

      integer setptr,last,max,num,type,stdout
      integer f_integer,f_real
      f_integer=0
      f_real=1

      stdout=6

      setptr=last
      if (type.eq.f_integer) then
        last=last+num
      else if (type.eq.f_real) then
        last=last+num*iintfp
      else
        write(stdout,1000)
 1000   format(/'ERROR: setptr received invalid data type.')
        stop
      endif

c To fix some minor alignment problems, make sure we always start on
c an even numbered element.
      last=last + mod(last,2)

      if (last.ge.max) call insmem('POTDENS',last,max)

      return
      end
