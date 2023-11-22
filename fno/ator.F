      function ator(str,err)

c This takes a string (which should contain a single valid real) and
c returns the value of that integer.  If an error occurs, err is set
c to 1 and ator returns:
c     1.0   if the string is empty
c     2.0   string contains a sign but no number
c     3.0   integer part contains non-digit characters
c     4.0   decimal part contains non-digit characters
c     5.0   exponent part doesn't start with 'E' or 'D'
c     6.0   exponent part empty
c     7.0   exponent part contains non-digit characters
c     8.0   a decimal point with no integer or decimal part entered

#include <aces.h>
#include <const.par>
      M_REAL ator
      integer err
      character*(*) str

      integer iln,sign,i,j,dec,dlen,elen,d,e,flag,strlen,atoi
      character*80 istr,dstr,estr
      character*1 sp,tab,c

      do 1 i=1,80
        istr(i:i)=char(0)
        dstr(i:i)=char(0)
        estr(i:i)=char(0)
    1 continue

      sp=char(32)
      tab=char(9)
      err=0
      sign=1
c i   : The pointer to str
c j   : The pointer to istr, dstr, or estr
      i=1
c dec : A flag saying a decimal point was entered (so ".e4" is treated
c       as an error)
      dec=0

c Error occurs if str is a null string.
      iln=strlen(str)
      if (iln.eq.0) then
        err=1
        ator=one
        return
      end if

c Leading spaces
   10 if (str(i:i).eq.sp .or. str(i:i).eq.tab) then
        i=i+1
        goto 10
      endif

c A sign
      if (str(i:i).eq.'+' .or. str(i:i).eq.'-') then
        if (str(i:i).eq.'-') sign=-1
        i=i+1
        if (i.gt.iln) then
          ator=two
          err=1
          return
        endif
      end if

c Spaces separating sign from number
   20 if (str(i:i).eq.sp .or. str(i:i).eq.tab) then
        i=i+1
        goto 20
      endif

c The integer section may consist only of digits.
      j=1
   30 c=str(i:i)
      if (c.eq.'.') then
        dec=1
        i=i+1
        goto 40
      endif
      if (c.eq.sp .or. c.eq.tab) goto 50
      if (c.eq.'e' .or. c.eq.'E' .or. c.eq.'d' .or. c.eq.'D')
     &    goto 60
      if (c.lt.'0' .or. c.gt.'9') then
        ator=three
        err=1
        return
      end if
      istr(j:j)=str(i:i)
      j=j+1
      i=i+1
      if (i.gt.iln) goto 90
      goto 30

c The decimal section may consist only of digits, at most 8 
   40 j=1
   41 c=str(i:i)
      if (c.eq.sp .or. c.eq.tab) goto 50
      if (c.eq.'e' .or. c.eq.'E' .or. c.eq.'d' .or. c.eq.'D')
     &    goto 60
      if (c.lt.'0' .or. c.gt.'9') then
        ator=four
        err=1
        return
      end if
      if (j .le. 8) then
        dstr(j:j)=str(i:i)
        j=j+1
      end if
      i=i+1
      if (i.gt.iln) goto 90
      goto 41

c Space before exponent
   50 if (str(i:i).eq.sp .or. str(i:i).eq.tab) then
        i=i+1
        goto 50
      endif
      c=str(i:i)
      if (c.eq.'e' .or. c.eq.'E' .or. c.eq.'d' .or. c.eq.'D')
     &    goto 60
      ator=five
      err=1
      return

c The exponent character
   60 i=i+1
      if (i.gt.iln) then
        ator=six
        err=1
        return
      endif

c Space after exponent character
   70 if (str(i:i).eq.sp .or. str(i:i).eq.tab) then
        i=i+1
        goto 70
      endif

c The exponent section may consist only of digits
      j=1
   80 c=str(i:i)
      if (c.lt.'0' .or. c.gt.'9') then
        ator=seven
        err=1
        return
      end if
      estr(j:j)=str(i:i)
      j=j+1
      i=i+1
      if (i.gt.iln) goto 90
      goto 80

c Parse the strings
   90 iln=strlen(istr)
      dlen=strlen(dstr)
      elen=strlen(estr)
      if (dec.eq.1 .and. iln.eq.0 .and. dlen.eq.0) then
        ator=eight
        return
      endif
      if (iln.ne.0) i=atoi(istr,flag)
      if (dlen.ne.0) d=atoi(dstr,flag)
      if (elen.ne.0) e=atoi(estr,flag)
      if (iln.eq.0 .and. dlen.eq.0) then
        ator=one
      else
        ator=zero
        if (iln.ne.0) ator=real(i)
        if (dlen.ne.0) ator=ator + real(d)/ten**dlen
      endif
      if (elen.ne.0) ator=ator*ten**e
      ator=ator*real(sign)

      return
      end

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End:
