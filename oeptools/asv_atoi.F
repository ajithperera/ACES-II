      subroutine asv_atoi (maxlen, str, int, err)
c
c A poor-man Fortran77 reimplementation of the strandrard STDLIB's ATOI.
c Converts a NULL-terminated string STR into the integer INT. If STR
c does not look like an integer ERR is non-zero on return. A valid
c string is one that contains no non-digit character (except for spaces
c " " and tabs "\t") after the sign. All the characters before the sign
c are ignored.
c
c Note that the integer is NOT checked to be withing the machine
c dependent limits.
c
c Igor Schweigert, Feb 2004 $Id: asv_atoi.FPP,v 1.2 2008/06/03 15:54:41 taube Exp $
c 
      implicit none
      character*(*) str
      character*(1) char
      character*1 asv_null
      integer maxlen, int, err, i, sign, ichar
c     
      asv_null=achar(0)
      sign = 0
      err = 0
      int = 0
      do i = 1, maxlen
         if (str (i:i) .eq. asv_null) then
            goto 10
         elseif (str (i:i).eq.char (9) .or. str (i:i).eq.char (32)) then
            continue
         elseif (sign.eq.0 .and. str (i:i) .eq. '-') then
            sign = -1
         elseif (sign.eq.0 .and. str (i:i) .eq. '+') then
            sign = 1
         elseif (str (i:i) .ge. '0' .and. str (i:i) .le. '9' ) then
            if (sign.eq.0) sign = 1
            int = int*10 + ichar (str (i:i)) - ichar ('0')
         else
            if (sign.ne.0) then
               err = 1
               return
            endif
         endif
      enddo
 10   int = sign * int
      return
      end
