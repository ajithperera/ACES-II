      subroutine asv_update_string (
     &     maxlen, npairs, keylst, vallst,
     &     key, val, valid, default)
c
c This rouinte updates the value VAL of a string-type keyword (KEY) from
c the list of keyword/value pairs (KEYLST/VALLST). If the keyword is not
c in the list, VAL is assigned the default value DEFAULT. This routine
c also checks if the value is among the valid values (which are passed
c in via a string (VALID) with values separated by commas).
c
c All the string must be terminated by a NULL character ('\0').
c
c Igor Schweigert, Feb 2004
c $Id: asv_update_string.FPP,v 1.2 2008/06/03 15:54:41 taube Exp $
c 
      implicit none
c
c     Arguments
c     
      integer
     &     maxlen, npairs
c     
      character*(*)
     &     keylst (npairs), vallst (npairs), default, valid, key, val
c
c     Local variables
c     
      logical isvalid
c     
      integer
     &     n, m, i, j, asv_strlen
     
      character*1 asv_null
c
c     The error message
c     
 5900 format (
     &     t3,'@ASV_UPDATE_STRING-F: Invalid string keyword value',/,
     &     t3,'(',a,',',a,'). A valid value is one the following:',/,
     &     t3,'(',a,')',/)
c
c     Update the string-type keyword: <> Look up th keyword in the list
c     and copy its value if found. <> Check if the value is valid. <>
c     Complain and die if it is not.
c 
      asv_null=achar(0)
      call asv_update_keyword (
     &     maxlen, npairs, keylst, vallst, default, key,val )
c     
      n = asv_strlen (maxlen, val)
      m = asv_strlen (maxlen, valid)
      j = 1
      isvalid = .true.
      do i = 1, m
         if (valid (i:i) .eq. ',') then
            if (isvalid .and. j.eq.n) goto 100
            isvalid = .true.
            j = 1 
         elseif (valid (i:i) .eq. asv_null) then
            if (isvalid .and. j.eq.n) goto 100
            isvalid = .false.
         elseif (valid (i:i) .ne. val (j:j) .or. j.eq.n) then
            isvalid = .false.
         else
            j = j + 1
         endif
      enddo
c     
 100  if (.not.isvalid) then
         n = asv_strlen (maxlen, val) - 1
         m = asv_strlen (maxlen, valid) - 1
         j = asv_strlen (maxlen, key) - 1
         write (6, 5900) key (1:j), val (1:n), valid (1:m)
         call errex
      endif
c     
      return
      end
