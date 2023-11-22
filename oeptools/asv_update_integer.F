      subroutine asv_update_integer (
     &     maxlen, npairs, keylst, vallst,
     &     key, val, minval, maxval, idefault, int)
c
c This rouinte finds the value (VAL) of a boolean-type keyword (KEY) in
c the list of keyword/value pairs (KEYLST/VALLST) . It checks if value
c is true ("TRUE", "1", "ON") or false ("FALSE", "0", "OFF") and set the
c flag FLAG correspondingly. 
c
c Igor Schweigert, Feb 2004
c $Id: asv_update_integer.FPP,v 1.3 2008/06/03 15:54:41 taube Exp $
c 
      implicit none
c
c     Arguments
c     
      integer
     &     maxlen, npairs, idefault, maxval, minval, int
c     
      character*(*)
     &     keylst (npairs), vallst (npairs), key, val
c
c     Local variables
c     
      integer
     &     n, m, asv_strlen
      character*1 asv_null
      
c
c     The error message
c     
 5900 format (
     &     t3,'@ASV_UPDATE_INTEGER-F: Invalid integer keyword value',/,
     &     t3,'(',a,',',a,'). A valid value must be no less than ', i3,/,
     &     t3,'and no more than ',i3,'.',/)
 5910 format (
     &     t3,'@ASV_UPDATE_INTEGER-F: Invalid integer keyword value',/,
     &     t3,'(',a,',',a,'). The value is not integer!',/)
c 
      asv_null=achar(0)
      call asv_update_keyword (
     &     maxlen, npairs, keylst, vallst,
     &     'default'//asv_null, key, val)
c 
      n = asv_strlen (maxlen, val)
      if (val (1:n) .eq.'default'//asv_null) then
       
         int = idefault
      else
         call asv_atoi (maxlen, val, int, m)
         if (m.ne.0) then
            write (6,*) key, val
            call errex
         endif
      endif
c 
      if (int.lt.minval .or. int.gt.maxval) then
         n = asv_strlen (maxlen, key) - 1
         m = asv_strlen (maxlen, val) - 1 
         write (6, 5900) key (1:n), val (1:m), minval, maxval
         call errex
      endif
c 
      return
      end
