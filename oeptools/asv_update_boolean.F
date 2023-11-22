      subroutine asv_update_boolean (
     &     maxlen, npairs, keylst, vallst,
     &     key, val, bdefault, flag)
c
c This rouinte finds the value (VAL) of a boolean-type keyword (KEY) in
c the list of keyword/value pairs (KEYLST/VALLST) . It checks if value
c is true ("TRUE", "1", "ON") or false ("FALSE", "0", "OFF") and set the
c flag FLAG correspondingly. 
c
c Igor Schweigert, Feb 2004
c $Id: asv_update_boolean.FPP,v 1.4 2008/06/03 15:54:41 taube Exp $
c 
      implicit none
c
c     Arguments
c     
      logical bdefault, flag
c     
      integer
     &     maxlen, npairs
c     
      character*(*)
     &     keylst (npairs), vallst (npairs), key, val
c
c     Local variables
c     
      integer
     &     n, m, asv_strlen     
      character*1 asv_null
      character*10 junk
c
c     The error message
c     
 5900 format (
     &     t3,'@ASV_UPDATE_BOOLEAN-F: Invalid boolean keyword value',/,
     &     t3,'(',a,',',a,'). Valid values are "TRUE", "ON", "1" or',/,
     &     t3,'"FALSE", "OFF", "0".',/)
c 
      asv_null=achar(0)
      call asv_update_keyword (
     &     maxlen, npairs, keylst, vallst,
     &     'default'//asv_null, key, val)
c 
      n = asv_strlen (maxlen, val)
      junk = 'default'//asv_null
 
      if (val (1:n) .eq.'default'//asv_null) then
         flag = bdefault
      elseif (
     &        val (1:n) .eq. 'TRUE'//asv_null .or.
     &        val (1:n) .eq. 'ON'//asv_null   .or.
     &        val (1:n) .eq. '1'//asv_null) then
         flag = .true.
      elseif (
     &        val (1:n) .eq. 'FALSE'//asv_null .or.
     &        val (1:n) .eq. 'OFF'//asv_null .or.
     &        val (1:n) .eq. '0'//asv_null) then
         flag = .false.
      else
         n = asv_strlen (maxlen, key) - 1
         m = asv_strlen (maxlen, val) - 1
         write (6, 5900) key (1:n), val (1:m)
         call errex
      endif
c 
      return
      end
      


