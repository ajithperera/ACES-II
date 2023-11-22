      subroutine asv_update_keyword (
     &     maxlen, npairs, keylst, vallst, default, key, val)
c
c This rouinte finds the value (VAL) of a generic type keyword (KEY) in
c the list of keyword/value pairs (KEYLST/VALLST) . If the keyword is
c not found in the list VAL is assign the default value DEFAULT.
c
c All the string must be terminated by a NULL character ("\0").
c
c Igor Schweigert, Feb 2004
c $Id: asv_update_keyword.FPP,v 1.1.1.1 2008/05/29 19:35:39 taube Exp $
c 
      implicit none
      integer maxlen, npairs, n, m, i, asv_strlen
      character*(*) keylst (npairs), vallst (npairs), default, key, val
c 
      n = asv_strlen (maxlen, default)
      val (1:n) = default (1:n)
c
      if (npairs.eq.0) return
      n = asv_strlen (maxlen, key)
      do i = 1, npairs
         if (keylst (i) (1:n) .eq. key (1:n)) then
            m = asv_strlen (maxlen, vallst (i))
            val (1:m) = vallst (i) (1:m)
            goto 100
         endif
      enddo
 100  continue
c 
      return
      end
