      subroutine asv_get_token (
     &     ncache, icache, cache, accumulate, maxlen,
     &     term, skip, comm,
     &     stat, last, token)
c
c This routine reads a token from the file UNIT. ACCUMULATE specifies if
c characters are actually written into the token and if they are MAXLEN
c specifies the maximum number of characters a token can have. A token
c is a series of characters that is terminated by one of the symbols
c specified by string TERM. Characters specified by string SKIP are not
c written to the token. Neither are characters between the two comment
c separators specified by string COMM. The status of the last reading
c operation is given by STAT, the last character read (which is always
c one of the terminating characters) is return in LAST, and the token is
c passed back through TOKEN.
c
c All the strings are terminated by the NULL character. Also the token
c is concatenated with NULL characters up to MAXLEN.
c
c For example, to read a keyword from ACES2 ZMAT-style file, one may
c specify TERM as '=,\n(\0', SKIP as ' \t', and COMM as '#\n'. Note that
c it is up to the caller, to check the terminating character and
c complain if it is not as expected to be ('=' in the example.)
c
c Note that TERM and SKIP have to be at list one character long and COMM
c at least two characters long.
c
c Note that FGETC is not portable in Fortran77, so instead we get
c characters from the cache (see also comments in ASV_PARSE_LIST)
c
c Igor Schweigert, Feb 2004
c $Id: asv_get_token.FPP,v 1.2 2008/06/03 15:54:41 taube Exp $
c 
      implicit none
c 
c     Arguments
c     
      logical
     &     accumulate
c 
      integer 
     &     ncache, icache, unit, maxlen, stat
c 
      character*(*)
     &     cache, term, skip, comm, token
      character*1
     &     last
c
c     Local variables
c     
      logical comment
      integer len, n, fgetc
      character*1 asv_null
c 
 5900 format (
     &     t3,'@ASV_GET_TOKEN-F: A token is too long. The limit is',/,
     &     t3,i4,' characters')
c 
      asv_null=achar(0)
      do n = 1, maxlen
         token (n:n) = asv_null
      enddo
c      
      len = 0
      comment = .false.
 100  if (icache.lt.ncache) then
         icache = icache + 1
         last (1:1) = cache (icache:icache) 
      else
         stat = -1
      endif
      if (stat.ne.0) then
         return
      elseif (comment) then
         if (last (1:1).eq.comm (2:2)) comment = .false.
         goto 100
      elseif (last (1:1).eq.comm (1:1)) then
         comment = .true.
         goto 100
      else
         if (last (1:1) .eq. asv_null) goto 100
         do n = 1, 256
            if (skip (n:n) .eq. asv_null) goto 110
            if (last (1:1) .eq. skip (n:n)) goto 100
         enddo
 110     do n = 1, 256
            if (term (n:n) .eq. asv_null) goto 120
            if (last (1:1) .eq. term (n:n)) return
         enddo
 120     if (accumulate) then
            if (len.eq.maxlen) then
               write (6, 5900) maxlen
               stat = 1
               return
            endif
            len = len + 1
            token (len:len) = last (1:1)
         endif
         goto 100
      endif
c     
      return
      end

