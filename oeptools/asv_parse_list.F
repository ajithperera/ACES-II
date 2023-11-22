      subroutine asv_parse_list (
     &     lucfg, ncache, maxlen, maxkey,
     &     cache, section, nkey, key, value)
c
c This routine parse the configuration file to get key/value pairs from
c a given section. The section is specified by "*TITLE(...)", the record
c delimitors are the comma (,) and newline (\n) characters, and the
c key/value delimitor is the equal (=) character. Stars (*) are not
c permitted in the record. Spaces and tabs are ignored. The key/values
c pairs for the given section are return as arrays of string.
c 
c The conceipt and implementation is as primitive as possible. The keys
c or values cannot be abbreviated, the type of the values cannot be
c specified (its up to the caller to decide what to do with the value
c strings.) For a more universal tool check out the PARSECFG or LIBCFG+
c library (both are distrbuted under the GPL.) The problem with these is
c that they use MALLOC in their low-level routines, so hooking them up
c with ACES2 (and ICORE!) will require considerable modifications.
c
c Note that the entire configuration file is loaded into the cache, and
c then parsed. Moreover, because of the traditional way Fortran77 reads
c formatted files, we have to assume a certain line length, so that the
c file size in the cache is the number of lines times the assumed line
c length. This is somewhat memory ineffecient, but the routine was
c originally written and debugged to work with streamed input. Then it
c turned out that there is no portable implementation of streams for
c Fortran77. So the simpliest workaround is to load the file in a buffer
c and loop over the buffer as we would loop over the stream.
c
c Igor Schweigert, Feb 2004
c $Id: asv_parse_list.FPP,v 1.7 2008/06/03 16:39:40 taube Exp $
c 
      implicit none
c
c     Argument
c     
      integer
     &     lucfg, ncache, maxlen, maxkey, nkey
c     
      character*(*)
     &     cache, section, key (maxkey), value (maxkey)
c
c     Local variables
c     
      integer stat, n, maxline, nchar, ichar, nsection, asv_strlen
      character*1 last,asv_null
c      parameter (asv_null=achar(0))
c
c     Messages to the user.
c   
 5905 format (
     &     t3,'@ASV_PARSE_LIST-F: The configuration file is too big.',/
     &     t3,'The cache size is ',i6,' bytes, the assumed line length'/
     &     t3,'is ',i3,' characters, thus the file cannot exceed ',i5,/
     &     t3,'lines.',/)
 5910 format (
     &     t3,'@ASV_PARSE_LIST-F: Cannot read the configuration file.'/)
 5920 format (
     &     t3,'@ASV_PARSE_LIST-F: Cannot find section ',a,'.',/)
 5930 format (
     &     t3,'@ASV_PARSE_LIST-F: Unexpected end of file while ',/,
     &     t3,'reading section ',a,'.',/)
 5940 format (
     &     t3,'@ASV_PARSE_LIST-F: Wrong record in section ',a,'.',/)
 5950 format (
     &     t3,'@ASV_PARSE_LIST-F: Too many keywords spesified.',/
     &     t3,'The limit is ',i3,'.',/)
 5960 format (
     &     t3,'@ASV_PARSE_LIST-F: A section name is too long.',/
     &     t3,'The limit is ',i2,' characters.',/)
 5963 format (
     &     t3,'@ASV_PARSE_LIST-F: A keyword is too long.',/
     &     t3,'The limit is ',i2,' characters.',/)
 5966 format (
     &     t3,'@ASV_PARSE_LIST-F: A value is too long.',/
     &     t3,'The limit is ',i2,' characters.',/)
c
c     Load the entire configuration file into the cache: <> Assume that
c     one line is no more than 100 or size of the cache characters
c     long. <> Rewind the file pointer. <> Copy the file into the cache
c     line by line.
c     
      asv_null=achar(0)
      maxline = 100
      if (maxline.gt.ncache) maxline = ncache
c     
      rewind (unit=lucfg)
c     
      n = 0
      stat=0
      do while (stat.eq.0)
         read (unit = lucfg, fmt = '(a)', iostat = stat)
     &        cache (n+1:n+maxline-1)
         cache (n+maxline:n+maxline) = '\n'
         n = n + maxline
         if (n+maxline.ge.ncache) then
            write (6, 5905) ncache, maxline, ncache/maxline
            call errex
         endif
         if (stat.gt.0) then
            write (6, 5910)
            call errex
         endif
      enddo
      nchar = n
c
c     Move the pointer to the begining of the section: <> Copy the
c     section name in a scratch array. <> Rewind the file pointer. <>
c     Loop for the next star (*). <> Read the name and check if it is
c     the section needed.
c
      nsection = asv_strlen (maxlen, section)
      key (1) (1:nsection) = section (1:nsection)
c
      ichar = 0
      stat = 0
      do while (value (1) (1:nsection) .ne. key (1) (1:nsection))
         call asv_get_token (
     &        nchar, ichar, cache, .false., maxlen,
     &        '*'//asv_null, asv_null, asv_null//asv_null,
     &        stat, last, value (1))
         if (stat.eq.0) call asv_get_token (
     &        nchar, ichar, cache, .true., maxlen,
     &        '(\n'//asv_null, ' \t,'//asv_null, '#\n',
     &        stat, last, value (1))
         if (stat.gt.0) then
            write (6, 5910)
            call errex
         elseif (stat.lt.0) then
            nkey = 0
            return
         endif
      enddo
c
c     Read the key/value pairs.
c     
      nkey = 0
      do while (last (1:1) .ne.')')
         nkey = nkey + 1
         if (nkey.ge.maxkey) then
            write (6, 5950) maxkey
            call errex
         endif
c     
         call asv_get_token (
     &        nchar, ichar, cache, .true., maxlen,
     &        '=,\n)'//asv_null, ' \t'//asv_null, '#\n',
     &        stat, last, key (nkey))
     
         if (stat.ne.0) then
            if (stat.gt.0) write (6, 5910)
            if (stat.lt.0) write (6, 5930) section (1:nsection-1)
            call errex
         endif
         if (last (1:1) .ne.'=') then
            if (asv_strlen (maxlen, key (nkey)).eq.1) then
               nkey = nkey - 1
               goto 200
            endif
            write (6, 5940) section (1:nsection-1)
            call errex
         endif
c     
         call asv_get_token (
     &        nchar, ichar, cache, .true., maxlen,
     &        '=,\n)'//asv_null, ' \t'//asv_null, '#\n',
     &        stat, last, value (nkey))
         if (stat.ne.0) then
            if (stat.gt.0) write (6, 5910)
            if (stat.lt.0) write (6, 5930) section (1:nsection-1)
            call errex
         endif
         if (last (1:1) .eq.'=') then
            write (6, 5940) section (1:nsection-1)
            call errex
         endif
 200     continue
      enddo
c     
      return
      end
