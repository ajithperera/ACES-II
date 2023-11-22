      integer function asv_strlen (maxlen, s)
c     
c A poor-man Fortran77 reimplementation of the C STRLEN function. Counts
c the amount of characters in the string S INCLUDING (unlike the real
c STRLEN!) the terminating NULL character. MAXLEN specifies the maximum
c number of characters possible, so we do not die if a string has no
c terminating character.
c
c This function is all we need to emulate C-style string handling within
c f77.
c
c Igor Schweigert, Feb 2004
c $Id: asv_strlen.FPP,v 1.2 2008/06/03 15:54:41 taube Exp $
c 
      character*(*) s
      character*1 asv_null
      integer maxlen
c
      asv_null=achar(0)
      do asv_strlen = 1, maxlen
         if (s (asv_strlen:asv_strlen).eq.asv_null) return
      enddo
      return
      end
      
      
            
         
      
