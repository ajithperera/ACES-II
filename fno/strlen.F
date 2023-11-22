      function strlen(str)
c
c The str_xxx routines perform a number of useful operations on strings.
c
c integer functions
c =================================================================
c strlen(str)
c    This returns the length of the string as it appears to a human (i.e.
c    it returns the location of the last non-space character).  0 is returned
c    if the string doesn't contain any non-blank characters.
c findex(string,substring)
c fcindex(string,char)
c fclindex(string,charlist)
c    Return the position of the first occurence of substring, char, or
c    any of the characters in charlist in string respectively.  If no
c    occurence is found, 0 is returned.
c lindex(string,substring)
c lcindex(string,char)
c lclindex(string,charlist)
c    Return the position of the last occurence.
c nindex(string,substring,ind)
c ncindex(string,char,ind)
c nclindex(string,charlist,ind)
c    Returns the position of the next occurence AFTER postion ind.  If
c    ind=0, the search starts at the start of the string.  If ind is
c    greater than the length of the string, 0 is returned.
c pindex(string,substring,ind)
c pcindex(string,char,ind)
c pclindex(string,charlist,ind)
c    Returns the position of the previous occurence BEFORE position.  If
c    ind=0 or ind is greater than the length of the string, the search
c    starts at the end of the string.
c fncindex(string,char)
c fnclindex(string,charlist)
c lncindex(string,char)
c lnclindex(string,charlist)
c nncindex(string,char,index)
c nnclindex(string,charlist,index)
c pncindex(string,char,index)
c pnclindex(string,charlist,index)
c    Return the positions of the occurences that are NOT the character
c    or one of the characters in the character list.
c
c logical functions
c =================================================================
c streq(str1,str2,insensitive)
c    This returns .true. if str1 and str2 are equal.  The comparison can
c    be case insensitive.
c
c subroutines
c =================================================================
c upcase(str)
c    This function converts str to all uppercase.
c delspc (instr, outstr)
c    This strips all the spaces and tabs from instr.  The resulting string
c    is returned as outstr (both instr and outstr may be the same string
c    in the calling program).


#include <aces.h>
      integer strlen
      character *(*) str
      character*1 sp,tab,zero
c      callstack_curr='STRLEN'
      strlen=len(str)
      sp=char(32)
      tab=char(9)
      zero=char(0)
   10 if (strlen.ne.0) then
        if (str(strlen:strlen).eq.sp .or. str(strlen:strlen).eq.tab
     &      .or. str(strlen:strlen).eq.zero) then
          strlen=strlen-1
          goto 10
        end if
      end if
      return
      end

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End:
