
      subroutine nl_init(nltitle,err,printdef0)
c
c The nl_xxx series of routines implements a fortran namelist.
c
c A namelist starts with an "*" in the first column followed immediately
c by the name of the namelist.  This is followed by any number of fields
c of the form KEY=VAL separated by spaces or newlines.  Any number of
c fields may appear on a single line provided the line is no more than
c 80 characters long.  The end of the namelist is signaled by the end of
c the file or a line with an "*" in the first column.
c
c Keys do not have to be a rigidly defined word.  Rather, they can be
c abbreviated.  A sample keyword might be "CALC*ULATION".  This means
c that the first four characters are required, and any remaining
c characters may be omitted, but if they are given, must match the
c keyword appropriately.  In other words, "CALC", "CALCUL" and
c "CALCULATION" all work but "CALCCC" won't.


c This opens the ZMAT file and searches for a namelist.  If found, it is
c read into the array namelist for parsing.
c
c If printdef is true, the namelist will default to printing (though it
c can be overridden by the PRINT_NL keyword in the namelist itself).
c Otherwise it will default to not printing.

#include <aces.h>
#include <callstack.com>
#include <namelist.com>
#include "./include/flags.com"
#include <io.par>
      character *(*) nltitle
      integer err
      logical printdef, printdef0

      character*(nllinelen) line,title
      character*(nllinelen) readline
      integer linenum,l,strlen
      logical eof,streq

      callstack_curr='NL_INIT'
c
      if (.not. multipoint .or. iprint .ge. 1) then 
        printdef = printdef0
      else
        printdef = .false.
      endif
c
c Open ZMAT file and search for the namelist
      open(unit=zio,file=zfil,form='formatted')
      linenum=0
      nlnumline=0
      title='*'//nltitle
      l=strlen(title)
      eof=.false.
   10 line=readline(linenum,zio,eof)
      if (eof) goto 999
      if (.not.streq(title,line(1:l),.true.)) goto 10

c Once the namelist is found, remove the '*TITLE' from it and if
c anything is left, put it as the first line of the namelist text
      line=line(l+1:strlen(line))
      if (strlen(line).gt.0) then
        call upcase(line)
        nltext(1)=line
        nlnumline=1
      endif

c Read until we are done with the file or find a line starting with '*'
   20 line=readline(linenum,zio,eof)
      if (eof) goto 999
      if (line(1:1).eq.'*') goto 999
      nlnumline=nlnumline+1
      if (nlnumline.gt.nlmaxline) then
        write(stdout,'(a)')
     &      '@NL_INIT-F:  maximum length of namelist exceeded'
        call errex
      end if
      call upcase(line)
      nltext(nlnumline)=line
      goto 20

  999 close(zio)
      call nl_log('print_nl',printdef,prt_nl)
      if (prt_nl) call nl_prttop(nltitle)
      if (nlnumline.eq.0) err=1
      return
      end

c Local Variables: c
c mode: fortran c
c fortran-do-indent: 2 c
c fortran-if-indent: 2 c
c fortran-continuation-indent: 4 c
c fortran-comment-indent-style: nil c
c End:
