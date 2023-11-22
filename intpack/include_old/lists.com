
#ifndef _LISTS_COM_
#define _LISTS_COM_

c See sym.com for more information related to this common block.  This
c common block should NEVER be used in any routines other then putlst,
c getlst, and routines called by these two.
c
c If you ever need to change the number of files involved in holding lists,
c change the appropriate parameter statements in sym.com and change the
c logic in the initio routine.
c
c Include this file AFTER sym.com.
c
c moio      : The physical record number on which (sublist,list) begins.
c moiowd    : The word address at which (sublist,list) begins.
c moiosz    : The number of total distributions in (sublist,list).
c moiods    : The size of the individual distribution (sublist,list).
c             This size is in FLOATING POINT WORDS!
c moiofl    : The file on which the list resides.
c totrec    : The total number of records included in each file
c totwrd    : The total number of words used in the last record of each file.
c maxrec    : The maximum number of records in each file (over the complete calculation)
c listfiles : The name of each list file.
c listio    : The unit number of each list file.
c buflenwrd : The length of each I/O cache buffer in integer words
c buflenbyt : The length of each I/O cache buffer in bytes
c numbuffer : The number of I/O cache buffers available
c oldestbuf : The oldest I/O cache buffer

c For bit packing cache information:
c ishfsz    :
c mask1     :
c mask2     :

      integer moio(numsublis,totlist),moiowd(numsublis,totlist),
     &    moiosz(numsublis,totlist),moiods(numsublis,totlist),
     &    moiofl(numsublis,totlist),
     &    totrec(totlistfile),totwrd(totlistfile),
     &    maxrec(totlistfile),
     &    listio(totlistfile),
     &    buflenwrd,buflenbyt,numbuffer,oldestbuf,
     &    ishfsz,mask1,mask2,buflen
      character*80 listfiles(totlistfile)

      common /lists/  moio,moiowd,moiosz,moiods,moiofl,totrec,totwrd,
     &    listio,ishfsz,mask1,mask2,buflenwrd,buflenbyt,numbuffer,
     &    oldestbuf, maxrec,buflen
      common /listsc/ listfiles

c Cache pointers

      integer pbufdir,pbufpos,pbuffil,pbufmod

      common /listp/ pbufdir,pbufpos,pbuffil,pbufmod



c quikget   :

      integer quikget(numsublis,totlist)
      common /where/ quikget

#endif /* _LISTS_COM_ */

