
#ifndef _JOBARC_COM_
#define _JOBARC_COM_

c This common block contains information about the JOBARC file.  The
c JOBARC file has two types of records associated with it: physical
c records and logical records.  Physical records are fixed-length records
c which actually make up the file.  They are used only in the actual I/O
c on the JOBARC file.  Logical records are variable-length records which
c contain the actual data used by the Aces3 programs.  Each logical
c record has a label, size, and location.  The getrec/putrec calls read
c or write a single logical record (whose length is completely unrelated
c to the length of the physical record).
c
c This include file should NEVER be included in a standard module.

c jamaxrec : the maximum number of logical records
c jarecwrd : the length of a physical record in integer words

c jalabel  : the label of each logical record
c jaloc    : the absolute word address of the first element of each
c            logical record in the JOBARC file
c jasize   : the size of each logical record in integer words

c janumrec : the number of physical records in the JOBARC file
c jarecbyt : the length of a physical record in bytes

      character*8 jalabel(jamaxrec)
      integer jaloc(jamaxrec),jasize(jamaxrec),
     &    janumrec,jarecbyt

      common /jobarcc/ jalabel
      common /jobarc/  jaloc,jasize,janumrec,jarecbyt
      save /jobarc/
      save /jobarcc/

#endif /* _JOBARC_COM_ */

