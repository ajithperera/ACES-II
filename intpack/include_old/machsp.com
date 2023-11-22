
#ifndef _MACHSP_COM_
#define _MACHSP_COM_

c This common block contains information about machine dependancies.

c iintln     : the length of an integer word in bytes
c ifltln     : the length of a floating point word in bytes
c iintfp     : the number of integer words in a floating point word
c            : (equal to ifltln/iintln -- either 1 or 2 typically)
c ialone     : mask for 1/4 word (255 or 65536) used for unpacking integral
c            : labels
c ibitwd     : number of bits in an integer word
c idirwrd    : 0 if direct access file records are in bytes, 1 if in words
c idirfact   : related to idirwrd, divide a record length (in bytes) by
c              this to get it in the correct units (bytes or words)
c ilnbuf     : the length of a buffer to read in chunks of integral files

      integer iintln,ifltln,iintfp,ialone,ibitwd,idirwrd,idirfact,ilnbuf
      common /machsp/ iintln,ifltln,iintfp,ialone,ibitwd,idirwrd,
     &    idirfact,ilnbuf
      save /machsp/

#endif /* _MACHSP_COM_ */

