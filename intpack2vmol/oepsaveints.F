      subroutine oepsaveints (
     &     luint, label, lnbuf, ibuf, dbuf,
     &     thres, nint, dint)
c
c This routine saves AO integrals on disk in the format similar to VMOL.
c
c To mimic the VMOL integral format first "********" and the label will
c be written to the file LUINT, so LIBR/locate.f can find the
c record. The file should already be opened with FORM = "UNFORMATTED"
c and the recording will start from from wherever the file pointer
c was. So it is up to the calling routine to figure out the position in
c the file.
c
c Following the label, the integrals will be stored in chunks of LNBUF
c doubles, each chunk is preceeded by LNBUF integer indeces, and
c followed by the actual size of the chunck (i.e. LNBUF most of the
c times, but if the last chunk size is less then LNBUF, then the size of
c the chunk.) This method allows to drop all the integrals below a
c certain threshold.
c
c Note that the record size is NINT/LNBUF+1 if NINT%LNBUF == 0 and
c NINT/LNBUF+2 if NINT%LNBUF != 0. The length of the last chunk is
c always LNBUF, also it is followed by an empty chunk with NBUF = -1 to
c be compatable with VMOL.
c
c Note also that the SAVEINTS implementation does not depend on the type
c of the integrals. The integral array is treated as a vector
c (1..NINT), and this vector index is the integral index stored on
c disk. So routines that supply a particular type of integrals to
c SAVEINTS or read these integrals from disk, should know how to pack
c and unpack the actual integral indeces.
c
c Igor Schweigert, Oct 2003
c $Id: oepsaveints.FPP,v 1.1.1.1 2008/05/29 19:35:39 taube Exp $
c
      implicit none
c
c     Arguments
c
      integer
     &     luint, lnbuf, ibuf (lnbuf), nint
c
      character*8 label
c
      double precision
     &     thres, dbuf (lnbuf), dint (nint)
c
c     Local variables
c
      integer
     &     nbuf, n, nsaved
c
c     Messages to the user
c
 5905 format (
     &     t3,'@SAVEINTS-I, ',i12,' integrals are passed to be',/,
     &     t3,'stored on disk.',/)
 5910 format (
     &     t3,'@SAVEINTS-I, ',i12,' integrals have been stored',/,
     &     t3,'on disk.',/)
 5915 format (
     &     t3,'@SAVEINTS-I, The record size is about',i8,'Kb',/,
     &     t3,'assuming 4 bytes per integer and 8 bytes per double.',/)
 5920 format (
     &     t3,'@SAVEINTS-I, ',f5.2'% of the integrals are below the',/,
     &     t3,'threshold of ',d7.1,/)
c
c     Store the integrals on disk: <> Write the "prelabel" and label. <>
c     Initialize the integer and double buffers. <> Loop over integral
c     index. For every integral index: <<>> Check if the integral is
c     beyond the thershold. <<>> If it is add the integral and its index
c     to the buffer. <<>> If buffer is full store it on disk and empty
c     the buffer. <> If there is something in the buffers, set the
c     unused buffer elements to zeros (so nobody gets confused while
c     inspecting the files) and store the entire buffer on disk. <>
c     Write another buffer filled with zeros and NBUF = -1 to indicate
c     the end of the record (see VMOL/TWOEL.F) <> Report the amount of
c     integrals stored and percentage of integrals below the threshold.
c
      write (luint) '********'
      write (luint) label
c 
      nbuf = 0
      call izero (ibuf, lnbuf)
      call zero (dbuf, lnbuf)
      nsaved = 0
c
      do n = 1, nint
         if (dabs (dint (n)) .gt. thres ) then
            nbuf = nbuf + 1
            ibuf (nbuf) = n
            dbuf (nbuf) = dint (n)
            nsaved = nsaved + 1
            if (nbuf .eq. lnbuf) then
               write (luint) dbuf, ibuf, nbuf
               nbuf = 0
            endif
         endif
      enddo
c
      if (nbuf .gt. 0) then
         call izero (ibuf (nbuf+1), lnbuf-nbuf)
         call zero (dbuf (nbuf+1), lnbuf-nbuf)
         write (luint) dbuf, ibuf, nbuf
      endif
c
      call izero (ibuf, lnbuf)
      call zero (dbuf, lnbuf)
      write (luint) dbuf, ibuf, -1
c
      write (6, 5905) nint
      write (6, 5910) nsaved
      write (6, 5920) (1.d0 - dble (nsaved)/dble (nint))*100.d0, thres
      write (6, 5915) (nsaved/lnbuf+2)*(lnbuf*12+4)/1024
c
      return
      end
