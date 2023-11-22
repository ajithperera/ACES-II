      subroutine plunk2(lufil,buf,ibuf,len,nut,nrec,jrec)
      implicit double precision(a-h,o-z)
      dimension buf(len),ibuf(len)
      write(lufil,rec=nrec) buf,ibuf,nut,jrec
      return
      end
