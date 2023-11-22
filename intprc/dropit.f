      SUBROUTINE DROPIT(iunit,vec1,ivec,length,label)
      implicit double precision (a-h,o-z)
      dimension vec1(length),ivec(length)
      write(iunit)vec1,ivec,label
      return
      end
