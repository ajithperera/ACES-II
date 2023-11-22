      block data
      implicit double precision (a-h,o-z)
      common/ggbase/nfst(4),nlst(4),nx(20),ny(20),nz(20)
      data nfst/1,2,5,11/,nlst/1,4,10,20/
      data nx/0,0,0,1,0,0,0,1,1,2,0,0,0,0,1,1,1,2,2,3/,
     1     ny/0,0,1,0,0,1,2,0,1,0,0,1,2,3,0,1,2,0,1,0/,
     2     nz/0,1,0,0,2,1,0,1,0,0,3,2,1,0,2,1,0,1,0,0/
      end
