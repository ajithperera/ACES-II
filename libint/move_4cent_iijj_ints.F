      Subroutine Move_4cent_iijj_ints(v,va1,va2,vb1,vb2,intblk,a1,a2,
     &                                b1,b2,Erd_index,Erd_scale,Nbasis)
      implicit none
      integer va1,va2,vb1,vb2
      integer a1,a2,b1,b2
      integer a,b,c,d
      integer brange1, brange2
      integer arange1, arange2
      integer i,j,k,l
      integer Nbasis 

      double precision v(va1:va2,vb1:vb2)
      double precision intblk(a1:a2,b1:b2) 
      double precision Erd_scale(Nbasis)
      double precision e1,e2,e3,e4
      double precision scale 
      Integer          Erd_index(Nbasis)

      brange1 = max(vb1, b1)
      brange2 = min(vb2, b2)
      arange1 = max(va1, a1)
      arange2 = min(va2, a2)
      Print*,intblk

      do b = brange1, brange2
      do a = arange1, arange2

         i = Erd_index(b)
         j = Erd_index(a)
      
         e2 = ERD_scale(b)
         e1 = ERD_scale(a)

         Scale = e1*e1*e2*e2

         v(j,i) = intblk(a,b)*scale 

      enddo
      enddo

      return
      end


