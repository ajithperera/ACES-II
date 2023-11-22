      Subroutine Move_4cent_integrals(v,va1,va2,vb1,vb2,vc1,vc2,vd1,vd2,
     &                                intblk,a1,a2,b1,b2,c1,c2,d1,d2,
     &                                Erd_index,Erd_scale,Nbasis)
      implicit none
      integer va1,va2,vb1,vb2,vc1,vc2,vd1,vd2
      integer a1,a2,b1,b2,c1,c2,d1,d2
      integer a,b,c,d
      integer drange1, drange2
      integer crange1, crange2
      integer brange1, brange2
      integer arange1, arange2
      integer i,j,k,l
      integer Nbasis 

      double precision v(va1:va2,vb1:vb2,vc1:vc2,vd1:vd2)
      double precision intblk(a1:a2,b1:b2,c1:c2,d1:d2) 
      double precision Erd_scale(Nbasis)
      double precision e1,e2,e3,e4
      double precision scale 
      Integer          Erd_index(Nbasis)

      drange1 = max(vd1, d1)
      drange2 = min(vd2, d2)
      crange1 = max(vc1, c1)
      crange2 = min(vc2, c2)
      brange1 = max(vb1, b1)
      brange2 = min(vb2, b2)
      arange1 = max(va1, a1)
      arange2 = min(va2, a2)

      do d = drange1, drange2
      do c = crange1, crange2
      do b = brange1, brange2
      do a = arange1, arange2

         i = Erd_index(d)
         j = Erd_index(c)
         k = Erd_index(b)
         l = Erd_index(a)
      
         e1 = ERD_scale(d)
         e2 = ERD_scale(c)
         e3 = ERD_scale(b)
         e4 = ERD_scale(a)

         Scale = e1*e2*e4*e4
         v(l,k,j,i) = intblk(a,b,c,d)*scale 

      enddo
      enddo
      enddo
      enddo

      return
      end


