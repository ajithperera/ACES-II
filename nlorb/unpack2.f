      subroutine unpack2(ab,a,b) 

      integer ab, a, b
        
      ioff(ab)=(ab*(ab-1))/2

      a = 1+(-1+aint(dsqrt(8.d0*ab+0.999d0)))/2
      b = ab-ioff(a)
 
      return
      end

