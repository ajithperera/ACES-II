      function mn_det3x3(A)
c
c explicitly calculate 3x3 determinant
c
      double precision mn_det3x3, det3x3, A(3,3)
c
      det3x3 = 0.0d0
      det3x3 = det3x3 + A(1,1) * (A(2,2) * A(3,3) - A(2,3)*A(3,2))
      det3x3 = det3x3 - A(1,2) * (A(2,1) * A(3,3) - A(3,1)*A(2,3))
      det3x3 = det3x3 + A(1,3) * (A(2,1) * A(3,2) - A(2,2)*A(3,1))
c
      mn_det3x3 = det3x3

      return
      end
      
