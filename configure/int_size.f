
      subroutine int_size(iIntLn)
      integer iIntLn, i
      integer i_a, i_b
      integer*1 is1_a, is1_b
      integer*2 is2_a, is2_b
      integer*4 is4_a, is4_b
      integer*8 is8_a, is8_b
      integer(kind=1) ik1_a, ik1_b
      integer(kind=2) ik2_a, ik2_b
      integer(kind=4) ik4_a, ik4_b
      integer(kind=8) ik8_a, ik8_b

      print *, "==> integer size test"
      print *, ""

      i_a = 1
      i_b = 0
      i = 0
      do while (i_a.ne.i_b)
         i_b = i_a
         i_a = ior(ishft(i_b,1),1)
         i = i + 1
      end do
      print *, 'integer   is ',i/8,' bytes'
      iIntLn = i/8

      is1_a = 1
      is1_b = 0
      i = 0
      do while (is1_a.ne.is1_b)
         is1_b = is1_a
         is1_a = ior(ishft(is1_b,1),1)
         i = i + 1
      end do
      print *, 'integer*1 is ',i/8,' bytes'

      is2_a = 1
      is2_b = 0
      i = 0
      do while (is2_a.ne.is2_b)
         is2_b = is2_a
         is2_a = ior(ishft(is2_b,1),1)
         i = i + 1
      end do
      print *, 'integer*2 is ',i/8,' bytes'

      is4_a = 1
      is4_b = 0
      i = 0
      do while (is4_a.ne.is4_b)
         is4_b = is4_a
         is4_a = ior(ishft(is4_b,1),1)
         i = i + 1
      end do
      print *, 'integer*4 is ',i/8,' bytes'

      is8_a = 1
      is8_b = 0
      i = 0
      do while (is8_a.ne.is8_b)
         is8_b = is8_a
         is8_a = ior(ishft(is8_b,1),1)
         i = i + 1
      end do
      print *, 'integer*8 is ',i/8,' bytes'

      ik1_a = 1
      ik1_b = 0
      i = 0
      do while (ik1_a.ne.ik1_b)
         ik1_b = ik1_a
         ik1_a = ior(ishft(ik1_b,1),1)
         i = i + 1
      end do
      print *, 'integer(kind=1) is ',i/8,' bytes'

      ik2_a = 1
      ik2_b = 0
      i = 0
      do while (ik2_a.ne.ik2_b)
         ik2_b = ik2_a
         ik2_a = ior(ishft(ik2_b,1),1)
         i = i + 1
      end do
      print *, 'integer(kind=2) is ',i/8,' bytes'

      ik4_a = 1
      ik4_b = 0
      i = 0
      do while (ik4_a.ne.ik4_b)
         ik4_b = ik4_a
         ik4_a = ior(ishft(ik4_b,1),1)
         i = i + 1
      end do
      print *, 'integer(kind=4) is ',i/8,' bytes'

      ik8_a = 1
      ik8_b = 0
      i = 0
      do while (ik8_a.ne.ik8_b)
         ik8_b = ik8_a
         ik8_a = ior(ishft(ik8_b,1),1)
         i = i + 1
      end do
      print *, 'integer(kind=8) is ',i/8,' bytes'

      print *, ""

      return
      end

