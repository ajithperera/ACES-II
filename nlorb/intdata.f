      subroutine intdata(integral,nocc,nvir)
      implicit none

      integer nocc, nvir, iii, aaa

      double precision integral(nvir,nvir,nocc,nocc)

      do iii = 1, nocc

         if (iii.eq.1) then
            open(unit = 100, file = 'integral_o1.out',
     & status = 'new')
               do aaa = 1, nvir
                  write(100,10) aaa,
     & abs(integral(aaa,aaa,iii,iii))
               end do
            close(100)
         end if

         if (iii.eq.2) then
            open(unit = 100, file = 'integral_o2.out',
     & status = 'new')
               do aaa = 1, nvir
                  write(100,10) aaa,
     & abs(integral(aaa,aaa,iii,iii))
               end do
            close(100)
         end if

         if (iii.eq.3) then
            open(unit = 100, file = 'integral_o3.out',
     & status = 'new')
               do aaa = 1, nvir
                  write(100,10) aaa,
     & abs(integral(aaa,aaa,iii,iii))
               end do
            close(100)
         end if

         if (iii.eq.4) then
            open(unit = 100, file = 'integral_o4.out',
     & status = 'new')
               do aaa = 1, nvir
                  write(100,10) aaa,
     & abs(integral(aaa,aaa,iii,iii))
               end do
            close(100)
         end if

         if (iii.eq.5) then
            open(unit = 100, file = 'integral_o5.out',
     & status = 'new')
               do aaa = 1, nvir
                  write(100,10) aaa,
     & abs(integral(aaa,aaa,iii,iii))
               end do
            close(100)
         end if

         if (iii.eq.6) then
            open(unit = 100, file = 'integral_o6.out',
     & status = 'new')
               do aaa = 1, nvir
                  write(100,10) aaa,
     & abs(integral(aaa,aaa,iii,iii))
               end do
            close(100)
         end if

         if (iii.eq.7) then
            open(unit = 100, file = 'integral_o7.out',
     & status = 'new')
               do aaa = 1, nvir
                  write(100,10) aaa,
     & abs(integral(aaa,aaa,iii,iii))
               end do
            close(100)
         end if

         if (iii.eq.8) then
            open(unit = 100, file = 'integral_o8.out',
     & status = 'new')
               do aaa = 1, nvir
                  write(100,10) aaa,
     & abs(integral(aaa,aaa,iii,iii))
               end do
            close(100)
         end if

         if (iii.eq.9) then
            open(unit = 100, file = 'integral_o9.out',
     & status = 'new')
               do aaa = 1, nvir
                  write(100,10) aaa,
     & abs(integral(aaa,aaa,iii,iii))
               end do
            close(100)
         end if

         if (iii.eq.10) then
            open(unit = 100, file = 'integral_o10.out',
     & status = 'new')
               do aaa = 1, nvir
                  write(100,10) aaa,
     & abs(integral(aaa,aaa,iii,iii))
               end do
            close(100)
         end if

         if (iii.eq.11) then
            open(unit = 100, file = 'integral_o11.out',
     & status = 'new')
               do aaa = 1, nvir
                  write(100,10) aaa,
     & abs(integral(aaa,aaa,iii,iii))
               end do
            close(100)
         end if

         if (iii.eq.12) then
            open(unit = 100, file = 'integral_o12.out',
     & status = 'new')
               do aaa = 1, nvir
                  write(100,10) aaa,
     & abs(integral(aaa,aaa,iii,iii))
               end do
            close(100)
         end if

         if (iii.eq.13) then
            open(unit = 100, file = 'integral_o13.out',
     & status = 'new')
               do aaa = 1, nvir
                  write(100,10) aaa,
     & abs(integral(aaa,aaa,iii,iii))
               end do
            close(100)
         end if

         if (iii.eq.14) then
            open(unit = 100, file = 'integral_o14.out',
     & status = 'new')
               do aaa = 1, nvir
                  write(100,10) aaa,
     & abs(integral(aaa,aaa,iii,iii))
               end do
            close(100)
         end if

         if (iii.eq.15) then
            open(unit = 100, file = 'integral_o15.out',
     & status = 'new')
               do aaa = 1, nvir
                  write(100,10) aaa,
     & abs(integral(aaa,aaa,iii,iii))
               end do
            close(100)
         end if

      end do 
 
 10   format(I10,F20.15) 

      return
      end

