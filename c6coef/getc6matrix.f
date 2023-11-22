      subroutine getc6matrix(numocca,numoccb,numpoints,omega0,
     & wghts,grdvl,polra,polrb,c6cff,totalc6)
      implicit none

      integer iii, jjj, qqq

      integer numocca, numoccb, numpoints
      
      double precision omega0, pi, totalc6

      double precision wghts(numpoints), grdvl(numpoints)

      double precision polra(numpoints,numocca),
     & polrb(numpoints,numoccb)

      double precision c6cff(numocca,numoccb)

      parameter (
     & pi = 3.14159265D0
     & )

      totalc6 = 0.0D0
      do iii = 1, numocca
         do jjj = 1, numoccb
            do qqq = 1, numpoints
           
               c6cff(iii,jjj) = c6cff(iii,jjj) + 
     & (2.0D0*3.0D0*omega0/pi)*wghts(qqq)*
     & polra(qqq,iii)*polrb(qqq,jjj)/(1+grdvl(qqq))**2

            end do 
            totalc6 = totalc6 + c6cff(iii,jjj)
         end do 
      end do

      return
      end
 
