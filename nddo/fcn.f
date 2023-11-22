      subroutine fcn( n, x, f, g ,natoms,q,internal,force,rnorm)
      integer n,i,ncycle,internal(2,n),natoms
      double precision x(n), f, g(n),q(natoms,3),qold,delta
     $,fnew,fold,rnorm
      logical stat1,stat2
      delta=.00001d0
          NCYCLE=NCYCLE+1
          print*,''
         PRINT*,'-----------------------------------------------'
         PRINT*,'GEOMETRY OPTIMIZATION CYCLE NUMBER: ',NCYCLE
         do 510 i=1,n
        q(  internal(1,i), internal(2,i)  ) = x(i)
510     continue
          stat1=.false.
          stat2=.true.
        call scfq(q,natoms,f,stat1,stat2)

          stat1=.true.
          stat2=.false.
         do 400 i=1,n
         qold=q(  internal(1,i), internal(2,i)  )

        q(  internal(1,i), internal(2,i)  ) = qold + delta
           call scfq(q,natoms,fnew,stat1,stat2)

             q(  internal(1,i), internal(2,i)  )= qold - delta
           call scfq(q,natoms,fold,stat1,stat2)
            g(i)=(fnew-fold)/(2.0d0*delta)
           
           q(  internal(1,i), internal(2,i)  ) = qold

 400      continue
         
c
c     calculate ||F|| of current point to monitor optimization
c
         rnorm=0.0D0
         do 80 i=1,n
         rnorm=rnorm+g(i)**2
 80      continue
         rnorm=dsqrt(rnorm)
         print*,''
         print*,'****************************************************'
         print*,'*       AT THE END OF CYCLE NUMBER ',NCYCLE,':      '
         PRINT*,'* GRADIENT NORM OF OPTIMIZED COORDINATES = ',RNORM
         PRINT*,'*            ENERGY = ',F,' eV                      '
          print*,'*****************************************************'
         



      
      

      return
      end

