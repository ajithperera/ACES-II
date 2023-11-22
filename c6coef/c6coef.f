      program c6coef
      implicit none

      integer icore(1), memscr, iiicrtop, iuhf,
     & iintln, ifltln, iintfp, ialone, ibitwd, jjjcrtop

      integer iii00000, iiiwghts, iiifreqs, iiigrdvl, iiipolra,
     & iiipolrb, iiic6cff

      integer jjj00000

      integer numpoints, numocca, numoccb
      
      double precision omega0, totalc6
     
      character*7 system
 
      parameter (
     & memscr   = 200000,
     & numpoints= 16,
     & omega0   = 0.4D0,
     & numocca  = 26,
     & numoccb  = 26,
     & system  = 'DIGLY  '
     & )

      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      common /istart/ iii00000, jjjcrtop 
      common /iending/ jjj00000, iiicrtop

      call aces_init(icore,iii00000,jjjcrtop,iuhf,.true.)
      iiicrtop = jjjcrtop - memscr
      jjj00000 = iiicrtop

      write(*,*)
      write(*,10) numpoints
      write(*,20) omega0
      write(*,30) numocca
      write(*,40) numoccb
      write(*,*)

      iiiwghts = iii00000
      call ckmem(iiiwghts,iiicrtop)
      iiifreqs = iiiwghts+numpoints*iintfp
      call ckmem(iiifreqs,iiicrtop)
      iiigrdvl = iiifreqs+numpoints*iintfp
      call ckmem(iiigrdvl,iiicrtop)

      call getgridinfo(numpoints,icore(iiiwghts),icore(iiifreqs),
     & icore(iiigrdvl),omega0)
      call prngridinfo(numpoints,icore(iiiwghts),icore(iiifreqs),
     & icore(iiigrdvl))
      
      iiipolra = iiigrdvl+numpoints*iintfp
      call ckmem(iiipolra,iiicrtop)
      iiipolrb = iiipolra+numpoints*numocca*iintfp
      call ckmem(iiipolrb,iiicrtop)

      call getpolarinfo(numpoints,numocca,icore(iiipolra),
     & numoccb,icore(iiipolrb),system)
      call prnmat('Freq Dep Polar A    ',icore(iiipolra),numpoints,
     & numocca,.true.)
      call prnmat('Freq Dep Polar B    ',icore(iiipolrb),numpoints,
     & numoccb,.true.)
      
      iiic6cff = iiipolrb+numpoints*numoccb*iintfp
      call ckmem(iiic6cff,iiicrtop)

      call dzerosec(numocca*numoccb,icore(iiic6cff))

      call getc6matrix(numocca,numoccb,numpoints,omega0,
     & icore(iiiwghts),icore(iiigrdvl),
     & icore(iiipolra),icore(iiipolrb),
     & icore(iiic6cff),totalc6)
      call prnmat('C6 coefficients AXB ',icore(iiic6cff),numocca,
     & numoccb,.true.)

      write(*,*)
      write(*,50) totalc6
      write(*,*)
      
C      call c6partition(numocca,numoccb,icore(iiic6cff),system)
      
      call aces_fin

 10   format('Total number of gridpoints = ',I6)
 20   format('W_0 value for frequencies  = ',F20.15)
 30   format('Total number of occupied A = ',I6)
 40   format('Total number of occupied B = ',I6)
 50   format('Total C6_AB                = ',F20.15)

      end

