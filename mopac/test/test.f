      program testy
      implicit double precision (a-h,o-z)
      parameter (natm = 11)
      integer anum(natm)
      dimension coords(3*natm)
      character*3 ham(2)

      data bohr / 0.52917706d0 /
      data anum / 16, 6, 6, 1, 1, 1, 8, 1, 1, 1, 1 /

      data coords /
     .   0.075890,    -0.078630,     0.015551, 
     .   0.072474,    -0.056662,     1.820287, 
     .   1.832392,     0.027081,    -0.240711, 
     .   2.296964,     0.857763,     0.287255, 
     .   2.370876,    -0.925985,     0.174272, 
     .   2.050611,     0.078347,    -1.303282, 
     .   2.834439,    -2.149219,     0.606271, 
     .   2.198956,    -2.689529,     0.104279, 
     .   0.463732,     0.884955,     2.199225, 
     .  -0.958286,    -0.164095,     2.147378, 
     .   0.658447,    -0.883334,     2.213998 /

      ham(1) = 'AM1'
      ham(2) = 'UHF'


      do 100 i=1, natm*3
         coords(i) = coords(i)/bohr
 100  continue


      call sehess1(natm, anum, coords, ham)

      end

