      subroutine gridlebedev(type, dim, wgt, crd)
c
c This routine constructs the point coordinates and weights for the
c Lebedev angular grid of a given type. If called with DIM = 0, then it
c returns the diminesion of the grid type TYPE.
c
c In:  type         : type of the grid
c      dim          : number of reducible (or uncontracted) points
c Out: wgt (dim)    : point weights
c      crd (3, dim) : point coordinates
c
c The correspondence between the grid type, the number of solutions of
c the Lebedev equation, and the numbers of irreducible and reducible
c points is
c
c      Grid type   1   2   3   4   5   6
c      Solutions   9  11  17  23  27  29
c    Irreducible   3   4   6   9  11  12
c      reducible  38  50 110 194 266 302
c
c
c The original version of this routine is KS.SRC/INTGRT/INITLEBED.F
c The orginal comments and credits are given below.
c
c Igor Schweigert, Feb 2005
c $Id: gridlebedev.FPP,v 1.1.1.1 2008/05/29 19:35:40 taube Exp $
c
c References:
c    Lebedev, V. I., Zh. vychisl. Mat. mat. Fiz., 16, 2, 293, 1976.
c    Lebedev, V. I., Zh. vychisl. Mat. mat. Fiz., 15, 1, 48, 1975.
c
c This subroutine sets up one of the angular grids (grid) to be a
c Lebedev grid of type.  There are several grids to choose from,
c obtained by solving the Lebedev equations.  Only certain solutions are
c possible.  The solutions included here are: N= 9, 11, 17, 23, 27, 29
c The parameter maxlebgrid is the number of grids implemented.
c
c The grid referred to by each of these solutions can be written in two
c forms: contracted or expanded.  For simplicity, each of the grids are
c included here in their contracted form (for ease of typing them in and
c to eliminate errors), but are then expanded for use in the rest of the
c program.  The contracted form of a grid consists of a number of
c weights and points.  For example, the contracted form for N=9 consists
c of 3 weights and points:
c
c     0.010   1.0  0.0  0.0
c     0.032   1.0  1.0  1.0
c     0.029   0.9  0.5  0.0
c
c The first column is the weight.  The last three columns are
c coordinates.  The expanded form consists of all possible combinations
c of the three coordinates (both positive and negative) in (x,y,z)
c coordinates, each with the given weight.  So, the N=9 solution
c consists of 6 points:
c
c   (0,0,1) (0,0,-1) (0,1,0) (0,-1,0) (1,0,0) (-1,0,0)
c
c all with the weight of 0.010, 8 points:
c
c   (1,1,1) (1,1,-1) (1,-1,1) (-1,1,1)
c   (1,-1,-1) (-1,1,-1) (-1,-1,1) (-1,-1,-1)
c
c all with the weight of 0.032, and 24 points:
c
c   (0.9,0.5,0) (0.9,-0.5,0) (-0.9,0.5,0) (-0.9,-0.5,0)
c   (0.9,0,0.5) (0.9,0,-0.5) (-0.9,0,0.5) (-0.9,0,-0.5)
c   (0,0.9,0.5) (0,0.9,-0.5) (0,-0.9,0.5) (0,-0.9,-0.5)
c   (0.5,0.9,0) (0.5,-0.9,0) (-0.5,0.9,0) (-0.5,-0.9,0)
c   (0.5,0,0.9) (0.5,0,-0.9) (-0.5,0,0.9) (-0.5,0,-0.9)
c   (0,0.5,0.9) (0,0.5,-0.9) (0,-0.5,0.9) (0,-0.5,-0.9)
c
c each with the weight of 0.029.  This gives a total of 6+8+24=38 points.
c
      implicit none
c
c     Arguments
c
      integer
     &     type, dim
c
      double precision
     &     wgt (dim), crd (3, dim)
c
c     Local variables
c     
      integer ntyp, maxipt
      parameter (ntyp   = 6)
      parameter (maxipt = 12)
c
      double precision zero, one, sqrt2,sqrt3
      parameter (zero=0d0)
      parameter (one=1d0)
      parameter (sqrt2=0.707106781186547462d0)
      parameter (sqrt3=0.577350269189625842d0)
c
      integer
     &     nipt (ntyp), nrpt (ntyp), ip, np, iipt, irpt,
     &     i, ix, nx, iy, ny, iz, nz
c     
      double precision
     &     wipt (ntyp, maxipt), cipt (ntyp, 3, maxipt),
     &     x, y, z, xp (6), yp (6), zp (6), sign (2)
c     
      data sign/1d0, -1d0/
c
      logical same, different
      same (x,y) = dabs (x-y) .lt. 1d-14
      different (x,y) = dabs (x-y) .gt. 1d-14
c
c     The actual definition of the Lebedev grids.
c     
c     Legend: nipt - number of irreducible points, nrpt - number of
c     reducible points, wipt - weights, cipt - coordinates of
c     irreducible points.
c     
      data nipt/3,4,6,9,11,12/
      data nrpt/38,50,110,194,266,302/
c
c     Grid #1: N = 9
c
      data (wipt(1,i), cipt(1,1,i), cipt(1,2,i), cipt(1,3,i),i=1,3)/
     &    0.00952380952397d0,one,zero,zero,
     &    0.0321428571429d0,sqrt3,sqrt3,sqrt3,
     &    0.0285714285714d0,0.888073833977d0,0.459700843381d0,zero/
c
c     Grid #2: N = 11
c
      data (wipt(2,i), cipt(2,1,i), cipt(2,2,i), cipt(2,3,i),i=1,4)/
     &    0.0126984126985d0,one,zero,zero,
     &    0.0225749559083d0,sqrt2,sqrt2,zero,
     &    0.02109357d0,sqrt3,sqrt3,sqrt3,
     &    0.0201733355379d0,
     &    0.301511344578d0,0.301511344578d0,0.904534033733d0/
c
c     Grid #3: N = 17, nipt = 6, nrpt = 110
c
      data (wipt(3,i), cipt(3,1,i), cipt(3,2,i), cipt(3,3,i),i=1,6)/
     &    0.003828270494940d0,one,zero,zero,
     &    0.009793737512490d0,sqrt3,sqrt3,sqrt3,
     &    0.008211737283190d0,
     &    0.1851156353450d0,0.1851156353450d0,0.9651240350870d0,
     &    0.009595471336070d0,
     &    0.3956894730560d0,0.3956894730560d0,0.8287699812530d0,
     &    0.00994281489118d0,
     &    0.690421048382d0,0.690421048382d0,0.215957291846d0,
     &    0.00969499636166d0,
     &    0.878158910604d0,0.478369028812d0,zero/
c
c     Grid #4: N = 23, nipt = 9, nrpt = 194
c
      data (wipt(4,i), cipt(4,1,i), cipt(4,2,i), cipt(4,3,i),i=1,9)/
     &    0.00178234044724d0,one,zero,zero,
     &    0.00571690594998d0,sqrt2,sqrt2,zero,
     &    0.00557338317884d0,sqrt3,sqrt3,sqrt3,
     &    0.00551877146727d0,
     &    0.444693317871d0,0.444693317871d0,0.777493219315d0,
     &    0.005158237711810d0,
     &    0.2892465627580d0,0.2892465627580d0,0.9125090968670d0,
     &    0.00560870408259d0,
     &    0.671297344270d0,0.671297344270d0,0.314196994183d0,
     &    0.00410677702817d0,
     &    0.129933544765d0,0.129933544765d0,0.982972302707d0,
     &    0.00505184606462d0,
     &    0.938319218138d0,0.345770219761d0,zero,
     &    0.00553024891623d0,
     &    0.836036015482d0,0.159041710538d0,0.525118572443d0/
c
c     Grid #5: N = 27, nipt = 11, nrpt = 266
c
      data (wipt(5,i), cipt(5,1,i), cipt(5,2,i), cipt(5,3,i),i=1,11)/
     &    -0.00131376912733d0,one,zero,zero,
     &    -0.00252272870489d0,sqrt2,sqrt2,zero,
     &    0.00418685388170d0,sqrt3,sqrt3,sqrt3,
     &    0.00531516797782d0,
     &    0.703937339159d0,0.703937339159d0,0.0945750764036d0,
     &    0.00425613135143d0,
     &    0.662033866370d0,0.662033866370d0,0.351315128565d0,
     &    0.00411248239441d0,
     &    0.464744872642d0,0.464744872642d0,0.753673939251d0,
     &    0.00359558489976d0,
     &    0.327742065497d0,0.327742065497d0,0.886098344997d0,
c The following point has been fixed.  The z value was originally
c 0.98 6948074629  (the space added for clarification) changed to
c 0.9896948074629
     &    0.00404714237709d0,
     &    0.101252624857d0,0.101252624857d0,0.9896948074629d0,
     &    0.00422958270065d0,
     &    0.850650808352d0,0.525731112119d0,zero,
     &    0.00407146759383d0,
     &    0.819343388819d0,0.524493924092d0,0.231479015871d0,
c The following point is suspect!  (X**2 + Y**2 + Z**2)=0.99999999950000
c which is off from 1.0 by a lot further then the other points.  Every
c other point is off by less then 1.0e-12.
     &    0.00408091422578d0,
     &    0.939227929750d0,0.323348454269d0,0.115311196765d0/
c
c     Grid #6: N = 29, nipt = 12, nrpt = 302
c
      data (wipt(6,i), cipt(6,1,i), cipt(6,2,i), cipt(6,3,i),i=1,12)/
     &    0.000854591172878d0,one,zero,zero,
     &    0.00359911928502d0,sqrt3,sqrt3,sqrt3,
     &    0.00365004580768d0,
     &    0.701176641609d0,0.701176641609d0,0.129238672710d0,
     &    0.00360482260142d0,
     &    0.656632941022d0,0.656632941022d0,0.371034178385d0,
     &    0.00357672966173d0,
     &    0.472905413258d0,0.472905413258d0,0.743452042987d0,
     &    0.00344978842429d0,
     &    0.351564034558d0,0.351564034558d0,0.867643624544d0,
     &    0.00310895312238d0,
     &    0.221964523631d0,0.221964523631d0,0.949454317226d0,
     &    0.00235210141366d0,
     &    0.0961830852303d0,0.0961830852303d0,0.990705621379d0,
     &    0.00360082093222d0,
     &    0.820326419828d0,0.571895589188d0,zero,
     &    0.00298234496317d0,
     &    0.964408914879d0,0.264415288706d0,zero,
     &    0.00357154055427d0,
     &    0.251003475177d0,0.800072749407d0,0.544867737258d0,
     &    0.00339231220501d0,
     &    0.902442529533d0,0.412772408317d0,0.123354853258d0/
c
c     Messages to the user.
c
 5905 format (
     &     t3,'@GRIDLEBEDEV-F, Invalid type of Lebedev grid: ',i3,'.',/,
     &     t3,'Valid types are 1-',i2,/)
 5910 format (
     &     t3,'@GRIDLEBEDEV-F: Arrays are too small. Type ',i2,' grid',/
     &     t3,'has ',i3,' reducible points.',/)
 5915 format (
     &     t3,'@GRIDLEBEDEV-F: The number of the reducible points does'/
     &     t3,'not match the definition of the grid type ',i2,'.',/)
c
c     Check if a valid grid is requested.
c
      if (type.gt.ntyp) then
         write (6, 5905) type, ntyp
         call errex
      endif
c
c     If DIM = 0, assign DIM the dimension of the requested grid and return.
c
      if (dim .eq. 0) then
         dim = nrpt (type)
         return
      endif
c
c     Check if the size of the arrays is large enough to keep the weight
c     and coordinates of the reducible set of points.
c
      if (dim .lt. nrpt (type)) then
         write (6, 5910) dim, type, nrpt (type)
         call errex
      endif
c
c     Calculate the reducible set of points: <> Figure out all
c     nontrivial combinations of the irreducible coordinates.
c
      irpt = 0
      do iipt = 1, nipt (type)
c
         x = cipt (type, 1, iipt)
         y = cipt (type, 2, iipt)
         z = cipt (type, 3, iipt)
c     
         np = 1
         xp (np) = x
         yp (np) = y
         zp (np) = z
c
         if (different (x,y)) then
            np = np + 1
            xp (np) = y
            yp (np) = x
            zp (np) = z
         endif
c
         if (different (x,z)) then
            np = np + 1
            xp (np) = z
            yp (np) = y
            zp (np) = x
         endif
c
         if (different (y,z)) then
            np = np + 1
            xp (np) = x
            yp (np) = z
            zp (np) = y
         endif
c
         if (different(x,y).and.different(x,z).and.different(y,z)) then
            np = np + 1
            xp (np) = z
            yp (np) = x
            zp (np) = y
            np = np + 1
            xp (np) = y
            yp (np) = z
            zp (np) = x
         endif
c
         do ip = 1, np
            nx = 2
            ny = 2
            nz = 2
            if (same (xp(ip),0d0)) nx = 1
            if (same (yp(ip),0d0)) ny = 1
            if (same (zp(ip),0d0)) nz = 1
            do ix = 1, nx
               do iy = 1, ny
                  do iz = 1, nz
                     irpt = irpt + 1
                     wgt (irpt) = wipt (type, iipt)
                     crd (1, irpt) = xp (ip) * sign (ix)
                     crd (2, irpt) = yp (ip) * sign (iy)
                     crd (3, irpt) = zp (ip) * sign (iz)
                  enddo
               enddo
            enddo
         enddo
      enddo
c
c     Check the total number of the reducible points.
c
      if (irpt .ne. nrpt (type)) then
         write (6, 5915) type
         call errex
      endif
c     
      return
      end
