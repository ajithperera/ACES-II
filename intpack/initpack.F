      subroutine initpack

c This program initializes parameters from the INTPACK namelist.

#include <aces.h>
#include <io.par>
#include <const.par>
#include <intgrtflags.com>

      integer err

      character*20 list(10)
      logical printdef

      printdef=.true.

      call callstack_push('INITPACK')
      call nl_init('INTPACK',err,printdef)

c Namelist entries:
c --------------------------------------------------------------------------
c PRINT_NL   : whether to print the values read in the namelist
c              Default: true
c --------------------------------------------------------------------------
c Integration flags
c
c NUMRADPTS  : the number of radial points to use in the integration around
c              each atomic center (default: 50)
c RADTYP     : the type of radial integration available are:
c                Handy           : (default)
c                Gauss-Legendre  :
c PARTPOLY   : partition polyhedra
c                equal     : equally regardless of atomic center
c                bsrad     : by Bragg-Slater radii (default)
c                dynamic   : dynamic (minimum in density) sizing
c RADSCAL    : how to do radial scaling
c                none      : no radial scaling
c                slater    : use Slater's rules for radial scaling (default)
c PARTTYP    : type of partitioning
c                rigid     : rigid partitioning
c                fuzzy     : fuzzy partitioning (default)
c FUZZYITER  : the number of iterations to use in fuzzy partitioning
c              (default: 4)
c RADLIMIT   : the length of the radial integration (default 3.00)
c
c ANGLGRID    : the angular grid to use by default at each radial point
c               (default: 4)
c GRIDTYPE    : the angular grid type to use by default at each radial point
c                lebedev    : default
c OVERIND     : the number of basis functions in each integral
c                default: 3
c --------------------------------------------------------------------------

      call nl_int('numradpts',50,int_numradpts)

      list(1)='Handy'
      list(2)='Gauss-Legendre'
      call nl_ele('radtyp',list,2,int_radtyp_handy,int_radtyp)

      list(1)='equal'
      list(2)='bsrad'
      list(3)='dynamic'
      call nl_ele('partpoly',list,3,int_partpoly_bsrad,int_partpoly)

      list(1)='none'
      list(2)='Slater'
      call nl_ele('radscal',list,2,int_radscal_slater,int_radscal)

      list(1)='rigid'
      list(2)='fuzzy'
      call nl_ele('parttyp',list,2,int_parttyp_fuzzy,int_parttyp)

      call nl_int('fuzzyiter',4,int_fuzzyiter)

      call nl_real('radlimit',three,int_radlimit)

      call nl_int('anglgrid',4,int_defanglgrid)

      call nl_int('overind',3, int_overlp)

      list(1)='lebedev'
      call nl_ele('gridtype',list,1,int_gridtype_leb,int_defgridtype)

      call nl_term
      call callstack_pop
      return
      end
