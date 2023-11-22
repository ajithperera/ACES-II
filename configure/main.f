
      program test
      integer iIntLn
      integer icore(2)
      common // icore

      print '()'
      print *, 'CONFIGURATION TESTS FOR PORTING ACES II'
      print '()'

c   o build-specific tests
      call int_size(iIntLn)

c   o pre-processor defines tests
      call recl_units(iIntLn)
      call mem_res(iIntLn)
      call fortran_names
      call args

c   o informational tests
      call c_sizeof
      call c_endian
      call c_heap(icore)
      call chars
      call shifts
      call namelist_format

      end

