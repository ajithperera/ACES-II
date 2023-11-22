
      subroutine namelist_format
      integer nml0, nml1
      namelist /struct/ nml0, nml1, nml2, nml3, nml4

      print *, "==> namelist formatting test"
      print *, ""
      nml0=0
      nml1=1
      nml2=2
      nml3=3
      nml4=4
      write(*,nml=struct)
      print *, ""

      return
      end

