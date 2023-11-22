      subroutine Put_ConvergInfo(OOmicroItr,gk_RMS,gk_MAX_elem,&
     &                          kappa_RMS,kappa_MAX_elem)
      integer,intent(in)::OOmicroItr
      double precision,intent(in)::gk_RMS,gk_MAX_elem
      double precision,intent(in)::kappa_RMS,kappa_MAX_elem

      logical::file_exists

      INQUIRE(FILE="ConvergInfo.txt", EXIST=file_exists)
      if (file_exists) then
        open(unit=1234,file="ConvergInfo.txt",status="old")
        write(1234,*)gk_RMS,gk_MAX_elem,kappa_RMS,kappa_MAX_elem
        close(1234)


      else
        open(unit=1234,file="ConvergInfo.txt",status="new")
        write(1234,*) "RMS MOgrad | MAX MOgrad | RMS Kappa | MAX Kappa"
        write(1234,*)gk_RMS,gk_MAX_elem,kappa_RMS,kappa_MAX_elem
        close(1234)

      endif
      end subroutine
