
c This routine initializes the info common block.

      subroutine aces_com_info_c4
      implicit none
      integer       nocco(2),nvrto(2)
      common /info/ nocco,   nvrto
      save   /info/
      call getrec_c4(-1,'JOBARC','NOCCORB',2,nocco)
      call getrec_c4(-1,'JOBARC','NVRTORB',2,nvrto)
      return
      end

