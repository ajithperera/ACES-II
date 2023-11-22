
      subroutine rmfiles3
      integer icore
      call aces_ja_init
      call aces_io_init(icore,1,0,.false.)
      call aces_io_remove(52,'MOABCD')
      call aces_io_remove(53,'DERINT')
      call aces_io_fin
      call aces_ja_fin
      return
      end

