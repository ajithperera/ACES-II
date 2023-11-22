      program main   
      implicit none
      Integer iuhf, i0, ICoreDim

#include "icore.com"

      call aces_init(iCore,i0,iCoreDim,iUHF,.true .)
      call lcc(iCore(i0),iCoreDim,iUHF)
      call aces_fin

      end

