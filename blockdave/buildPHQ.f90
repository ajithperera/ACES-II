subroutine buildPHQ(V,cisEvecs,Nsize,Nblocks,PHQ)
  integer,intent(in)::NSize,Nblocks
  double precision,intent(in)::V(Nsize,Nsize),cisEvecs(Nsize,Nblocks)
  double precision,intent(inout)::PHQ(1,Nblocks-1)


! Assumes that we are studying the lowest root
! Computes PHQ == QHP
!  PHQ=
!print*,matmul(matmul(transpose(cisEvecs(:,1)),V),cisEvecs(:,2:Nblocks))
PHQ=matmul(matmul(transpose(cisEvecs(:,1:1)),V),cisEvecs(:,2:Nblocks))
!print*,matmul(cisEvecs(:,1),V)
end subroutine
