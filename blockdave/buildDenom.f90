subroutine buildDenom(Denom,NSize,nbas,nocc,nvirt)
implicit integer(A-Z)
integer,intent(in)::NSize,nbas,nocc,nvirt
double precision,intent(inout)::Denom(NSize,NSize)

double precision::evals(nbas)

CALL GETREC(20,'JOBARC','SCFEVALA',NBAS,evals)

count_i=1
do i=1,nocc
  do a=1,nvirt
    count_j=1
    do j=1,nocc
      do b=1,nvirt
      Denom(count_i,count_j)=1.0d0/(evals(i)+evals(j)-evals(nocc+a)-evals(nocc+b))
      print*,'val:',Denom(count_i,count_j)
        count_j=count_j+1
      enddo
     enddo
     count_i=count_i+1
  enddo
enddo
print*,'final it:',count_i,count_j

end subroutine
  

