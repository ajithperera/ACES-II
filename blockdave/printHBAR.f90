subroutine printHBAR(Work,NSize, Maxcor,Iuhf,Irrepx)
  integer,intent(in)::Nsize,Maxcor,Iuhf,Irrepx
  double precision, intent(inout)::Work(Nsize)
  integer::k,z


 open(1000,file='HBARout')
 print*,'writing EOM HBAR'
 do k=1,Nsize
  Work=0.0d0
  Work(k)=1.0d0 
  Call LANCZOS_DUMP_VEC(Irrepx,Work,NSize,&
                      490,0,0,443,0,IUHF,.False.)
  Call Hbarxc(Work,Maxcor,Iuhf,1,Irrepx)
  Call Loadvec1(Irrepx,Work,Maxcor,Iuhf,490,2,460,&
                           Nsize,.False.)
  
  do z=1,NSize
    write(1000,*)Work(z)
  enddo 
 enddo
close(1000)
end subroutine
