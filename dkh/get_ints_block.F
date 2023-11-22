      Subroutine get_ints_block(Coreham,Shell_block,Ishell,Jshell,
     +                           Nprims,Prim_off,Cont_off,
     +                           Iprimfns,jprimfns,Nshells)

      Implicit Double Precision (A-H, O-Z)

      Dimension Coreham(Nprims, Nprims)
      Dimension Shell_block(iprimfns,jprimfns)
      Integer Prim_off(Nshells),Cont_off(Nshells)

      ioff = Prim_off(Ishell) - 1
      Joff = Prim_off(Jshell) - 1

#ifdef _DEBUG_LVLM
      write(6,*)
      write(6,"(a,2(1x,i4))")"Ishell,Jshell",Ishell,Jshell
      write(6,"(a,2(1x,i4))")"Ioff,Joff",ioff, joff
      write(6,"(a,2(1x,i4))")"Iprimfns,jprimfns",Iprimfns,jprimfns
      Write(6,*)
#endif 
      Do J = 1, Jprimfns
         Do I= 1, Iprimfns
            shell_block(I,J) = Coreham(Ioff+I,Joff+J)
         Enddo 
      Enddo 

#ifdef _DEBUG_LVLM
      Write(6,"(a,2i4)") "The shell-block",ishell,jshell
      write(6,"(6(1x,E12.5))") ((shell_block(i,j),i=1,Iprimfns),
     +                           j=1,jprimfns)
#endif 
      
      Return
      End 
     
