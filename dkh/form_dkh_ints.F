      Subroutine form_dkh_ints(Dkh_block,Dkh_int,Ishell,Jshell,
     +                         Icntrfns,Jcntrfns,Ncnfns,
     +                         Prim_off,Cont_off,Nshells)

      Implicit Double Precision (A-H, O-Z)

      Dimension Dkh_block(Icntrfns,Jcntrfns)
      Dimension Dkh_int(Ncnfns,Ncnfns)
      Integer Prim_off(Nshells),Cont_off(Nshells)

      ioff = Cont_off(Ishell) - 1
      Joff = Cont_off(Jshell) - 1

#ifdef _DEBUG_LVLM
      write(6,*)
      write(6,"(a,2(1x,i4))")"Ishell,Jshell",Ishell,Jshell
      write(6,"(a,2(1x,i4))")"Ioff,Joff",ioff, joff
      write(6,"(a,2(1x,i4))")"Icntrfns,Jcntrfns",Icntrfns,Jcntrfns
      write(6,"(6(1x,E12.5))") ((Dkh_block(i,j),i=1,Icntrfns),j=1,
     +                          Jcntrfns)
      Write(6,"(a,2i4)") "The DKH integral accumulation:before" 
      write(6,"(6(1x,E12.5))") ((Dkh_int(i,j),i=1,Ncnfns),
     +                           j=1,Ncnfns)
#endif 

      Jct = 0
      Do J = Joff+1, Jcntrfns + Joff
         Jct = Jct + 1
         Ict = 0
         Do I= Ioff+1,Icntrfns + Ioff
            Ict = Ict + 1
            Dkh_int(I,J) = Dkh_block(Ict,Jct)
         Enddo
      Enddo

#ifdef _DEBUG_LVLM
      Write(6,"(a,2i4)") "The DKH integral accumulation:after" 
      write(6,"(6(1x,E12.5))") ((Dkh_int(i,j),i=1,Ncnfns),
     +                           j=1,Ncnfns)
#endif 
      
      Return
      End 
     
