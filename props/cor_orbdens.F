      Subroutine Cor_orbdens(Evecs,Docc,Dens,Norbs,Nbfns,Naobfns,Ispin)

      Implicit Double Precision(A-H,O-Z)

      Dimension Evecs(Naobfns,Nbfns), Dens(Naobfns*Naobfns,Norbs),
     +          Docc(Nbfns,Nbfns)

      Data Done,Dnull /1.0D0,0.0D0/

      Call Dzero(Dens,Naobfns*Naobfns*Norbs)

      Do Iorbs = 1, Norbs 
         Call Mk_corden(Evecs(1,Iorbs),Dens(1,Iorbs),Nbfns,Naobfns,
     +                  Docc(Iorbs,Iorbs))
#ifdef _DEBUG_LVLM
        Write(6,*) 
        If (Ispin .EQ. 1) Then
        write(6,"(a)") "The alpha orbital density matrices"
        call output(Dens(1,Iorbs),1,Naobfns,1,Naobfns,Naobfns,Naobfns,1)
        else 
        write(6,"(a)") "The beta orbital density matrices"
        call output(Dens(1,Iorbs),1,Naobfns,1,Naobfns,Naobfns,Naobfns,1)
        endif 
#endif

      Enddo 

      Return
      End 


 
